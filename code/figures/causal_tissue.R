library(ggplot2)
library(dplyr)
library(pander)
library(scales)
panderOptions('table.split.table', Inf)
theme_set(theme_bw(base_size=15))  

traits = read.delim2('data/misc/gwas_metadata.txt', stringsAsFactors = F) %>% select(Tag, new_Phenotype, Category, color)                      
traits = traits[order(paste(traits$Category, traits$new_Phenotype)), ]
rownames(traits) = NULL
# 
trait_abbrev = read.csv('data/misc/GWAS_metadata.csv', stringsAsFactors = F)
traits = traits %>% left_join(trait_abbrev %>% select(Tag, new_abbreviation), by = 'Tag')
traits %>% pander

df = data.frame()
for(i in 1 : nrow(traits)) {
  filename = paste0('data/causal_tissue/eQTL-mashr-EUR--', traits$Tag[i], '__enrichment.combined.rds')
  if(file.exists(filename)) {
    combined = readRDS(filename)$test
  } else {
    next
  }
  combined = combined %>% mutate(trait = traits$Tag[i], color = traits$color[i], name = traits$new_abbreviation[i])
  df = rbind(df, as.data.frame(combined))
}

df$merged_membership[df$merged_membership == "Uterus-vargina-ovary"] = "Uterus-vagina-ovary"
df_clean = df %>% group_by(trait) %>% summarize(ntissue = sum(!is.na(odds_ratio))) # %>% filter(ntissue >= 8)
df_clean = df %>% filter(trait %in% df_clean$trait)
df_clean$name = factor(df_clean$name, levels = traits$new_abbreviation[order(paste(traits$Category, traits$new_abbreviation))])
color_guide = df_clean$color[order(df_clean$name)][!duplicated(df_clean$trait)]
names(color_guide) = df_clean$name[order(df_clean$name)][!duplicated(df_clean$trait)]
temp = df_clean %>% mutate(signed_log10p = -log10(pval) * ((odds_ratio > 1) - 0.5) * 2, is_signif = ifelse(pval < 0.05, 'dot', 'no_dot'))
maxr = max(temp$signed_log10p, na.rm = T)
minr = min(temp$signed_log10p, na.rm = T)
p = temp %>% ggplot(aes(x = name, y = merged_membership)) + geom_tile(aes(fill = signed_log10p), color = 'gray') + geom_point(aes(size = is_signif)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_gradientn(colors = c("#0000FF", 'white', "#808000"), na.value = 'white', limits = c(minr, maxr), values =  rescale(c(minr,0,maxr)))+
     scale_size_manual(values=c(dot=1, no_dot=NA), guide="none") + theme(legend.position="bottom")
p = p + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_rect(fill = "lightgray"))

p = p + theme(
  plot.title = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.ticks.x = element_blank(),
  axis.line.x = element_blank()
)
p = p + geom_point(aes(x = name, y = 0.2, color = name), show.legend = F, shape = 15, size = 4) + scale_color_manual(values = as.character(color_guide[as.character(df_clean$name[order(df_clean$name)][!duplicated(df_clean$trait)])])) + coord_cartesian(ylim = c(0.5, 15), expand = T)
ggsave('../../output/causal_tissue_complete-mashr-EUR.png', p, width = 13, height = 6, unit = 'in')
