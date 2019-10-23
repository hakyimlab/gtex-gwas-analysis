library(dplyr)
library(ggplot2)
theme_set(theme_bw(base_size=18))
source('../helpers/rlib_doc.R')
library(stringr)
regulations = c('expression-fixed', 'splicing')
traits = read.csv('data/misc/Ayellet-Exclude-List_UKBB_GWAS_hg38_v8_Harmon_Imputed_deflation.csv') %>% mutate(trait = get_trait(File.name)) %>% filter(Deflated.imputed.GWAS.p.value.distribution..1.strong.deflation..2.mild.deflation..0.no.deflation. == 0)
tissues = readRDS('data/misc/tissue_abbr_sample_size.rds')

df_cor = data.frame()
for(r in regulations) {
  sub = read.table(paste0('data/mediation_analysis_new_setting/', r, '/cor_gwas_qtl.txt.gz'), stringsAsFactors = F, header = T) %>% filter(!is.na(cor), trait %in% traits$trait) %>% mutate(regulation = r)
  df_cor = rbind(df_cor, sub)
}
df_cor = df_cor %>% mutate(cor_se = (cor_ci95upper - cor_ci95lower) / 2 / 1.96)

temp = df_cor %>% filter(type != 'negative') %>% group_by(regulation, trait, type) %>% summarize(median_across_tissue = median(cor))

library(ggridges)
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/8151c6fe70e3d4ee43d9ce340ecc0eb65172e616/my_ggplot_theme.R')
temp$regulation[temp$regulation == 'expression-fixed'] = 'expression'
temp$type[temp$type == 'positive'] = 'observed'
temp$type[temp$type == 'shuffled_within_ldscore_bin'] = 'shuffled'
temp = temp %>% mutate(merge = paste(regulation, type))
temp$merge = factor(temp$merge, levels = c('expression observed', 'expression shuffled', 'splicing observed', 'splicing shuffled')[4:1])
mymixer = c(
  'splicing observed' = rgb(200, 90, 40, maxColorValue = 255),
  'expression observed' = rgb(42, 155, 204, maxColorValue = 255),
  'splicing shuffled' = 'gray',
  'expression shuffled' = 'gray'
)
p = temp %>% ggplot() + geom_density_ridges(aes(y = merge, fill = merge, color = merge, x = median_across_tissue, height=..density..), alpha = .6) # + geom_boxplot(aes(x = regulation, y = median_across_tissue, fill = type_white), width = 0.05, position = dodge) 
p = p + scale_fill_manual(values = mymixer) + scale_color_manual(values = mymixer) + theme(legend.position='none') + xlab('Correlation between QTL and GWAS effects')
p = p + th
p = p + geom_text(data = data.frame(
  x = c(0.2, 0.2), y = c(1.5, 3.5), text = c('splicing', 'expression')
), aes(x = x, y = y, label = text), size = 9
)
p = p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

ggsave('../../output/SFIG-COR-DELTA-GAMMA.png', p, width = 7, height = 5.5)

p = temp %>% ggplot() + geom_density_ridges(aes(y = merge, fill = merge, color = merge, x = median_across_tissue, height=..density..), alpha = .6, jittered_points = T) # + geom_boxplot(aes(x = regulation, y = median_across_tissue, fill = type_white), width = 0.05, position = dodge) 
p = p + scale_fill_manual(values = mymixer) + scale_color_manual(values = mymixer) + theme(legend.position='none') + xlab('Correlation between QTL and GWAS effects')
p = p + th
p = p + geom_text(data = data.frame(
  x = c(0.2, 0.2), y = c(1.5, 3.5), text = c('splicing', 'expression')
), aes(x = x, y = y, label = text), size = 9
)
p = p + theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())

ggsave('../../output/SFIG-COR-DELTA-GAMMA-with-jitter.png', p, width = 7, height = 5.5)







# sigma plot
df_sigma = data.frame()
df_stat = data.frame()
models = 'model3'  # c('model3', 'model4', 'model8')
for(r in regulations) {
  for(m in models) {
    sub = read.table(paste0('data/mediation_analysis_new_setting/', r, '/sigma_gene-', m, '.RE.txt.gz'), stringsAsFactors = F, header = T) %>% mutate(regulation = r, model = m) %>% filter(trait %in% traits$trait)
    sub2 = read.table(paste0('data/mediation_analysis_new_setting/', r, '/sigma_gene-', m, '.stat.txt.gz'), stringsAsFactors = F, header = T) %>% filter(trait %in% traits$trait) %>% mutate(regulation = r, model = m)
    # sub = inner_join(sub, sub2, by = c('trait', 'tissue'))
    df_sigma = rbind(df_sigma, sub)
    df_stat = rbind(df_stat, sub2)
  }
}


temp2 = df_sigma %>% filter(grp == 'gene') %>% group_by(trait, type, regulation, model) %>% summarize(median_across_tissues = median(vcov)) %>% ungroup()

# temp2 %>% ggplot() + geom_boxplot(aes(x = merge, y = median_across_tissues))
temp2$regulation[temp2$regulation == 'expression-fixed'] = 'expression'
temp2$type[temp2$type == 'positive'] = 'observed'
temp2$type[temp2$type == 'negative'] = 'shuffled'
temp2 = temp2 %>% mutate(merge = paste(regulation, type))
temp2$merge = factor(temp2$merge, levels = c('expression observed', 'expression shuffled', 'splicing observed', 'splicing shuffled')[4:1])

p = temp2 %>% ggplot() + geom_density_ridges(aes(y = merge, fill = merge, color = merge, x = median_across_tissues, height=..density..), alpha = .6) # + geom_boxplot(aes(x = regulation, y = median_across_tissue, fill = type_white), width = 0.05, position = dodge) 
p = p + th
p = p + theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
p = p + scale_fill_manual(values = mymixer) + scale_color_manual(values = mymixer) + theme(legend.position='none') + xlab(expression(paste('Estimated ', sigma[gene]^2)));p
p = p + scale_x_sqrt(expand = c(0, 0))
ggsave('../../output/FIG-DOSE-RESPONSE-CONCORDANCE-B.png', p, width = 5.7, height = 5.5)


p = temp2 %>% ggplot() + geom_density_ridges(aes(y = merge, fill = merge, color = merge, x = median_across_tissues, height=..density..), alpha = .6, jittered_points = T) # + geom_boxplot(aes(x = regulation, y = median_across_tissue, fill = type_white), width = 0.05, position = dodge) 
p = p + th

p = p + theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
p = p + scale_fill_manual(values = mymixer) + scale_color_manual(values = mymixer) + theme(legend.position='none') + xlab(expression(paste('Estimated ', sigma[gene]^2)));p
p = p + scale_x_sqrt(expand = c(0, 0))
ggsave('../../output/FIG-DOSE-RESPONSE-CONCORDANCE-B-with-jitter.png', p, width = 5.7, height = 5.5)
