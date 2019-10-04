suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(argparse)
  library(glue)
  library(grid)
  library(gridExtra)
  library(lattice)
})

df_master <- read_tsv("../../data/minimal_data/eur_genes_tissue_spec_rcp0.1.txt")

df <- df_master[df_master$tissue=="Whole_Blood",]
chisq <- chisq.test(table(df[names(df) %in% c("effect_size_share_group", "effect_size_concordance")]))
  
pp <- ggplot(df, aes(x=effect_size_rank1_beta, y=effect_size_rank2_beta, color=effect_size_concordance))
pp <- pp + geom_abline(slope=1,intercept=0) + geom_vline(xintercept=0) + geom_hline(yintercept=0) 
pp <- pp + theme_bw() + theme(axis.text=element_text(size=20), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), legend.position=0) 
pp <- pp + theme(panel.grid=element_blank())
pp <- pp + xlab(expression(beta[prim])) + ylab(expression(beta[sec])) + labs(color="Tissue Sharing")
pp <- pp + geom_point(alpha=.6, size=2)
pp <- pp + coord_fixed(ratio=1, xlim=c(-4,4), ylim=c(-4,4))
  
pvalue <- textGrob(paste("p-value: ",round(chisq$p.value,7)))
pp2 <- ggplot(data = as.data.frame(chisq$residuals), aes(x=effect_size_concordance, y=effect_size_share_group, fill=Freq))
pp2 <- pp2 + geom_tile(color = "black")
pp2 <- pp2 + scale_fill_gradient2(low = "blue", high = "red", mid = "white", limits=c(-3,3), name="Pearson\nResidual") 
pp2 <- pp2 + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  
lay <- rbind(c(1,1,1,1,NA,NA),
             c(1,1,1,1,2,2),
             c(1,1,1,1,2,2),
             c(1,1,1,1,3,NA))
  
filename <- glue::glue("../../output/SFIG-CONCORDANCE-MEDIATING-EFFECTS-RANK-BY-EFFECT-SIZE.png")
print(filename)
png(filename, width = 1500, height = 1000, res = 150)
grid.arrange(grobs = list(pp, pp2, pvalue), layout_matrix = lay)
dev.off()

