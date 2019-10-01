suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(argparse)
  library(glue)
  library(grid)
  library(gridExtra)
  library(lattice)
})

tissue_lst <- readLines("../data/tissue_list.txt")
phenotype_lst <- readLines("../data/phenotype_list.txt")
metadata <- read_tsv("../data/gwas_metadata.txt") %>% rename(phenotype=Tag)
df_master <- read_tsv("../data/eur_pvs_analysis/eur_genes_tissue_spec_rcp0.1.txt")

for (tissue in tissue_lst)
{
  #short_phenotype <- as.character(metadata[metadata$phenotype==phenotype, "new_Phenotype"])
  df <- df_master[df_master$tissue==tissue,]
  #tmp <- df[names(df) %in% c("zscore_rank1_beta","zscore_rank2_beta")]
  chisq <- chisq.test(table(df[names(df) %in% c("zscore_share_group", "zscore_concordance")]))
  plot_title <- "{gsub(pattern='_', replacement=' ', x=tissue)} - All Phenotypes - Enloc - EUR\nPrimary vs. secondary eQTL by Z-score" %>% glue::glue()
  
  pp <- ggplot(df, aes(x=zscore_rank1_beta, y=zscore_rank2_beta))
  pp <- pp + geom_abline(slope=1,intercept=0,color="darkgrey",size=1) + geom_hline(yintercept=0,size=1,color="darkgrey") + geom_vline(xintercept=0,size=1,color="darkgrey")
  pp <- pp + theme_bw() + theme(axis.text=element_text(size=20), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30)) 
  pp <- pp + ggtitle(plot_title)
  pp <- pp + xlab(expression(beta[prim])) + ylab(expression(beta[sec])) + labs(color="Tissue Sharing")
  pp <- pp + geom_point(alpha=.6, size=6, color="black")
  pp <- pp + coord_fixed(ratio=1, xlim=c(-5,5), ylim=c(-5,5))
  
#  pvalue <- textGrob(paste("p-value: ",round(chisq$p.value,7)))
#  pp2 <- ggplot(data = as.data.frame(chisq$residuals), aes(x=zscore_concordance, y=zscore_share_group, fill=Freq))
#  pp2 <- pp2 + geom_tile(color = "black")
#  pp2 <- pp2 + scale_fill_gradient2(low = "blue", high = "red", mid = "white", limits=c(-3,3), name="Pearson\nResidual") 
#  pp2 <- pp2 + theme_bw() + theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  
#  lay <- rbind(c(1,1,1,1,NA,NA),
#               c(1,1,1,1,2,2),
#               c(1,1,1,1,2,2),
#               c(1,1,1,1,3,NA))
  
  filename <- glue::glue("../docs/figure/prim_vs_sec_eur.Rmd/FIG-DOSE-RESPONSE-CONCORDANCE-C.png")
  print(filename)
  png(filename, width = 1500, height = 1500, res = 150)
# grid.arrange(grobs = list(pp, pp2, pvalue), layout_matrix = lay)
  print(pp)
  dev.off()
}
