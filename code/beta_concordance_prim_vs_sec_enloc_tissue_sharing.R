suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(argparse)
  library(glue)
})

dapg_file_pattern <- "../data/dapg_selected_variants/expression/gwas_and_eqtl/DAPG_with_mashr__{tissue}.rds"
tissue_lst <- readLines("../data/tissue_list.txt")
metadata <- read_tsv("../data/gwas_metadata.txt") %>% rename(phenotype=Tag)
df_master <- read_tsv("../data/coloc_enloc_gene_tissue_specificity.txt")

for (tissue in tissue_lst)
{
  df <- df_master[df_master$tissue==tissue,]
  
  plot_title <- "{gsub(pattern='_', replacement=' ', x=tissue)} - All Phenotypes - Enloc\nPrimary vs. secondary eQTL by Z-score" %>% glue::glue()
  
  pp <- ggplot(df, aes(x=zscore_rank1_beta, y=zscore_rank2_beta, color=zscore_share_group)) #color parameter to show genes that are tissue sharing/specific
  pp <- pp + geom_point(alpha=.6, size=4) #+ geom_encircle(aes(group=cluster), s_shape = 1, expand = 0, show.legend = FALSE)    needs library(ggalt)
  pp <- pp + theme_bw() + theme(axis.text=element_text(size=20), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), plot.title=element_text(size=30))
  pp <- pp + geom_abline(slope=1,intercept=0)
  pp <- pp + stat_density2d(color = alpha(colour="grey",alpha=0.6), contour = TRUE, size=1.25)
  pp <- pp + ggtitle(plot_title)
  pp <- pp + xlab(expression(beta[prim])) + ylab(expression(beta[sec])) + labs(color="Tissue Sharing")
  pp <- pp + coord_fixed(ratio=1, xlim=c(-5,5), ylim=c(-5,5))
  
  filename <- glue::glue("../docs/figure/prim_vs_sec_tissues_zscore.Rmd/prim_vs_sec_enloc_zscore_{tissue}.png")
  print(filename)
  png(filename, width = 1500, height = 1500, res = 150)
  print(pp)
  dev.off()
}