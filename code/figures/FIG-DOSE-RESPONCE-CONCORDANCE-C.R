suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(argparse)
  library(glue)
})

metadata <- read_tsv("../../data/gwas_metadata.txt") %>% rename(phenotype=Tag)
df_master <- read_tsv("../../data/minimal_data/eur_genes_tissue_spec_rcp0.1.txt")

df <- df_master[df_master$tissue=="Whole_Blood",]

pp <- ggplot(df, aes(x=zscore_rank1_beta, y=zscore_rank2_beta))
pp <- pp + geom_abline(slope=1,intercept=0)
pp <- pp + geom_vline(xintercept=0) + geom_hline(yintercept=0) 
pp <- pp + theme_bw() + theme(axis.text=element_text(size=20), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), panel.grid=element_blank())
pp <- pp + xlab(expression(beta[prim])) + ylab(expression(beta[sec])) + labs(color="Tissue Sharing")
pp <- pp + geom_point(alpha=.6, size=6)
pp <- pp + coord_fixed(ratio=1, xlim=c(-3.9,3.9), ylim=c(-3.9,3.9))
  
filename <- glue::glue("../../output/FIG-DOSE-RESPONSE-CONCORDANCE-C.png")
print(filename)
png(filename, width = 1500, height = 1500, res = 150)
print(pp)
dev.off()
