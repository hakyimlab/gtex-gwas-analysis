library(tidyverse)
library(yaml)


dapg_file_pattern <- "../data/dapg_variants/with_mashr_results/DAPG_with_mashr__{tissue}.rds"

tissues <- readLines("../data/tissue_list.txt")
phenotypes <- readLines("../data/phenotypes.txt")

rows <- vector(mode="list", length=0)
for (tissue in tissues) {
  dfs <- readRDS(glue::glue(dapg_file_pattern))
  for (phenotype in phenotypes) {
    df <- dfs[[phenotype]]
    corr <- cor(abs(df$eqtl_effect_size), abs(df$gwas_effect_size), use = "complete.obs")
    corr_shuffled <- cor(abs(df$eqtl_effect_size), sample(abs(df$gwas_effect_size)), use = "complete.obs")
    corr_sp <- cor(abs(df$eqtl_effect_size), abs(df$gwas_effect_size), use = "complete.obs", method="spearman")
    corr_shuffled_sp <- cor(abs(df$eqtl_effect_size), sample(abs(df$gwas_effect_size)), use = "complete.obs", method="spearman")
    rows <- c(
      rows,
      list(data.frame(
        "tissue"=c(tissue,tissue), 
        "phenotype"=c(phenotype,phenotype), 
        pearson=c(corr, corr_shuffled), 
        spearman=c(corr_sp, corr_shuffled_sp),
        data=c("real", "shuffled"))
      )
    )
  }
}

suppressWarnings({
  corr_df <- bind_rows(rows)
})

saveRDS(corr_df, "../data/eqtl_gwas_correlation.rds")
