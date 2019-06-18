library(tidyverse)
library(yaml)

dapg_file_pattern <- "../data/dapg_selected_variants/splicing/sqtl/{tissue}.v8.sqtl_filtered_DAPG.txt" # DAPG_with_mashr__{tissue}.rds"
gwas_file_pattern <- "../data/dapg_selected_variants/splicing/gwas/{phenotype}__filtered_by_dapg_sqtl.txt"

tissues <- readLines("../data/tissue_list.txt")
phenotypes <- readLines("../data/phenotype_list.txt")

rows <- vector(mode="list", length=0)
for (phenotype in phenotypes) {
  message(glue::glue("Processing {phenotype}..."))
  gwas <- read.table(glue::glue(gwas_file_pattern), header=TRUE) %>% 
    select(panel_variant_id, effect_size, standard_error) %>%
    rename(variant_id=panel_variant_id, gwas_effect_size=effect_size, gwas_se=standard_error) 
  
  for (tissue in tissues) {
    sqtl <- read.table(glue::glue(dapg_file_pattern), header=TRUE) %>% 
            select(gene_id, variant_id, slope) %>%
            rename(eqtl_effect_size=slope)
    
    df <- inner_join(sqtl, gwas, by="variant_id")
    
    corr <- cor(
      abs(df$eqtl_effect_size), 
      abs(df$gwas_effect_size), 
      use = "complete.obs")
    
    corr_shuffled <- cor(
      abs(df$eqtl_effect_size), 
      sample(abs(df$gwas_effect_size)), 
      use = "complete.obs")

    corr_sp <- cor(
      abs(df$eqtl_effect_size), 
      abs(df$gwas_effect_size), 
      use = "complete.obs", method="spearman")
    
    corr_shuffled_sp <- cor(
      abs(df$eqtl_effect_size), 
      sample(abs(df$gwas_effect_size)), 
      use = "complete.obs", method="spearman")
    
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

saveRDS(corr_df, "../output/sqtl_gwas_correlation.rds")
