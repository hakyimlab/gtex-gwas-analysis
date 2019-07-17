library(tidyverse)
library(yaml)

# Author: Rodrigo Bonazzola
# This script computes correlation (Pearson and Spearman) between GWAS and sQTL effect sizes.

setwd("~/Genomica/GTEx/gtex-gwas-analysis/code/")
type <- "splicing"
dapg_file_pattern <- "../data/dapg_selected_variants/splicing/sqtl/{tissue}.v8.sqtl_filtered_DAPG.txt" # DAPG_with_mashr__{tissue}.rds"
gwas_file_pattern <- "../data/dapg_selected_variants/splicing/gwas/{phenotype}__filtered_by_dapg_sqtl.txt"
tissues <- readLines("../data/tissue_list.txt")
phenotypes <- readLines("../data/phenotype_list.txt")
output_file <- "../output/sqtl_gwas_correlation_with_ci.tsv"

rows <- vector(mode="list", length=0)

# boolean indicating whether to use all the 
use_top_n <- c(FALSE, TRUE)
shuffle <- c(FALSE, TRUE)
worst_tissue <- "Kidney_Cortex"
worst_tissue_df <- read.table(glue::glue(gsub("tissue", "worst_tissue", dapg_file_pattern)), header = TRUE)
n_variants <- nrow(worst_tissue_df)


for (phenotype in phenotypes) {
  message(glue::glue("Processing {phenotype}..."))
  
  gwas <- read.table(glue::glue(gwas_file_pattern), header=TRUE) %>% 
    select(
      panel_variant_id, 
      effect_size, 
      standard_error
    ) %>% 
    dplyr::rename(
      variant_id = panel_variant_id, 
      gwas_effect_size = effect_size, 
      gwas_se = standard_error
    )
  
  for (tissue in tissues) {
    for (use_top_n_ in use_top_n) {
      for (shuffle_ in shuffle) {
        for (dapg_file_pattern_ in dapg_file_pattern) {
        
          # message(glue::glue("  Processing {tissue}..."))
          qtl_df <- read.table(glue::glue(dapg_file_pattern), header=TRUE)
          
          if (use_top_n_)
            qtl_df <- qtl_df %>% top_n(n_variants, wt=-pval_nominal)
      
          qtl_df <- qtl_df %>% 
            select(gene_id, variant_id, slope) %>%
            dplyr::rename(qtl_effect_size=slope)
          
          suppressWarnings({
            df <- inner_join(qtl_df, gwas, by="variant_id")
          })
      
          if (shuffle_) {
            corr_test <- cor.test(
              abs(df$qtl_effect_size), 
              sample(abs(df$gwas_effect_size)), 
              use = "complete.obs"
            )
          } else {
            corr_test <- cor.test(
              abs(df$qtl_effect_size), 
              abs(df$gwas_effect_size), 
              use = "complete.obs"
            )
          }
          
          new_row <- data.frame(
            "tissue" = tissue,
            "phenotype" = phenotype,
            "pearson" = corr_test$estimate,
            "ci95_lower" = corr_test$conf.int[1],
            "ci95_upper" = corr_test$conf.int[2],
            "pvalue" = corr_test$p.value,
            "data" = ifelse(shuffle_, "shuffled", "real"),
            "type" = type,
            "top_n" = nrow(df)
          )
          
          rows <- c(rows, list(new_row))
          
        }
      }
    }
  }
}


suppressWarnings({
  corr_df <- bind_rows(rows)
})


write.table(
  corr_df, output_file, 
  quote = FALSE, sep = "\t", row.names = FALSE
)