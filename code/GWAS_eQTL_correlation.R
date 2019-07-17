library(tidyverse)
library(yaml)

data_type <- "expression"
dapg_file_pattern <- "../data/dapg_selected_variants/expression/gwas_and_eqtl/DAPG_with_mashr__{tissue}.rds"
metadata <- read_tsv("../data/gwas_metadata.txt")
tissues <- readLines("../data/tissue_list.txt")
phenotypes <- readLines("../data/phenotype_list.txt")
use_top_n <- c(FALSE, TRUE)
shuffle <- c(FALSE, TRUE)

worst_tissue <- "Kidney_Cortex"
worst_tissue_df <- readRDS(glue::glue(gsub("tissue", "worst_tissue", dapg_file_pattern)))

output_file <- "../output/eqtl_gwas_correlation_with_ci.tsv"

rows <- vector(mode="list", length=0)

for (tissue in tissues) {
  message(tissue)
  dfs <- readRDS(glue::glue(dapg_file_pattern))
  
  for (phenotype in phenotypes) {
    for (use_top_n_ in use_top_n) {
      for (shuffle_ in shuffle) {
        df <- dfs[[phenotype]]
        
        if (use_top_n_) {
          n_variants <-nrow(worst_tissue_df[[phenotype]])
          df <- df %>% top_n(n_variants, abs(eqtl_effect_size/slope_se))
        }
            
        if (shuffle_) {
          corr_test <- cor.test(
            abs(df$eqtl_effect_size),
            sample(abs(df$gwas_effect_size)),
            use = "complete.obs"
          )
          corr_sp <- cor(abs(df$eqtl_effect_size), sample(abs(df$gwas_effect_size)), use = "complete.obs", method="spearman")
        } else {
          corr_test <- cor.test(
            abs(df$eqtl_effect_size),
            abs(df$gwas_effect_size),
            use = "complete.obs"
          )
          corr_sp <- cor(abs(df$eqtl_effect_size), abs(df$gwas_effect_size), use = "complete.obs", method="spearman")
        }
        
        new_row <- data.frame(
          "tissue" = tissue,
          "phenotype" = phenotype,
          "pearson" = corr_test$estimate,
          "ci95_lower" = corr_test$conf.int[1],
          "ci95_upper" = corr_test$conf.int[2],
          "pvalue" = corr_test$p.value,
          "data" = ifelse(shuffle_, "shuffled", "real"),
          "type" = data_type,
          "top_n" = nrow(df)
        )
        
        rows <- c(rows, list(new_row))
        
      }
    }
  }
}

suppressWarnings({
  corr_df <- bind_rows(rows)
})


write.table(corr_df, output_file, quote = FALSE, sep = "\t", row.names = FALSE)