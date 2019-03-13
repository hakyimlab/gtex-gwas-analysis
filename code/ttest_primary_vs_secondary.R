
# Author: Rodrigo Bonazzola
# This scripts runs t-tests to assess whether primary eQTLs (meaning those with greater effect size)
# have a larger contribution to phenotype than secondary eQTLs.

# install.packages("bigrquery")
library(bigrquery)

gtexv8repo_dir <- "~/Genomica/GTEx/gtex-v8"
bigquery_subdir <- "GWAS-subgroup/how-to/BigQuery"

setwd(file.path(gtexv8repo_dir, bigquery_subdir))
source("BigQuery.R")
source("Definitions.R")

myrepo_dir <- "~/Genomica/GTEx/gtex-gwas-analysis/"
work_subdir <- "data/"
workdir <- file.path(repo_dir, work_subdir)
setwd(workdir)

data_directory <- file.path(work_dir, "dapg_variants/rds_with_ldscore")
file_pattern <- "{data_directory}/DAPG_with_ldscore_by_phenotype_{tissue}.rds"

t.test_df <- data.frame("tissue"=character(), "phenotype"=character(), "t.test_pvalue"=numeric())

for (tissue in tissues) {
  
  message(glue::glue("Processing {tissue}..."))
  df <- readRDS(glue::glue(file_pattern))
  
  for (phenotype in traits) {

    message(glue::glue("  Processing {phenotype}..."))
    df_pheno <- df[[phenotype]] %>% mutate(gwas_effect_size_est=gwas_zscore/sqrt(maf*(1-maf)))
    
    n_eqtl <- table(df_pheno$gene_id)
    
    df_pheno <- df_pheno %>%
                .[,c("phenotype", "tissue", "gene_id", "eqtl_effect_size", "gwas_effect_size_est")] %>%
                group_by(gene_id) %>%
                mutate(rank_by_beta=order(order(-abs(eqtl_effect_size)))) %>%
                as.data.frame()
    
    n_eqtl <- table(df_pheno$gene_id)
    
    df_pheno <- subset(df_pheno, gene_id %in% names(n_eqtl[n_eqtl >= 2]))
    
    prim_eqtl <- df_pheno %>% .[, c("gene_id", "gwas_effect_size_est", "rank_by_beta")] %>% 
                              filter(rank_by_beta == 1) # %>% 
    prim_eqtl <- setNames(prim_eqtl[, "gwas_effect_size_est"], prim_eqtl[, "gene_id"])
    
    sec_eqtl <- df_pheno %>% .[, c("gene_id", "gwas_effect_size_est", "rank_by_beta")] %>% 
                             filter(rank_by_beta == 2) # %>% 
                             # setNames(.[, "gene_id"], .[, "gwas_effect_size"])
    
    sec_eqtl <- setNames( sec_eqtl[, "gwas_effect_size_est"], sec_eqtl[, "gene_id"])

    # to make sure the genes are in the same order in both
    sec_eqtl <- sec_eqtl[names(prim_eqtl)]
        
    tt <- t.test(abs(prim_eqtl), abs(sec_eqtl))
    n_row <- data.frame("tissue"=tissue, "phenotype"=phenotype, "t.test_pvalue"=tt$p.value)
    
    t.test_df <- rbind(t.test_df, n_row)
    
  }
}

saveRDS(t.test_df, "t_test_primary_vs_secondary.rds")