source("code/BigQuery.R", chdir = TRUE)

df <- readRDS("data/dapg_variants/with_mashr_results/DAPG_with_mashr__Whole_Blood.rds")

gwas_metadata <- gtex_gwas_metadata() %>% 
                 select("Tag", "Sample_Size", "Deflation") %>%
                 rename("phenotype"="Tag")

sigma_y_est_df <- data.frame("phenotype"=character(),
                             "mean_sigma_y"=numeric(),
                             "sd_sigma_y"=numeric(),
                             "mean_log10_sigma_y"=numeric(),
                             "sd_log10_sigma_y"=numeric())


for (phenotype_ in names(df)) {
  
  sample_size <- gwas_metadata %>% filter(phenotype==phenotype_) %>% .[1,"Sample_Size"]
  is_deflated <- (gwas_metadata %>% filter(phenotype==phenotype_) %>% .[1,"Deflation"]) != 0
  
  if (is_deflated)
    next
  
  df_pheno <- df[[phenotype_]]
  
  df_pheno_ <- df_pheno %>% filter(maf > 0.4 & !is.na(gwas_se) & is.finite(gwas_se) & abs(gwas_zscore) < 2)
  
  # EAGLE ECZEMA 
  if (phenotype_ == "EAGLE_Eczema") {
    df_pheno_ <- df_pheno_ %>% filter(gwas_se < 0.1)
  }
  
  sigma_y <- sqrt(2 * sample_size * df_pheno_$maf * (1-df_pheno_$maf) * df_pheno_$gwas_se**2)
  mean_sigma_y <- mean(sigma_y)
  sd_sigma_y <- sd(sigma_y)

  log10_sigma_y <- log10(2 * sample_size) + log10(df_pheno_$maf) + log10(1-df_pheno_$maf) + 2 * log10(df_pheno_$gwas_se)
  mean_log10_sigma_y <- mean(log10_sigma_y)
  sd_log10_sigma_y <- sd(log10_sigma_y)
  
  new_row <- data.frame("phenotype" = phenotype_,
                        "mean_sigma_y" = mean_sigma_y,
                        "sd_sigma_y" = sd_sigma_y,
                        "ratio" = mean_sigma_y/sd_sigma_y,
                        "mean_log10_sigma_y" = mean_log10_sigma_y,
                        "sd_log10_sigma_y" = sd_log10_sigma_y)
  
  sigma_y_est_df <- rbind(sigma_y_est_df, new_row)
  
  # n_total <- nrow(df_pheno)
  # n_used <- nrow(df_pheno_)
  # print(glue::glue("{phenotype_}: {n_used}/{n_total}"))
}


saveRDS(sigma_y_est_df, "output/estimated_sigma_y.rds")