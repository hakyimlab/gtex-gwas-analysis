suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

FOLDER <-"output"
dir.create(FOLDER, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(FOLDER, p)

options(gargle_oauth_email = TRUE)

gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as phenotype, Deflation as deflation, new_abbreviation as abbreviation
                       FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}
                       WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

trait_whitelist <- "('GIANT_HEIGHT', 'CARDIoGRAM_C4D_CAD_ADDITIVE', 'UKB_20002_1111_self_reported_asthma')"

gwas_s <- (function(){
  query_ <- glue::glue("
SELECT f.phenotype, f.panel_variant_id, f.chromosome, f.zscore as formatted, i.zscore as imputed, i.imputation_status as status
FROM {gwas_formatted_tbl$dataset_name}.{gwas_formatted_tbl$table_name} as f
FULL OUTER JOIN {gwas_imputation_verification_tbl$dataset_name}.{gwas_imputation_verification_tbl$table_name} as i
ON f.phenotype = i.phenotype AND f.panel_variant_id = i.panel_variant_id
WHERE f.phenotype in {trait_whitelist} and i.imputation_status = 'imputed'
AND ((REGEXP_CONTAINS(f.panel_variant_id, 'C_G|G_C|A_T|T_A')) OR (REGEXP_CONTAINS(f.panel_variant_id, 'C_G|G_C|A_T|T_A')))
AND (f.chromosome = 'chr1' or i.chromosome =  'chr1')
") %>% glue::glue()
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
gwas_s <- gwas_s %>% mutate(phenotype=factor(phenotype, levels=c("GIANT_HEIGHT", "CARDIoGRAM_C4D_CAD_ADDITIVE", "UKB_20002_1111_self_reported_asthma")))


p_ <- ggplot(gwas_s) + theme_bw(base_size=20) +
  geom_point(aes(formatted, imputed), alpha=0.5) +
  geom_density_2d(aes(formatted, imputed)) +
  geom_abline(slope=1, intercept = 0) +
  facet_wrap(~phenotype) +
  xlab("Original") + ylab("Imputed") +
  ggtitle("Imputation agrees with original values", subtitle="Variants on chromosome 1")
save_plot(p_, fp_("GWAS_IMPUTATION_QUALITY.png"),500,1500)

p_ <- gwas_s %>% filter(phenotype %in% c("GIANT_HEIGHT", "CARDIoGRAM_C4D_CAD_ADDITIVE")) %>%
  ggplot() + theme_bw(base_size=20) +
  geom_point(aes(formatted, imputed), alpha=0.5) +
  geom_density_2d(aes(formatted, imputed)) +
  geom_abline(slope=1, intercept = 0) +
  facet_wrap(~phenotype, ncol = 1) +
  xlab("Original") + ylab("Imputed") +
  ggtitle("Imputation agrees with original values", subtitle="Variants on chromosome 1")
save_plot(p_, fp_("GWAS_IMPUTATION_QUALITY_2.png"),1000,500)
