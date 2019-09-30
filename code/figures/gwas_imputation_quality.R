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

gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as phenotype, Deflation as deflation, new_abbreviation as abbreviation
                       FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}
                       WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

trait_whitelist <- "('GIANT_HEIGHT', 'CARDIoGRAM_C4D_CAD_ADDITIVE', 'UKB_20002_1111_self_reported_asthma')"

# gwas_s <- (function(){
#   query_ <- glue::glue("
# SELECT f.phenotype, f.panel_variant_id, f.chromosome, f.zscore as formatted, i.zscore as imputed
# FROM {gwas_formatted_tbl$dataset_name}.{gwas_formatted_tbl$table_name} as f
# FULL OUTER JOIN {gwas_tbl$dataset_name}.{gwas_tbl$table_name} as i
# ON f.phenotype = i.phenotype AND f.panel_variant_id = i.panel_variant_id
# WHERE f.phenotype in {trait_whitelist}
# AND ((REGEXP_CONTAINS(f.panel_variant_id, 'C_G|G_C|A_T|T_A')) OR (REGEXP_CONTAINS(f.panel_variant_id, 'C_G|G_C|A_T|T_A')))
# AND (f.chromosome = 'chr1' or i.chromosome =  'chr1' or f.chromosome = 'chr6' or i.chromosome = 'chr6')
# ") %>% glue::glue()
#   query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
# })()
gwas_s <- (function(){
  query_ <- glue::glue("
SELECT phenotype, panel_variant_id, zscore, imputation_status
FROM {gwas_imputation_verification_tbl $dataset_name}.{gwas_imputation_verification_tbl $table_name} as f
WHERE f.phenotype in {trait_whitelist}
") %>% glue::glue()
  d <- query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
