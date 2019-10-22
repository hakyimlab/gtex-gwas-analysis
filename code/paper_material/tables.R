suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(xtable))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

options(gargle_oauth_email = TRUE)

RESULT<-"output/tables"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

###############################################################################

count <- (function(){
  query <- glue::glue("
SELECT phenotype, tissue, count(*) as count FROM `{predixcan_mashr_tbl_eqtl$dataset_name}.{predixcan_mashr_tbl_eqtl$table_name}`
group by phenotype, tissue order by phenotype, tissue")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as phenotype, new_abbreviation as abbreviation, Category as category, Sample_Size as sample_size
                       FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}
                       WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

order_ <- gwas_metadata %>% filter(!grepl("UKB", phenotype)) %>% group_by(category) %>% arrange(-sample_size) %>% slice(1) %>% ungroup %>% arrange(-sample_size) %>% .$category
all_ <- unique(gwas_metadata$category)
order_ <- order_ %>% c(all_[!(all_ %in% order_)])

t <- gwas_metadata %>% mutate(category=factor(category, levels=order_)) %>% arrange(category, phenotype) %>%
  mutate(category=gsub("_", "-", category), phenotype=gsub("_", " ", phenotype), abbreviation=gsub("_", "\\_", abbreviation), sample_size=as.integer(sample_size))
t %>% select(Category=category, Phenotype=phenotype, Abbreviation=abbreviation, `Sample Size`=sample_size) %>% xtable %>% print.xtable(file = fp_("gwas_table.tex"), include.rownames=FALSE)

###############################################################################

tc_ <- function(tbl) {
  query <- glue::glue("
SELECT tissue, count(*) as count FROM `{tbl$dataset_name}.{tbl$table_name}`
group by tissue")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}

gtex_ <- (function(){
  query <- glue::glue("
SELECT tissue, tissue_name as name, v8_eur as european_samples, tissue_abbrv as abbreviation FROM `{gtex_tissue_metadata_tbl$dataset_name}.{gtex_tissue_metadata_tbl$table_name}`")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})() %>% rename(`european samples` = european_samples)

eqtl_n  <- tc_(prediction_mashr_models_extra_tbl_eqtl) %>% rename(`expression models` = count)
sqtl_n  <- tc_(prediction_mashr_models_extra_tbl_sqtl) %>% rename(`splicing models` = count)

gtex_ <- gtex_ %>%
  inner_join(eqtl_n, by="tissue") %>%
  inner_join(sqtl_n, by="tissue")

gtex_%>% arrange(tissue) %>% select(-tissue) %>%  xtable %>% print.xtable(file = fp_("gtex_models_table.tex"), include.rownames=FALSE)
