suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
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


gwas_data_policies <- "data/gwas_data_policies.txt" %>% r_tsv_(col_types=cols_only(Phenotypes="c", Pubmed_paper_id="c")) %>%
  rename(phenotype= Phenotypes, id = Pubmed_paper_id)



gwas_metadata_ <- (function(){
  query <- glue::glue("
SELECT Tag as phenotype, new_abbreviation as abbreviation, Category as category,
Sample_Size as sample_size, Deflation as deflation, new_Phenotype as pheno_short, PUBMED_Paper_link as pubmed
FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

(function(){
  gwas_metadata_ %>% select(phenotype, id=pubmed) %>%  filter(!grepl("UKB", phenotype)) %>%
    left_join(gwas_data_policies %>% select(phenotype) %>% mutate(present=TRUE), by="phenotype") %>%
    filter(is.na(present)) %>% save_delim(fp_("_gwas_missing_pubmed.txt"))

})()
