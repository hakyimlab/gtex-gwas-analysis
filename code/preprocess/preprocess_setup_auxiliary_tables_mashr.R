###############################################################################
# WARNING! Don't run this script. It is meant to be run once per project, o set up tables.
# If you run it again, tables will get duplicate data.
# USE MINICONDA
###############################################################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(RCurl))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

options(gargle_oauth_email = TRUE)

###############################################################################

RESULT<-"results/kk"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

options(gargle_oauth_email = TRUE)

# (function(){
#     query <- glue::glue("
# SELECT p.gene, a.gene_type, count(*)
# FROM {predixcan_mashr_tbl_eqtl$dataset_name}.{predixcan_mashr_tbl_eqtl$table_name} as p
# INNER JOIN {gencode_all_annotation_tbl$dataset_name}.{gencode_all_annotation_tbl$table_name} AS a
# ON p.gene = a.gene_id
# GROUP BY p.gene, a.gene_type")
#     df <- query_exec(query, project = predixcan_mashr_tbl_eqtl$project, max_pages = Inf) %>% suppressWarnings()
# })()

#We avoid bq_table_create because it converts `snake_case` to `camelCase`
#Create a table with bonferroni thresholds
(function(){
  ds <- bq_dataset(predixcan_mashr_tbl_eqtl$project, predixcan_mashr_tbl_eqtl$dataset_name)

  (function(){
  query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {predixcan_mashr_tbl_eqtl$dataset_name}.{predixcan_mashr_tbl_eqtl$table_name}
                      WHERE pvalue IS NOT NULL
                      GROUP BY phenotype")
  df <- query_exec(query, project = predixcan_mashr_tbl_eqtl$project, max_pages = Inf) %>% suppressWarnings() %>%
    mutate(b=0.05/n) %>% select(phenotype, n, b)
  df %>% save_delim(fp_("expression_b_mashr.txt"))

  bq_ <- bq_table(predixcan_mashr_tbl_count_eqtl$project,predixcan_mashr_tbl_count_eqtl$dataset_name, predixcan_mashr_tbl_count_eqtl$table_name)
  bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
  bq_table_upload(bq_, df)
  })()

  (function(){
  query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {predixcan_mashr_tbl_sqtl$dataset_name}.{predixcan_mashr_tbl_sqtl$table_name}
                      WHERE pvalue IS NOT NULL
                      GROUP BY phenotype")
  df <- query_exec(query, project = predixcan_mashr_tbl_sqtl$project, max_pages = Inf) %>% suppressWarnings() %>%
    mutate(b=0.05/n) %>% select(phenotype, n, b)
  df %>% save_delim(fp_("splicing_b_mashr.txt"))

  bq_ <- bq_table(predixcan_mashr_tbl_count_sqtl$project, predixcan_mashr_tbl_count_sqtl$dataset_name, predixcan_mashr_tbl_count_sqtl$table_name)
  bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
  bq_table_upload(bq_, df)
  })()
})()

(function(){
  ds <- bq_dataset(multixcan_mashr_tbl_eqtl$project, multixcan_mashr_tbl_eqtl$dataset_name)

  (function(){
  query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {multixcan_mashr_tbl_eqtl$dataset_name}.{multixcan_mashr_tbl_eqtl$table_name}
                      WHERE pvalue IS NOT NULL
                      GROUP BY phenotype")
  df <- query_exec(query, project = multixcan_mashr_tbl_eqtl$project, max_pages = Inf) %>% suppressWarnings() %>%
    mutate(b=0.05/n) %>% select(phenotype, n, b)
  df %>% save_delim(fp_("expression_m_b.txt"))

  bq_ <- bq_table(multixcan_mashr_tbl_count_eqtl$project, multixcan_mashr_tbl_count_eqtl$dataset_name, multixcan_mashr_tbl_count_eqtl$table_name)
  bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
  bq_table_upload(bq_, df)
  })()

  (function(){
  query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {multixcan_mashr_tbl_sqtl$dataset_name}.{multixcan_mashr_tbl_sqtl$table_name}
                       WHERE pvalue IS NOT NULL
                       GROUP BY phenotype")
  df <- query_exec(query, project = multixcan_mashr_tbl_sqtl$project, max_pages = Inf) %>% suppressWarnings() %>%
    mutate(b=0.05/n) %>% select(phenotype, n, b)
  df %>% save_delim(fp_("splicing_m_b.txt"))

  bq_ <- bq_table(multixcan_mashr_tbl_count_sqtl$project, multixcan_mashr_tbl_count_sqtl$dataset_name, multixcan_mashr_tbl_count_sqtl$table_name)
  bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
  bq_table_upload(bq_, df)
  })()
})()
