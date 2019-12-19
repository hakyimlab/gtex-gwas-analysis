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

sp_count <- function(sp_tbl, spc_tbl, store) {
  query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {sp_tbl$dataset_name}.{sp_tbl$table_name}
                        WHERE pvalue IS NOT NULL
                        GROUP BY phenotype")
  df <- query_exec(query, project = sp_tbl$project, max_pages = Inf) %>% suppressWarnings() %>%
    mutate(b=0.05/n) %>% select(phenotype, n, b)
  df %>% save_delim(fp_(store))

  bq_ <- bq_table(spc_tbl$project, spc_tbl$dataset_name, spc_tbl$table_name)
  bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
  bq_table_upload(bq_, df)
}

sp_count(predixcan_ctimp_np_tbl_eqtl, predixcan_ctimp_np_tbl_count_eqtl, "expression_b_ctimpnp.txt")
sp_count(predixcan_en_np_tbl_eqtl, predixcan_en_np_tbl_count_eqtl, "expression_b_ennp.txt")
