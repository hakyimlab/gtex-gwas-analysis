suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))


sqtl <- (function(){
  query <- "
SELECT COUNT(*) FROM  {enloc_tbl_sqtl_eur$dataset_name}.{enloc_tbl_sqtl_eur$table_name}" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

eqtl <- (function(){
  query <- "
SELECT COUNT(*) FROM  {enloc_tbl_eqtl_eur$dataset_name}.{enloc_tbl_eqtl_eur$table_name}" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
