suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))


sqtl <- (function(){
  query <- "
SELECT *
FROM  {enloc_tbl_sqtl_eur$dataset_name}.{enloc_tbl_sqtl_eur$table_name} as e
WHERE phenotype =  'UKB_20022_Birth_weight' and tissue='Artery_Coronary' and molecular_qtl_trait='intron_18_47930446_47930800'" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
