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

###############################################################################

RESULT<-"results/kk"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

# ### Check if BigQuery has the same I can load. Answer: no.
# ### Both BigQuery and my code seem to lose something in the way.
# (function(){
#   eqtl_predixcan <- get_predixcan_results_logic()
#   read_ <- function(trait) {
#     get_predixcan_result(eqtl_predixcan, trait, col_types=cols_only(gene="c", tissue="c",pvalue="d")) %>% filter(!is.na(pvalue))
#   }
#   p_ <- read_("UKB_20002_1265_self_reported_migraine")
#
#   query <- glue::glue("SELECT gene,tissue, pvalue FROM {predixcan_tbl_eqtl$dataset_name}.{predixcan_tbl_eqtl$table_name}
#                       WHERE pvalue IS NOT NULL and phenotype = 'UKB_20002_1265_self_reported_migraine'")
#   df <- query_exec(query, project = predixcan_tbl_eqtl$project, max_pages = Inf) %>% suppressWarnings()
#
#   k <- p_ %>% full_join(df, by=c("gene", "tissue"))
# })()

# #We avoid bq_table_create because it converts `snake_case` to `camelCase`
# # Create a table with bonferroni thresholds
# (function(){
#   ds <- bq_dataset(predixcan_tbl_eqtl$project, predixcan_tbl_eqtl$dataset_name)
#
#   (function(){
#   query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {predixcan_tbl_eqtl$dataset_name}.{predixcan_tbl_eqtl$table_name}
#                       WHERE pvalue IS NOT NULL
#                       GROUP BY phenotype")
#   df <- query_exec(query, project = predixcan_tbl_eqtl$project, max_pages = Inf) %>% suppressWarnings() %>%
#     mutate(b=0.05/n) %>% select(phenotype, n, b)
#   df %>% save_delim(fp_("expression_b.txt"))
#
#   bq_ <- bq_table(predixcan_tbl_count_eqtl$project,predixcan_tbl_count_eqtl$dataset_name, predixcan_tbl_count_eqtl$table_name)
#   bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
#   bq_table_upload(bq_, df)
#   })()
#
#   (function(){
#   query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {predixcan_tbl_sqtl$dataset_name}.{predixcan_tbl_sqtl$table_name}
#                       WHERE pvalue IS NOT NULL
#                       GROUP BY phenotype")
#   df <- query_exec(query, project = predixcan_tbl_sqtl$project, max_pages = Inf) %>% suppressWarnings() %>%
#     mutate(b=0.05/n) %>% select(phenotype, n, b)
#   df %>% save_delim(fp_("splicing_b.txt"))
#
#   bq_ <- bq_table(predixcan_tbl_count_sqtl$project, predixcan_tbl_count_sqtl$dataset_name, predixcan_tbl_count_sqtl$table_name)
#   bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
#   bq_table_upload(bq_, df)
#   })()
# })()

# (function(){
#   ds <- bq_dataset(multixcan_tbl_eqtl$project, multixcan_tbl_eqtl$dataset_name)
#
#   (function(){
#   query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {multixcan_tbl_eqtl$dataset_name}.{multixcan_tbl_eqtl$table_name}
#                       WHERE pvalue IS NOT NULL
#                       GROUP BY phenotype")
#   df <- query_exec(query, project = multixcan_tbl_eqtl$project, max_pages = Inf) %>% suppressWarnings() %>%
#     mutate(b=0.05/n) %>% select(phenotype, n, b)
#   df %>% save_delim(fp_("expression_m_b.txt"))
#
#   bq_ <- bq_table(multixcan_tbl_count_eqtl$project, multixcan_tbl_count_eqtl$dataset_name, multixcan_tbl_count_eqtl$table_name)
#   bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
#   bq_table_upload(bq_, df)
#   })()
#
#   (function(){
#   query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {multixcan_tbl_sqtl$dataset_name}.{multixcan_tbl_sqtl$table_name}
#                        WHERE pvalue IS NOT NULL
#                        GROUP BY phenotype")
#   df <- query_exec(query, project = multixcan_tbl_sqtl$project, max_pages = Inf) %>% suppressWarnings() %>%
#     mutate(b=0.05/n) %>% select(phenotype, n, b)
#   df %>% save_delim(fp_("splicing_m_b.txt"))
#
#   bq_ <- bq_table(multixcan_tbl_count_sqtl$project, multixcan_tbl_count_sqtl$dataset_name, multixcan_tbl_count_sqtl$table_name)
#   bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
#   bq_table_upload(bq_, df)
#   })
# })()

# (function(){
#   ds <- bq_dataset(gwas_tbl$project, gwas_tbl$dataset_name)
#
#   query <- glue::glue("SELECT count(pvalue) as n, phenotype FROM {gwas_tbl$dataset_name}.{gwas_tbl$table_name}
#                      WHERE pvalue IS NOT NULL
#                      GROUP BY phenotype")
#   df <- query_exec(query, project = gwas_tbl$project, max_pages = Inf) %>% suppressWarnings() %>%
#     mutate(b=0.05/n) %>% select(phenotype, n, b)
#   df %>% save_delim(fp_("gwas_b.txt"))
#
#   bq_ <- bq_table(gwas_tbl$project, gwas_tbl$dataset_name, "gwas_results_count")
#   bq_ <- bq_table_create(bq_, fields=as_bq_fields(df), friendly_name = "Number of results per trait")
#   bq_table_upload(bq_, df)
#
# })()
#
# (function(){
#   files <- list.files("data/splicing_introns")
#   r <- list()
#   for (i in 1:length(files)) {
#     message(files[i])
#     d <- file.path("data/splicing_introns", files[i]) %>% r_tsv_(col_types=cols_only(intron_id="c", chromosome="c", start_location="i", end_location="i"))
#     r[[i]] <- d
#   }
#   r <- do.call(rbind,r) %>% unique
#
#   ds <- bq_dataset(intron_annotation_tbl$project, intron_annotation_tbl$dataset_name)
#
#
#   bq_ <- bq_table(intron_annotation_tbl$project, intron_annotation_tbl$dataset_name, intron_annotation_tbl$table_name)
#   bq_ <- bq_table_create(bq_, fields=as_bq_fields(r), friendly_name = "intron specification")
#   bq_table_upload(bq_, r)
# })()
#
