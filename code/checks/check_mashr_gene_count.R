suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

options(gargle_oauth_email = TRUE)
suppressWarnings(source("code/helpers/_helpers_big_query_getters.R"))

selected_genes <- get_selected_gene_count()
selected_genes_in_mashr <- get_selected_gene_count_in_sp()
genes_in_mashr <- get_gene_count_in_sp()
