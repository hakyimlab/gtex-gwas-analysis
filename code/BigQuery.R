source("TableInfo.R")

#install.packages("bigrquery")
library(bigrquery)

#install.packages("tidyverse")
library(tidyverse)

##########################################################################################################

getListOfPhenotypes <- function(table, filter_deflated=T) {
  if (filter_deflated)
    phenoLst <- query_exec(glue::glue("SELECT Tag as phenotype FROM {table$dataset_name}.{table$table_name} WHERE Deflation = 0"), table$project)
  else
    phenoLst <- query_exec(glue::glue("SELECT Tag as phenotype FROM {table$dataset_name}.{table$table_name}"), table$project)
  phenoLst[[1]]
}

getListOfTissues <- function(table) {
  tissueLst <- query_exec(glue::glue("SELECT tissue FROM {table$dataset_name}.{table$table_name} GROUP by tissue"), table$project)
  tissueLst[[1]]
}

# perform a basic select query on a BigQuery table, retrieving all (default) or some of the columns
basicQuery <- function(table, columns=NULL, fake_submit=FALSE, max_pages = 1e4) {
  if (!is.null(columns)) { 
    columns <- paste(columns, collapse = ", ")
  } else {
    columns <- "*"
  }
  query <- glue::glue("SELECT {columns} FROM {table$dataset_name}.{table$table_name}")
  if (fake_submit) {
    return(query)
  } else {
    return(query_exec(query, table$project, max_pages = max_pages))
  }
}

# select a set of rows, matching the given phenotype and/or tissue
queryByPhenotypeAndTissue <- function(table, columns=NULL, pheno_pattern=NULL, tissue_pattern=NULL, additional_condition=NULL, show_query=FALSE, fake_submit=FALSE, max_pages=1e4) {
  select_part <- basicQuery(table, columns, fake_submit = TRUE)
  if (!is.null(pheno_pattern) || !is.null(tissue_pattern) || !is.null(additional_condition)) {
    # phenotypes
    where_part <- paste0("WHERE")
    if (!is.null(pheno_pattern)) {
      pheno_conditions <- character()
      for (pheno in pheno_pattern) {
        pheno_conditions <- c(pheno_conditions, glue::glue("phenotype LIKE '{pheno}'"))
      }
      pheno_conditions <- paste0(pheno_conditions, collapse = " OR ")
      where_part <- glue::glue("{where_part} ({pheno_conditions})")
    }
    # tissues
    if (!is.null(tissue_pattern)) {
      if (where_part != "WHERE") 
        where_part <- glue::glue("{where_part} AND")
      tissue_conditions <- character()
      for (tissue in tissue_pattern) {
        tissue_conditions <- c(tissue_conditions, glue::glue("tissue LIKE '{tissue}'"))
      }
      tissue_conditions <- paste0(tissue_conditions, collapse = " OR ")
      where_part <- glue::glue("{where_part} ({tissue_conditions})")
    }
    # Additional condition
    if (!is.null(additional_condition)) {
      if (where_part != "WHERE") 
        where_part <- glue::glue("{where_part} AND")
      where_part <- glue::glue("{where_part} {additional_condition}")
    }
  }
  
  sql_query <- paste(select_part, where_part)
  if (show_query)
    message(sql_query)
  
  if (!fake_submit) {
    df <- query_exec(sql_query, table$project, max_pages = max_pages)
    return(df)
  } else {
    return(sql_query)
  }
}

queryModels <- function(table=prediction_models_tbl, genes, gene_column = "gene", tissue_pattern=NULL, show_query=FALSE, fake_submit=FALSE, max_pages=1e4) {
  select_part <- glue::glue("SELECT * FROM {table$dataset_name}.{table$table_name}")
  where_part <- "WHERE"
  
  if (!is.null(tissue_pattern)) {
    if (where_part != "WHERE") where_part <- glue::glue("{where_part} AND ")
    where_part <- glue::glue("{where_part} tissue LIKE '{tissue_pattern}'")
    tissue_conditions <- character()
    for (tissue in tissue_pattern) {
      tissue_conditions <- c(tissue_conditions, glue::glue("tissue LIKE '{tissue}'"))
    }
    tissue_conditions <- paste0(tissue_conditions, collapse = " OR ")
    where_part <- glue::glue("{where_part} {tissue_conditions}'")
  }
  
  if(!is.null(genes)) {
    if (where_part != "WHERE") where_part <- glue::glue("{where_part} AND ")
    gene_conditions <- character()
    for (gene in genes) {
      gene_conditions <- c(gene_conditions, glue::glue("{gene_column} LIKE '{gene}'"))
    }
    gene_conditions <- paste0(gene_conditions, collapse = " OR ")
    where_part <- glue::glue("{where_part} {gene_conditions}")
  }
  
  sql_query <- paste(select_part, where_part)
  
  if (show_query) {
    message(sql_query)
  }
  
  if (!fake_submit) {
    df <- query_exec(sql_query, table$project, max_pages=max_pages)
    return(df)
  } else {
    return(sql_query)
  }
}

get_gtex_tissue_metadata <- function() {
  basicQuery(gtex_tissue_metadata_tbl)
}

gtex_gwas_metadata <- function() {
  basicQuery(gwas_metadata_tbl)
}

# source("Definitions.R")