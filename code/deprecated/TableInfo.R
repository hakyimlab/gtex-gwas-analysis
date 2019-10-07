#######################################################################################
############## Definition of the BigQuery tables available to be queried ##############
#######################################################################################

# This script defines some objects containing the information 
# necessary to connect to each of the BQ tables.

tableInfo <- function(dataset="GTEx_V8_ElasticNet_EUR_2018_07_05", table="predixcan_results", project="gtex-awg-im") {
  info <- list()
  info$project <- project
  info$dataset_name <- dataset
  info$table_name <- table
  # add column names here?
  info
}

########################## DEFINITION OF BigQuery TABLES ########################

# elastic net models and gene-level associations
prediction_models_tbl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "prediction_models")
prediction_models_extra_tbl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "extra")
predixcan_tbl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "predixcan_results")
multixcan_tbl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "multixcan_results")

# conditional analysis (LDACC)
CA_eqtl_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "eqtl_analysis")
CA_gwas_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "gwas_results")
CA_eqtl_and_gwas_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "gwas_and_eqtl")

# DAPG
DAPG_eqtl_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "eqtl_analysis")
DAPG_gwas_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "gwas_results")
DAPG_eqtl_and_gwas_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "gwas_and_eqtl")

# colocalization results
coloc_tbl <- tableInfo("GTEx_V8_COLOC", "coloc_with_enloc_priors")
enloc_tbl <- tableInfo("GTEx_V8_ENLOC", "enloc_all_results")

# annotations and other metadata
ensembl_collapsed_annotations_tbl <- tableInfo("annotations", "ensembl_collapsed")
gene_essentiality_annotation_tbl <- tableInfo("annotations", "human_gene_essentiality_scores")
gencode_annotation_tbl <- tableInfo("annotations", "gencode_v26")
gtex_sample_size_tbl <- tableInfo("annotations", "sample_size")
gwas_metadata_tbl <- tableInfo("GTEx_V8_metadata", "gwas_metadata")
gtex_tissue_metadata_tbl <- tableInfo("GTEx_V8_metadata", "gtex_tissue_metadata")

# miscellaneous
ld_independent_regions_tbl <- tableInfo("miscellaneous", "ld_independent_regions")
gwas_catalog_tbl <- tableInfo("miscellaneous", "gwas_catalog_v102")