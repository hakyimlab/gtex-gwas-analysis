# derived from Rodrigo's code

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
gwas_tbl <- tableInfo("GWAS_all", "gwas")

prediction_models_tbl_eqtl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "prediction_models")
prediction_models_extra_tbl_eqtl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "extra")

prediction_models_tbl_sqtl <- tableInfo("GTEX_V8_ElasticNet_EUR_Splicing_2018_11_19", "models_weights")
prediction_models_extra_tbl_sqtl <- tableInfo("GTEX_V8_ElasticNet_EUR_Splicing_2018_11_19", "models_extra")

predixcan_tbl_eqtl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "predixcan_results")
predixcan_tbl_count_eqtl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "predixcan_results_count")

predixcan_tbl_sqtl <- tableInfo("GTEX_V8_ElasticNet_EUR_Splicing_2018_11_19", "predixcan_results")
predixcan_tbl_count_sqtl <- tableInfo("GTEX_V8_ElasticNet_EUR_Splicing_2018_11_19", "predixcan_results_count")

multixcan_tbl_eqtl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "multixcan_results")
multixcan_tbl_count_eqtl <- tableInfo("GTEx_V8_ElasticNet_EUR_2018_07_05", "multixcan_results_count")

multixcan_tbl_sqtl <- tableInfo("GTEX_V8_ElasticNet_EUR_Splicing_2018_11_19", "multixcan_results")
multixcan_tbl_count_sqtl <- tableInfo("GTEX_V8_ElasticNet_EUR_Splicing_2018_11_19", "multixcan_results_count")

# conditional analysis (LDACC)
CA_eqtl_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "eqtl_analysis")
CA_gwas_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "gwas_results")
CA_eqtl_and_gwas_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "gwas_and_eqtl")

# DAPG
DAPG_eqtl_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "eqtl_analysis")
DAPG_gwas_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "gwas_results")
DAPG_eqtl_and_gwas_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "gwas_and_eqtl")

# colocalization results
coloc_tbl_eqtl <- tableInfo("GTEx_V8_COLOC", "coloc_with_enloc_priors")
enloc_tbl_eqtl <- tableInfo("GTEx_V8_ENLOC", "enloc_all_results")
enloc_tbl_eqtl_eur <- tableInfo("GTEx_V8_ENLOC", "enloc_eqtl_eur")
enloc_tbl_sqtl_eur <- tableInfo("GTEx_V8_ENLOC", "merged_enloc_sqtl")

# annotations and other metadata
ensembl_collapsed_annotations_tbl <- tableInfo("annotations", "ensembl_collapsed")
gene_essentiality_annotation_tbl <- tableInfo("annotations", "human_gene_essentiality_scores")
gencode_all_annotation_tbl <- tableInfo("annotations", "gencode_v26_all")
gencode_annotation_tbl <- tableInfo("annotations", "gencode_v26")
intron_annotation_tbl <- tableInfo("annotations", "introns")
gtex_sample_size_tbl <- tableInfo("annotations", "sample_size")
gwas_metadata_tbl <- tableInfo("GTEx_V8_metadata", "gwas_metadata")
gtex_tissue_metadata_tbl <- tableInfo("GTEx_V8_metadata", "gtex_tissue_metadata")

# miscellaneous
ld_independent_regions_tbl <- tableInfo("miscellaneous", "ld_independent_regions")
gwas_catalog_tbl <- tableInfo("miscellaneous", "gwas_catalog_v102")
