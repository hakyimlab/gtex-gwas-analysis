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
gwas_tbl <-                         tableInfo("GWAS_all", "gwas")
gwas_tbl_count <-                   tableInfo("GWAS_all", "gwas_results_count")
gwas_formatted_tbl <-               tableInfo("GWAS_all", "formatted_gwas")
gwas_imputation_verification_tbl <- tableInfo("GWAS_all", "gwas_imputation_verification")

#elastic net predixcan GTEX v7
v7_prediction_en_models_tbl_eqtl <-       tableInfo("GTEx_V7_HapMap_2017_11_29", "weights")
v7_prediction_en_models_extra_tbl_eqtl <- tableInfo("GTEx_V7_HapMap_2017_11_29", "extra")

#elastic net predixcan
prediction_en_models_tbl_eqtl <-       tableInfo("GTEx_V8_ElasticNet_EUR_v1", "weights_eqtl")
prediction_en_models_extra_tbl_eqtl <- tableInfo("GTEx_V8_ElasticNet_EUR_v1", "extra_eqtl")

prediction_en_models_tbl_sqtl <-       tableInfo("GTEx_V8_ElasticNet_EUR_v1", "weights_sqtl")
prediction_en_models_extra_tbl_sqtl <- tableInfo("GTEx_V8_ElasticNet_EUR_v1", "extra_sqtl")

predixcan_en_tbl_eqtl <-               tableInfo("GTEx_V8_ElasticNet_EUR_v1", "spredixcan_eqtl")
predixcan_en_tbl_count_eqtl <-         tableInfo("GTEx_V8_ElasticNet_EUR_v1", "spredixcan_eqtl_count")

predixcan_en_tbl_sqtl <-               tableInfo("GTEx_V8_ElasticNet_EUR_v1", "spredixcan_sqtl")
predixcan_en_tbl_count_sqtl <-         tableInfo("GTEx_V8_ElasticNet_EUR_v1", "spredixcan_sqtl_count")

multixcan_en_tbl_eqtl <-               tableInfo("GTEx_V8_ElasticNet_EUR_v1", "smultixcan_eqtl")
multixcan_en_tbl_count_eqtl <-         tableInfo("GTEx_V8_ElasticNet_EUR_v1", "smultixcan_eqtl_count")

multixcan_en_tbl_sqtl <-               tableInfo("GTEx_V8_ElasticNet_EUR_v1", "smultixcan_sqtl")
multixcan_en_tbl_count_sqtl <-         tableInfo("GTEx_V8_ElasticNet_EUR_v1", "smultixcan_sqtl_count")

#elastic net predixcan without palindromic
predixcan_en_np_tbl_eqtl <-               tableInfo("GTEx_V8_ElasticNet_EUR_v1", "spredixcan_eqtl_np")
predixcan_en_np_tbl_count_eqtl <-         tableInfo("GTEx_V8_ElasticNet_EUR_v1", "spredixcan_eqtl_np_count")

#mashr predixcan
prediction_mashr_models_tbl_eqtl <-       tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "weights_eqtl")
prediction_mashr_models_extra_tbl_eqtl <- tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "extra_eqtl")

prediction_mashr_models_tbl_sqtl <-       tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "weights_sqtl")
prediction_mashr_models_extra_tbl_sqtl <- tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "extra_sqtl")

predixcan_mashr_tbl_eqtl <-               tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "spredixcan_eqtl")
predixcan_mashr_tbl_count_eqtl <-         tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "spredixcan_eqtl_count")

predixcan_mashr_tbl_sqtl <-               tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "spredixcan_sqtl")
predixcan_mashr_tbl_count_sqtl <-         tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "spredixcan_sqtl_count")

multixcan_mashr_tbl_eqtl <-               tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "smultixcan_eqtl")
multixcan_mashr_tbl_count_eqtl <-         tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "smultixcan_eqtl_count")

multixcan_mashr_tbl_sqtl <-               tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "smultixcan_sqtl")
multixcan_mashr_tbl_count_sqtl <-         tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "smultixcan_sqtl_count")

#mashr with harmonized(unimputed) gwas
predixcan_mashr_hq_tbl_eqtl <-               tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "spredixcan_eqtl_hq")
predixcan_mashr_hq_tbl_count_eqtl <-         tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "spredixcan_eqtl_hq_count")

predixcan_mashr_hn_tbl_eqtl <-               tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "spredixcan_eqtl_hn")
predixcan_mashr_hn_tbl_count_eqtl <-         tableInfo("GTEx_V8_PF_MASHR_EUR_v1", "spredixcan_eqtl_hn_count")

#EN-DAPGW predixcan
prediction_en_dapgw_models_tbl_eqtl <-       tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "weights_eqtl")
prediction_en_dapgw_models_extra_tbl_eqtl <- tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "extra_eqtl")

prediction_en_dapgw_models_tbl_sqtl <-       tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "weights_sqtl")
prediction_en_dapgw_models_extra_tbl_sqtl <- tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "extra_sqtl")

predixcan_en_dapgw_tbl_eqtl <-               tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "spredixcan_eqtl")
predixcan_en_dapgw_tbl_count_eqtl <-         tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "spredixcan_eqtl_count")

predixcan_en_dapgw_tbl_sqtl <-               tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "spredixcan_sqtl")
predixcan_en_dapgw_tbl_count_sqtl <-         tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "spredixcan_sqtl_count")

multixcan_en_dapgw_tbl_eqtl <-               tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "smultixcan_eqtl")
multixcan_en_dapgw_tbl_count_eqtl <-         tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "smultixcan_eqtl_count")

multixcan_en_dapgw_tbl_sqtl <-               tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "smultixcan_sqtl")
multixcan_en_dapgw_tbl_count_sqtl <-         tableInfo("GTEx_V8_ENDAPGW_EUR_v1", "smultixcan_sqtl_count")

#CTIMP
prediction_ctimp_models_tbl_eqtl <-       tableInfo("GTEx_V8_PF_CTIMP_EUR_v1", "weights_eqtl")
prediction_ctimp_models_extra_tbl_eqtl <- tableInfo("GTEx_V8_PF_CTIMP_EUR_v1", "extra_eqtl")

predixcan_ctimp_tbl_eqtl <-               tableInfo("GTEx_V8_PF_CTIMP_EUR_v1", "spredixcan_eqtl")
predixcan_ctimp_tbl_count_eqtl <-         tableInfo("GTEx_V8_PF_CTIMP_EUR_v1", "spredixcan_eqtl_count")

#CTIMP without palindromic
predixcan_ctimp_np_tbl_eqtl <-               tableInfo("GTEx_V8_PF_CTIMP_EUR_v1", "spredixcan_eqtl_np")
predixcan_ctimp_np_tbl_count_eqtl <-         tableInfo("GTEx_V8_PF_CTIMP_EUR_v1", "spredixcan_eqtl_np_count")

# conditional analysis (LDACC)
CA_eqtl_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "eqtl_analysis")
CA_gwas_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "gwas_results")
CA_eqtl_and_gwas_tbl <- tableInfo("GTEx_V8_ConditionalAnalysis_2018_10_05", "gwas_and_eqtl")

# DAPG
#DAPG_eqtl_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "eqtl_analysis")
DAPG_eqtl_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "eqtl_analysis")
DAPG_gwas_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "gwas_results")
DAPG_eqtl_and_gwas_tbl <- tableInfo("GTEx_V8_DAPG_2018_10_05", "gwas_and_eqtl")

DAPG_eqtl_clusters <- tableInfo("GTEx_V8_DAPG_EUR_v1", "clusters_eqtl")
DAPG_sqtl_clusters <- tableInfo("GTEx_V8_DAPG", "clusters_sqtl")

# colocalization results
coloc_tbl_eqtl <- tableInfo("GTEx_V8_COLOC", "coloc_with_enloc_priors")
enloc_tbl_eqtl <- tableInfo("GTEx_V8_ENLOC", "enloc_all_results")
enloc_tbl_eqtl_eur <- tableInfo("GTEx_V8_ENLOC_v1", "enloc_eqtl_eur")
enloc_tbl_sqtl_eur <- tableInfo("GTEx_V8_ENLOC_v1", "enloc_sqtl_eur")

# annotations and other metadata
ensembl_collapsed_annotations_tbl <-  tableInfo("annotations", "ensembl_collapsed")
gene_essentiality_annotation_tbl <-   tableInfo("annotations", "human_gene_essentiality_scores")
gencode_all_annotation_tbl <-         tableInfo("annotations", "gencode_v26_all")
gencode_annotation_tbl <-             tableInfo("annotations", "gencode_v26")
intron_annotation_tbl <-              tableInfo("annotations", "introns")
gtex_sample_size_tbl <-               tableInfo("annotations", "sample_size")
intron_gene_mapping_tbl <-            tableInfo("annotations", "intron_gene_map")
gwas_metadata_tbl <-                  tableInfo("GTEx_V8_metadata", "gwas_metadata")
trait_metadata_tbl <-                 tableInfo("GTEx_V8_metadata", "phenotype_classes_colors")
gtex_tissue_metadata_tbl <-           tableInfo("GTEx_V8_metadata", "gtex_tissue_metadata")

# miscellaneous
ld_independent_regions_tbl <-   tableInfo("miscellaneous", "ld_independent_regions")
ld_independent_regions_2_tbl <- tableInfo("annotations", "ld_independent_regions_2")
gwas_catalog_tbl <-             tableInfo("miscellaneous", "gwas_catalog_v102")
