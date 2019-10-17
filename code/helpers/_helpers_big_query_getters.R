get_gwas_metadata <- function() {
  query <- glue::glue("SELECT Tag as phenotype, Deflation as deflation, new_abbreviation as abbreviation, Category as category
                       FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}
                       WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}

sql_pheno_whitelist <- function(phenos) {
  sprintf("(%s)", toString(sprintf("'%s'", phenos)))
}

get_gencode <- function(){
  query <- glue::glue("SELECT *
                       FROM {gencode_all_annotation_tbl$dataset_name}.{gencode_all_annotation_tbl$table_name}
                       WHERE gene_type in ('lincRNA', 'protein_coding', 'pseudogene')")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}

get_introns <- function(){
  query <- glue::glue("SELECT * FROM {intron_annotation_tbl$dataset_name}.{intron_annotation_tbl$table_name}")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}

###############################################################################

gwas_metadata <- get_gwas_metadata()
pheno_whitelist_ <- gwas_metadata$phenotype
pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist_)


################################################################################################################################

get_e_count <- function(e_tbl = enloc_tbl_eqtl_eur, threshold = NULL,  pheno_whitelist_=pheno_whitelist_){

  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist_)

  query <- "
SELECT COUNT(DISTINCT molecular_qtl_trait)
FROM {e_tbl$dataset_name}.{e_tbl$table_name} as e
WHERE e.phenotype in {pheno_whitelist}
" %>% glue::glue()

  if (!is.null(threshold)) {
    if (threshold < 0) stop("need non-negative rcp threshold")
    if (threshold > 1) stop("need rcp threshold below one")
    if (threshold > 0) {
      query <- query %>% glue::glue(" AND e.locus_rcp > {threshold}")
    }
  }

  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf) %>% as.integer
}


get_e_count_for_thresholds <- function(e_tbl, thresholds, e_threshold=NULL, pheno_whitelist=pheno_whitelist_) {
  spg <- list()
  for (i in 1:length(thresholds)) {
    spg[[i]] <- get_e_count(e_tbl=e_tbl, threshold=thresholds[i], pheno_whitelist=pheno_whitelist)
  }
  unlist(spg)
}

###############################################################################

get_spe_count <- function(sp_tbl = predixcan_mashr_tbl_eqtl, e_tbl = enloc_tbl_eqtl_eur,  s_threshold = NULL, e_threshold = NULL, pheno_whitelist=pheno_whitelist_){
  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist)
  query <- "
  SELECT COUNT(DISTINCT sp.gene)
  FROM {sp_tbl$dataset_name}.{sp_tbl$table_name} as sp" %>% glue::glue()

  if (!is.null(e_threshold)) {
    if (e_threshold < 0) stop("need non-negative enloc threshold")
    if (e_threshold > 1) stop("need enloc threshold below 1")
    if (e_threshold > 0) {
      query <- query %>% glue::glue("
                                    INNER JOIN ( SELECT phenotype, tissue, molecular_qtl_trait as gene, locus_rcp as rcp FROM {e_tbl$dataset_name}.{e_tbl$table_name}
                                    ) as e
                                    ON sp.phenotype = e.phenotype AND sp.tissue=e.tissue AND sp.gene = e.gene
                                    ")
    }

    }

  query <- query %>% glue::glue("
                                WHERE sp.phenotype in {pheno_whitelist}")


  if (!is.null(s_threshold)) {
    if (s_threshold < 0) stop("need non-negative significance threshold")
    if (s_threshold > 1) stop("need significance threshold below one")
    if (s_threshold < 1) {
      query <- query %>% glue::glue(" AND sp.pvalue < {s_threshold}")
    }
  }

  if (!is.null(e_threshold) && e_threshold > 0) {
    query <- query %>% glue::glue(" AND rcp > {e_threshold}")
  }

  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf) %>% as.integer
}

get_spe_count_for_thresholds <- function(sp_tbl, e_tbl, thresholds, e_threshold=NULL, pheno_whitelist=pheno_whitelist_) {
  spg <- list()
  for (i in 1:length(thresholds)) {
    spg[[i]] <- get_spe_count(sp_tbl=sp_tbl, e_tbl=e_tbl, s_threshold=thresholds[i], e_threshold=e_threshold, pheno_whitelist=pheno_whitelist)
  }
  unlist(spg)
}

###############################################################################

get_bonferroni <- function(c_tbl, pheno_whitelist=pheno_whitelist_){
  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist_)
  query <- "
SELECT *
FROM {c_tbl$dataset_name}.{c_tbl$table_name} as b
WHERE b.phenotype in {pheno_whitelist}
" %>% glue::glue()
 query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}


