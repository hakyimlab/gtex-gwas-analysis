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

get_gtex_metadata <- function(){
  query <- glue::glue("SELECT tissue,
v6_all, v6_eur, v7_all, v7_eur, v8_all, v8_eur,
tissue_name  as name, tissue_abbrv as abbreviation, tissue_color_hex as color
FROM {gtex_tissue_metadata_tbl$dataset_name}.{gtex_tissue_metadata_tbl$table_name}")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}

###############################################################################

gwas_metadata <- get_gwas_metadata()
pheno_whitelist_ <- gwas_metadata$phenotype
pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist_)


###############################################################################

get_models_count <- function(e_tbl) {
  query <- "
SELECT COUNT(DISTINCT gene) as models, tissue
FROM {e_tbl$dataset_name}.{e_tbl$table_name}
GROUP BY tissue
" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}

###############################################################################

get_gene_count_in_sp <- function(sp_tbl = predixcan_en_tbl_eqtl) {
  query <- "
SELECT COUNT(DISTINCT gene)
FROM {sp_tbl$dataset_name}.{sp_tbl$table_name}
" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}

get_selected_gene_count <- function(g_tbl = gencode_all_annotation_tbl) {
  query <- "
SELECT COUNT(DISTINCT g.gene_id)
FROM {g_tbl$dataset_name}.{g_tbl$table_name} as g
WHERE g.gene_type  in ('lincRNA', 'protein_coding', 'pseudogene')
" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}

get_selected_gene_count_in_sp <- function(g_tbl = gencode_all_annotation_tbl, sp_tbl = predixcan_en_tbl_eqtl) {
  query <- "
SELECT COUNT(DISTINCT g.gene_id)
FROM {g_tbl$dataset_name}.{g_tbl$table_name} as g
INNER JOIN (
  SELECT gene
  FROM {sp_tbl$dataset_name}.{sp_tbl$table_name}
  GROUP BY gene
) as sp ON sp.gene = g.gene_id
WHERE g.gene_type  in ('lincRNA', 'protein_coding', 'pseudogene')
" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
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

################################################################################################################################

get_e_count <- function(e_tbl = enloc_tbl_eqtl_eur, threshold = NULL,  pheno_whitelist_=pheno_whitelist_, restrict_to_g_tbl=NULL){
  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist_)

  query <- "
SELECT COUNT(DISTINCT e.molecular_qtl_trait)
FROM {e_tbl$dataset_name}.{e_tbl$table_name} as e
" %>% glue::glue()

  if (!is.null(restrict_to_g_tbl)) {
    query <- glue::glue(query,"
INNER JOIN
  {restrict_to_g_tbl$dataset_name}.{restrict_to_g_tbl$table_name} as g
ON
  g.gene_id = e.molecular_qtl_trait
")
  }

  query <- query %>% glue::glue("
WHERE e.phenotype in {pheno_whitelist}")

  if (!is.null(restrict_to_g_tbl)) {
    query <- query %>% glue::glue(" AND
g.gene_type  in ('lincRNA', 'protein_coding', 'pseudogene')"
    )
  }

  if (!is.null(threshold)) {
    if (threshold < 0) stop("need non-negative rcp threshold")
    if (threshold > 1) stop("need rcp threshold below one")
    if (threshold > 0) {
      query <- query %>% glue::glue(" AND
e.locus_rcp > {threshold}")
    }
  }

  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf) %>% as.integer
}


get_e_count_for_thresholds <- function(e_tbl, thresholds, e_threshold=NULL, pheno_whitelist=pheno_whitelist_, restrict_to_g_tbl=NULL) {
  d <- list()
  for (i in 1:length(thresholds)) {
    d[[i]] <- get_e_count(e_tbl=e_tbl, threshold=thresholds[i], pheno_whitelist=pheno_whitelist, restrict_to_g_tbl=restrict_to_g_tbl)
  }
  unlist(d)
}

################################################################################################################################

get_sp_count <- function(g_tbl, sp_tbl, threshold = NULL,  pheno_whitelist_=pheno_whitelist_){
  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist_)

  query <- "
SELECT COUNT(DISTINCT g.gene_id)
FROM
  {g_tbl$dataset_name}.{g_tbl$table_name} as g
INNER JOIN
  {sp_tbl$dataset_name}.{sp_tbl$table_name} as sp
ON
  g.gene_id = sp.gene
WHERE
  g.gene_type  in ('lincRNA', 'protein_coding', 'pseudogene') AND
  sp.phenotype in {pheno_whitelist}
" %>% glue::glue()

  if (!is.null(threshold)) {
    if (threshold < 0) stop("need non-negative significance threshold")
    if (threshold > 1) stop("need significance threshold below one")
    if (threshold < 1) {
      query <- query %>% glue::glue(" AND sp.pvalue < {threshold}")
    }
  }

  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf) %>% as.integer
}


get_sp_count_for_thresholds <- function(g_tbl, sp_tbl, thresholds, pheno_whitelist=pheno_whitelist_) {
  d <- list()
  for (i in 1:length(thresholds)) {
    d[[i]] <- get_sp_count(g_tbl=g_tbl, sp_tbl=sp_tbl, threshold=thresholds[i], pheno_whitelist=pheno_whitelist)
  }
  unlist(d)
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

get_count_genes_with_sp_intron <- function(g_tbl, igm_tbl, sp_tbl, threshold = NULL, pheno_whitelist=pheno_whitelist_) {
  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist)
  query <- "
SELECT COUNT(DISTINCT g.gene_id)
FROM {g_tbl$dataset_name}.{g_tbl$table_name} as g
INNER JOIN {igm_tbl$dataset_name}.{igm_tbl$table_name} as igm
  ON g.gene_id = igm.gene_id
INNER JOIN {sp_tbl$dataset_name}.{sp_tbl$table_name} as sp
  ON sp.gene = igm.intron_id
WHERE
  g.gene_type  in ('lincRNA', 'protein_coding', 'pseudogene') AND
  sp.phenotype in {pheno_whitelist}
" %>% glue::glue()

  if (!is.null(threshold)) {
    if (threshold < 0) stop("need non-negative significance threshold")
    if (threshold > 1) stop("need significance threshold below one")
    if (threshold < 1) {
      query <- query %>% glue::glue(" AND
sp.pvalue < {threshold}")
    }
  }

  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)

}

get_count_genes_with_sp_intron_for_thresholds <- function(g_tbl, igm_tbl, sp_tbl, thresholds, pheno_whitelist=pheno_whitelist_) {
  d <- list()
  for (i in 1:length(thresholds)) {
    d[[i]] <- get_count_genes_with_sp_intron(g_tbl, igm_tbl, sp_tbl, threshold=thresholds[i], pheno_whitelist=pheno_whitelist)
  }
  unlist(d)
}

###############################################################################

get_count_genes_with_e_intron <- function(g_tbl, igm_tbl, e_tbl, threshold = NULL, pheno_whitelist=pheno_whitelist_) {
  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist)
  query <- "
SELECT COUNT(DISTINCT g.gene_id)
FROM {g_tbl$dataset_name}.{g_tbl$table_name} as g
INNER JOIN {igm_tbl$dataset_name}.{igm_tbl$table_name} as igm
  ON g.gene_id = igm.gene_id
INNER JOIN {e_tbl$dataset_name}.{e_tbl$table_name} as e
  ON e.molecular_qtl_trait = igm.intron_id
WHERE
  g.gene_type  in ('lincRNA', 'protein_coding', 'pseudogene') AND
  e.phenotype in {pheno_whitelist}
" %>% glue::glue()

  if (!is.null(threshold)) {
    if (threshold < 0) stop("need non-negative significance threshold")
    if (threshold > 1) stop("need significance threshold below one")
    if (threshold > 0) {
      query <- query %>% glue::glue(" AND
e.locus_rcp > {threshold}")
    }
  }

  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)

}

get_count_genes_with_e_intron_for_thresholds <- function(g_tbl, igm_tbl, e_tbl, thresholds, pheno_whitelist=pheno_whitelist_) {
  d <- list()
  for (i in 1:length(thresholds)) {
    d[[i]] <- get_count_genes_with_e_intron(g_tbl, igm_tbl, e_tbl, threshold=thresholds[i], pheno_whitelist=pheno_whitelist)
  }
  unlist(d)
}

###############################################################################

get_count_genes_with_spe_intron <- function(g_tbl, igm_tbl, sp_tbl, e_tbl, s_threshold = NULL, c_threshold, pheno_whitelist=pheno_whitelist_) {
  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist)
  query <- "
SELECT COUNT(DISTINCT g.gene_id)
FROM {g_tbl$dataset_name}.{g_tbl$table_name} as g
INNER JOIN {igm_tbl$dataset_name}.{igm_tbl$table_name} as igm
  ON g.gene_id = igm.gene_id
INNER JOIN {sp_tbl$dataset_name}.{sp_tbl$table_name} as sp
  ON sp.gene = igm.intron_id
INNER JOIN {e_tbl$dataset_name}.{e_tbl$table_name} as e
  ON sp.gene = e.molecular_qtl_trait AND sp.phenotype = e.phenotype AND sp.tissue = e.tissue
WHERE
  g.gene_type  in ('lincRNA', 'protein_coding', 'pseudogene') AND
  sp.phenotype in {pheno_whitelist}
" %>% glue::glue()

  if (!is.null(s_threshold)) {
    if (s_threshold < 0) stop("need non-negative significance threshold")
    if (s_threshold > 1) stop("need significance threshold below one")
    if (s_threshold < 1) {
      query <- query %>% glue::glue(" AND
sp.pvalue < {s_threshold}")
    }
  }

  if (!is.null(c_threshold)) {
    if (c_threshold < 0) stop("need non-negative colocalization threshold")
    if (c_threshold > 1) stop("need colocalization threshold below one")
    if (c_threshold < 1) {
      query <- query %>% glue::glue(" AND
e.locus_rcp > {c_threshold}")
    }
  }

  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)

}

get_count_genes_with_spe_intron_for_thresholds <- function(g_tbl, igm_tbl, sp_tbl, e_tbl, s_thresholds, c_threshold, pheno_whitelist=pheno_whitelist_) {
  d <- list()
  for (i in 1:length(s_thresholds)) {
    d[[i]] <- get_count_genes_with_spe_intron(g_tbl, igm_tbl, sp_tbl, e_tbl, s_threshold=s_thresholds[i], c_threshold= c_threshold, pheno_whitelist=pheno_whitelist)
  }
  unlist(d)
}

#######################################################################################

get_spredixcan_significant_genes <- function(sp_tbl, sp_count_tbl, pheno_whitelist=pheno_whitelist_) {
    pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist)
  query <- "
SELECT COUNT(*)
FROM (
  SELECT DISTINCT gene
  FROM ( SELECT px.gene
         FROM {sp_tbl$dataset_name}.{sp_tbl$table_name} as px
         JOIN {sp_count_tbl$dataset_name}.{sp_count_tbl$table_name} as px_count
         ON px_count.phenotype=px.phenotype
         WHERE px.pvalue < px_count.b and px.phenotype in {pheno_whitelist}
  )
)" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf) %>% as.integer
}

get_spredixcan_significant_genes_by_trait <- function(sp_tbl, sp_count_tbl, pheno_whitelist=pheno_whitelist_) {
  pheno_whitelist <- sql_pheno_whitelist(pheno_whitelist)
  query <- "
SELECT px.phenotype,  COUNT(*) as count
FROM {sp_tbl$dataset_name}.{sp_tbl$table_name} as px
JOIN {sp_count_tbl$dataset_name}.{sp_count_tbl$table_name} as px_count
ON px_count.phenotype=px.phenotype
WHERE px.pvalue < px_count.b and px.phenotype in {pheno_whitelist}
GROUP BY px.phenotype" %>% glue::glue()
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}
