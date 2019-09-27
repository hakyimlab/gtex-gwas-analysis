###############################################################################
# USE MINICONDA
###############################################################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

DATA<-"data/summaries"
dir.create(DATA, showWarnings = FALSE, recursive=TRUE)
dp_ <- function(p) file.path(DATA, p)

#########################################################################

gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as phenotype, Deflation as deflation from
                      {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name} WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)

})()
pheno_whitelist <- sprintf("(%s)", toString(sprintf("'%s'", gwas_metadata$phenotype)))

#########################################################################

get_if_not_exists <- function(path, method) {
  # if (file.exists(path)) {
  #   r_tsv_(path)
  # } else {
    d <- method()
    save_delim(d, path)
    d
  # }
  # method()
}

#########################################################################
get_multixcan_regions_eqtl <- function() {
  query_ <- glue::glue(
"WITH ranked_genes AS (
  SELECT *,
  COUNT(*) OVER (PARTITION BY region, phenotype) as count,
  ROW_NUMBER() OVER (PARTITION BY region, phenotype ORDER BY pvalue) as rk
    FROM ( SELECT r.region, m.*
      FROM ( SELECT m.*
             FROM {multixcan_tbl_eqtl$dataset_name}.{multixcan_tbl_eqtl$table_name} as m
             JOIN {multixcan_tbl_count_eqtl$dataset_name}.{multixcan_tbl_count_eqtl$table_name} as m_count
             ON m_count.phenotype=m.phenotype
             WHERE m.pvalue < m_count.b and m.phenotype in {pheno_whitelist}
           ) as m
      JOIN ( SELECT a.*, ld.region
        FROM {gencode_all_annotation_tbl$dataset_name}.{gencode_all_annotation_tbl$table_name} AS a
        JOIN {ld_independent_regions_2_tbl$dataset_name}.{ld_independent_regions_2_tbl$table_name} as ld ON (ld.start_location < a.start_location AND a.start_location <= ld.end_location AND a.chromosome = ld.chromosome)
        ) AS r
      ON m.gene = r.gene_id
    )
  )
  SELECT region, phenotype, count, pvalue as best_pvalue FROM ranked_genes  WHERE rk = 1"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}
m_eqtl_ <- get_if_not_exists(dp_("regions_multixcan_eqtl.txt"), get_multixcan_regions_eqtl)

get_multixcan_regions_sqtl <- function() {
  query_ <- glue::glue(
"WITH ranked_introns AS (
  SELECT *,
  COUNT(*) OVER (PARTITION BY region, phenotype) as count,
  ROW_NUMBER() OVER (PARTITION BY region, phenotype ORDER BY pvalue) as rk
    FROM ( SELECT r.region, m.*
      FROM ( SELECT m.*
             FROM {multixcan_tbl_sqtl$dataset_name}.{multixcan_tbl_sqtl$table_name} as m
             JOIN {multixcan_tbl_count_sqtl$dataset_name}.{multixcan_tbl_count_sqtl$table_name} as m_count
             ON m_count.phenotype=m.phenotype
             WHERE m.pvalue < m_count.b and m.phenotype in {pheno_whitelist}
           ) as m
      JOIN ( SELECT a.*, ld.region
        FROM {intron_annotation_tbl$dataset_name}.{intron_annotation_tbl$table_name} AS a
        JOIN {ld_independent_regions_2_tbl$dataset_name}.{ld_independent_regions_2_tbl$table_name}  as ld ON (ld.start_location < a.start_location AND a.start_location <= ld.end_location AND a.chromosome = ld.chromosome)
        ) AS r
      ON m.gene = r.intron_id
    )
  )
  SELECT region, phenotype, count, pvalue as best_pvalue FROM ranked_introns  WHERE rk = 1"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}
m_sqtl_ <- get_if_not_exists(dp_("regions_multixcan_sqtl.txt"), get_multixcan_regions_sqtl)

#########################################################################
get_predixcan_regions_eqtl <- function() {
  query_ <- glue::glue(
"WITH ranked_genes AS (
  SELECT *,
  COUNT(*) OVER (PARTITION BY region, phenotype) as count,
  ROW_NUMBER() OVER (PARTITION BY region, phenotype ORDER BY pvalue) as rk
    FROM (
      SELECT r.region, px.*
      FROM ( SELECT px.*
             FROM {predixcan_tbl_eqtl$dataset_name}.{predixcan_tbl_eqtl$table_name} as px
             JOIN {predixcan_tbl_count_eqtl$dataset_name}.{predixcan_tbl_count_eqtl$table_name} as px_count
             ON px_count.phenotype=px.phenotype
             WHERE px.pvalue < px_count.b and px.phenotype in {pheno_whitelist}
           ) as px
      JOIN ( SELECT a.*, ld.region
        FROM {gencode_all_annotation_tbl$dataset_name}.{gencode_all_annotation_tbl$table_name} AS a
        JOIN {ld_independent_regions_2_tbl$dataset_name}.{ld_independent_regions_2_tbl$table_name} as ld ON (ld.start_location < a.start_location AND a.start_location <= ld.end_location AND a.chromosome = ld.chromosome)
        ) AS r
      ON px.gene = r.gene_id
    )
  )
  SELECT region, phenotype, count, pvalue as best_pvalue FROM ranked_genes  WHERE rk = 1"
)
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}
p_eqtl_ <- get_if_not_exists(dp_("regions_predixcan_eqtl.txt"), get_predixcan_regions_eqtl)

get_predixcan_regions_sqtl <- function() {
  query_ <- glue::glue(
"WITH ranked_genes AS (
  SELECT *,
  COUNT(*) OVER (PARTITION BY region, phenotype) as count,
  ROW_NUMBER() OVER (PARTITION BY region, phenotype ORDER BY pvalue) as rk
    FROM (
      SELECT r.region, px.*
      FROM ( SELECT px.*
             FROM {predixcan_tbl_sqtl$dataset_name}.{predixcan_tbl_sqtl$table_name} as px
             JOIN {predixcan_tbl_count_sqtl$dataset_name}.{predixcan_tbl_count_sqtl$table_name} as px_count
             ON px_count.phenotype=px.phenotype
             WHERE px.pvalue < px_count.b and px.phenotype in {pheno_whitelist}
           ) as px
      JOIN ( SELECT a.*, ld.region
        FROM {intron_annotation_tbl$dataset_name}.{intron_annotation_tbl$table_name} AS a
        JOIN {ld_independent_regions_2_tbl$dataset_name}.{ld_independent_regions_2_tbl$table_name} as ld ON (ld.start_location < a.start_location AND a.start_location <= ld.end_location AND a.chromosome = ld.chromosome)
        ) AS r
      ON px.gene = r.intron_id
    )
  )
  SELECT region, phenotype, count, pvalue as best_pvalue FROM ranked_genes  WHERE rk = 1"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}
p_sqtl_ <- get_if_not_exists(dp_("regions_predixcan_sqtl.txt"), get_predixcan_regions_sqtl)

###############################################################################
get_enloc_regions_eqtl <- function() {
  query_ <- glue::glue(
"WITH ranked_genes AS (
  SELECT *,
  ROW_NUMBER() OVER (PARTITION BY region, phenotype ORDER BY rcp DESC) as rk,
  COUNT(*) OVER (PARTITION BY region, phenotype) as count
    FROM (
      SELECT r.region, e.*
      FROM ( SELECT molecular_qtl_trait as gene_id, locus_rcp as rcp, phenotype, tissue
             FROM {enloc_tbl_eqtl_eur$dataset_name}.{enloc_tbl_eqtl_eur$table_name} as e
             WHERE locus_rcp>0.5 ) as e
      JOIN ( SELECT a.*, ld.region
        FROM {gencode_all_annotation_tbl$dataset_name}.{gencode_all_annotation_tbl$table_name} AS a
        JOIN {ld_independent_regions_2_tbl$dataset_name}.{ld_independent_regions_2_tbl$table_name} as ld ON (ld.start_location < a.start_location AND a.start_location <= ld.end_location AND a.chromosome = ld.chromosome)
        ) AS r
      ON e.gene_id = r.gene_id
    )
  )
  SELECT region, phenotype, count, rcp as best_rcp FROM ranked_genes WHERE rk = 1"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}
e_eqtl_ <- get_if_not_exists(dp_("regions_enloc_eqtl.txt"), get_enloc_regions_eqtl)

get_enloc_regions_sqtl <- function() {
  query_ <- glue::glue(
"WITH ranked_genes AS (
  SELECT *,
  COUNT(*) OVER (PARTITION BY region, phenotype) as count,
  ROW_NUMBER() OVER (PARTITION BY region, phenotype ORDER BY rcp DESC) as rk
    FROM (
      SELECT r.region, e.*
      FROM ( SELECT molecular_qtl_trait as intron_id, locus_rcp as rcp, phenotype, tissue
             FROM {enloc_tbl_sqtl_eur$dataset_name}.{enloc_tbl_sqtl_eur$table_name} as e
             WHERE locus_rcp>0.5 ) as e
      JOIN ( SELECT a.*, ld.region
        FROM {intron_annotation_tbl$dataset_name}.{intron_annotation_tbl$table_name} AS a
        JOIN  {ld_independent_regions_2_tbl$dataset_name}.{ld_independent_regions_2_tbl$table_name} as ld ON (ld.start_location < a.start_location AND a.start_location <= ld.end_location AND a.chromosome = ld.chromosome)
        ) AS r
      ON e.intron_id = r.intron_id
    )
  )
  SELECT region, phenotype, count, rcp as best_rcp FROM ranked_genes WHERE rk = 1"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}
e_sqtl_ <- get_if_not_exists(dp_("regions_enloc_sqtl.txt"), get_enloc_regions_sqtl)

###############################################################################

get_gwas_regions <- function() {
  query_ <- glue::glue(
"WITH ranked_gwas AS (
  SELECT *,
  ROW_NUMBER() OVER (PARTITION BY region, phenotype ORDER BY pvalue) as rk,
  COUNT(*) OVER (PARTITION BY region, phenotype) as count
    FROM ( SELECT ld.region, g.*
           FROM ( SELECT g.* FROM `{gwas_tbl$dataset_name}.{gwas_tbl$table_name}` as g
                  JOIN `{gwas_tbl_count$dataset_name}.{gwas_tbl_count$table_name}` as g_count ON g_count.phenotype=g.phenotype
                  WHERE g.pvalue < g_count.b and g.phenotype in {pheno_whitelist}
                ) as g
          JOIN
            `{ld_independent_regions_2_tbl$dataset_name}.{ld_independent_regions_2_tbl$table_name}` as ld ON (ld.start_location < g.position AND g.position <= ld.end_location AND g.chromosome = ld.chromosome
          )
    )
  )
  SELECT region, phenotype, count, pvalue as best_pvalue FROM ranked_gwas WHERE rk = 1"
)
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
}
g_ <- get_if_not_exists(dp_("regions_gwas.txt"), get_gwas_regions)
message("##\nNumber of regions containing GWAS-Significant dectections: ", g_ %>% count(region) %>% nrow, "\n##")

gwas_only_ <- g_ %>% filter(count > 0) %>% select(region, phenotype)
tidy_up <- function(d, method) {
  d %>% inner_join(gwas_only_, by=c("region", "phenotype")) %>% select(region, phenotype, count) %>% mutate(method=method)
}

results <- tidy_up(g_, method="gwas") %>%
    rbind(tidy_up(m_eqtl_, "smultixcan_eqtl")) %>%
    rbind(tidy_up(m_sqtl_, "smultixcan_sqtl")) %>%
    rbind(tidy_up(p_eqtl_, "spredixcan_eqtl")) %>%
    rbind(tidy_up(p_sqtl_, "spredixcan_sqtl")) %>%
    rbind(tidy_up(e_eqtl_, "enloc_eqtl")) %>%
    rbind(tidy_up(e_sqtl_, "enloc_sqtl"))
results %>% save_delim(dp_("regions.txt"))

metadata_ <- (function(){
  gm_ <- (function(){
    query <- glue::glue(
"SELECT Tag as trait, new_abbreviation as abbreviation, Category
 FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}  WHERE Deflation=0")
    query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
  })()
  gm_ <- gm_ %>%
    mutate(Category=ifelse(Category == "Endocrine system disease", "Endocrine system", Category)) %>%
    mutate(Category=ifelse(Category == "Morphology", "Hair morphology", Category)) %>%
    mutate(Category=ifelse(Category == "Psychiatric_neurologic", "Psychiatric-neurologic", Category))

  tm_ <- (function(){
    query <- glue::glue(
"SELECT *
 FROM {trait_metadata_tbl$dataset_name}.{trait_metadata_tbl$table_name}")
    query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
  })()
  gm_ %>% inner_join(tm_, by=c("Category"="phenotype_class"))
})()

r_ <- results %>% select(region, n=count, trait=phenotype, results=method) %>%
  inner_join(metadata_ %>% rename(category=Category), by="trait")
r_ %>% save_delim(dp_("ld_region_count.txt"))
