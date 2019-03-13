#install.packages("ggplot2")
library(ggplot2)

#install.packages("dplyr")
library(dplyr)

#install.packages("RCurl")
library(RCurl)

#install.packages("bigrquery")
library(bigrquery)

repo_dir <- "~/Genomica/GTEx/gtex-v8"
bigquery_subdir <- "GWAS-subgroup/how-to/BigQuery/"
work_dir <- "~/Genomica/GTEx/v8/analysis"

setwd(file.path(repo_dir, bigquery_subdir))
source("BigQuery.R")
setwd(work_dir)

gtex_tissue_metadata_df <- get_gtex_tissue_metadata() %>%
  mutate(tissue_color=sapply(tissue_color_hex, function(x) paste0("#", x))) %>%
  select(-tissue_color_hex) %>%
  mutate(tissue_color=as.character(tissue_color))

gtex_gwas_metadata_df <- gtex_gwas_metadata() %>%
  select(Tag, new_Phenotype, abbreviation, Deflation) %>%
  rename("phenotype"="Tag")

# Read pre-computed counts of LD blocks with GWAS hits from a file
read_file_from_link <- function(link, sep="\t") read.delim(text=getURL(link), sep=sep, header=TRUE, stringsAsFactors=FALSE)

# gwas_ld_blocks_counts_link <- "https://storage.googleapis.com/gtex-exchange/miscellanea/gwas_pred_1e-6_ldblock_count.txt"
# 
# gwas_ld_blocks_counts <- read_file_from_link(gwas_ld_blocks_counts_link) %>%
#   select(trait_id, tissue_id, is_gwas_ldblock_count) %>%
#   rename("phenotype"="trait_id", "tissue"="tissue_id", "gwas_ld_block_count"="is_gwas_ldblock_count")

get_gwas_blocks <- function(phenotype) {
  phenotype <- gsub("\\.|-", "_", phenotype)
  query <- "WITH ranked_genes AS (
      SELECT *, ROW_NUMBER() OVER (PARTITION BY region_name ORDER BY pvalue) as rk
      FROM (
        SELECT ld.region_name, gwas.*
        FROM `gtex-awg-im.GWAS.{phenotype}` as gwas
        JOIN `gtex-awg-im.miscellaneous.ld_independent_regions` as ld
        ON (ld.start < gwas.position AND gwas.position < ld.end AND ld.chromosome = gwas.chromosome)
      )
    )
    SELECT region_name, pvalue as best_pvalue FROM ranked_genes WHERE rk = 1"
  
  df <- query_exec(glue::glue(query), 
                   project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
  df <- df %>% mutate(phenotype=phenotype)
  df
}

# For each trait, tissue and region, get the best gene-level p-value
# The download could take around 10 minutes
get_predixcan_blocks <- function() {
  query_exec(
    "WITH ranked_genes AS (
      SELECT *, ROW_NUMBER() OVER (PARTITION BY region_name, phenotype, tissue ORDER BY pvalue) as rk
      FROM (
      SELECT 
        r.region_name, 
        px.*
        FROM `gtex-awg-im.GTEx_V8_ElasticNet_EUR_2018_07_05.predixcan_results` as px 
        JOIN (
          SELECT a.*, ld.region_name 
          FROM `gtex-awg-im.annotations.gencode_v26` AS a
          JOIN `gtex-awg-im.miscellaneous.ld_independent_regions` as ld
          ON (ld.start < a.start AND a.start < ld.end AND a.chromosome = ld.chromosome)
        ) AS r
        ON SUBSTR(px.gene, 1, 15) = r.gene_id
        ) 
    )
    SELECT region_name, phenotype, tissue, pvalue as best_pvalue FROM ranked_genes WHERE rk = 1", 
    project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf
  )
}

get_multixcan_blocks <- function() {
  query_exec(
    "
    WITH ranked_genes AS (
      SELECT *, ROW_NUMBER() OVER (PARTITION BY region_name, phenotype ORDER BY pvalue) as rk
      FROM (
        SELECT r.region_name, mx.*
        FROM `gtex-awg-im.GTEx_V8_ElasticNet_EUR_2018_07_05.multixcan_results` as mx
        JOIN (
          SELECT a.*, ld.region_name 
          FROM `gtex-awg-im.annotations.gencode_v26` AS a
          JOIN `gtex-awg-im.miscellaneous.ld_independent_regions` as ld
          ON (ld.start < a.start AND a.start < ld.end AND a.chromosome = ld.chromosome)
        ) AS r
        ON SUBSTR(mx.gene, 1, 15) = r.gene_id
      ) 
    )

    SELECT region_name, phenotype, pvalue as best_pvalue FROM ranked_genes WHERE rk = 1", 
    project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf
  )
}

get_enloc_blocks <- function() {
  query_exec(
    "
    WITH ranked_genes AS (
      SELECT *, ROW_NUMBER() OVER (PARTITION BY region_name, phenotype ORDER BY rcp DESC) as rk
      FROM (
        SELECT 
        r.region_name, 
        enloc.*
        FROM `gtex-awg-im.GTEx_V8_ENLOC.enloc_all_results` as enloc 
        JOIN (
          SELECT a.*, ld.region_name 
          FROM `gtex-awg-im.annotations.gencode_v26` AS a
          JOIN `gtex-awg-im.miscellaneous.ld_independent_regions` as ld
          ON (ld.start < a.start AND a.start < ld.end AND a.chromosome = ld.chromosome)
        ) AS r
        ON SUBSTR(enloc.gene_id, 1, 15) = r.gene_id
      ) 
    )
    SELECT region_name, phenotype, rcp as best_rcp FROM ranked_genes WHERE rk = 1", 
    project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf
  )
}

# get the gene count per tissue to estimate a Bonferroni p-value threshold for each tissue
get_gene_count_per_tissue <- function() {
  query_exec(
    "SELECT tissue, COUNT(*) as gene_count
    FROM `gtex-awg-im.GTEx_V8_ElasticNet_EUR_2018_07_05.extra`
    GROUP BY tissue", project = "gtex-awg-im", use_legacy_sql = FALSE
  )
}


################################################
####  GWAS: best hits per Pickrell region.  ####
################################################
pheno_lst <- getListOfPhenotypes(predixcan_tbl) #[1:5]
if ( !("gwas_best_hit_per_region" %in% ls() && nrow(gwas_best_hit_per_region) == 194004)) {
  gwas_best_hit_per_region <- data.frame(region_name=character(), best_pvalue=numeric(), phenotype=character())
  for(pheno in pheno_lst) {
    dd <- get_gwas_blocks(pheno)  
    gwas_best_hit_per_region <- rbind(gwas_best_hit_per_region, dd)
  }
  gwas_best_hit_per_region <- gwas_best_hit_per_region %>% rename(best_gwas_pvalue=best_pvalue)
} 
########################

########################
if (!"px_df" %in% ls()) # to avoid making the query if we want to re-run the whole script
  px_df <- get_predixcan_blocks()
gene_count_per_tissue_df <- get_gene_count_per_tissue()
########################

########################
mx_df <- get_multixcan_blocks()
gene_count_across_tissues <- query_exec("SELECT COUNT(DISTINCT gene) FROM GTEx_V8_ElasticNet_EUR_2018_07_05.extra", 
                                        project = "gtex-awg-im", max_pages = Inf, use_legacy_sql = F)[1,1] # %>% .$gene %>% unique() %>% length()
mx_pvalue_threshold <- 0.05/gene_count_across_tissues
mx_df <- mx_df %>% rename("best_mx_pvalue"="best_pvalue")
########################

########################
enloc_df <- get_enloc_blocks()
########################

best_per_block_df <- inner_join(gwas_best_hit_per_region, mx_df)
best_per_block_df <- inner_join(best_per_block_df, enloc_df)

# px_ld_blocks_counts <- px_df %>%
#   inner_join(gene_count_per_tissue_df, by="tissue") %>%
#   mutate(pval_threshold=0.05/gene_count) %>%
#   mutate(is_predixcan_signif = best_pvalue<pval_threshold) %>%
#   group_by(phenotype, tissue) %>%
#   summarise(predixcan_signif_ld_block_count=sum(is_predixcan_signif))

gwas_hits_counts <- best_per_block_df %>% group_by(phenotype) %>% 
                    filter(best_gwas_pvalue < 1e-6) %>% 
                    summarise(gwas_counts=n())

gwas_mx_counts <- best_per_block_df %>% group_by(phenotype) %>% 
                  filter(best_gwas_pvalue < 1e-6 & best_mx_pvalue < 0.05/gene_count_across_tissues) %>%
                  summarise(gwas_mx_counts=n())

gwas_mx_enloc_counts <- best_per_block_df %>% group_by(phenotype) %>% 
                        filter(best_gwas_pvalue < 1e-6 & best_mx_pvalue < 0.05/gene_count_across_tissues & best_rcp > 0.5) %>% 
                        summarise(gwas_mx_enloc_counts=n())

all_counts <- full_join(gwas_hits_counts, gwas_mx_counts)
all_counts <- full_join(all_counts, gwas_mx_enloc_counts)
all_counts <- inner_join(all_counts, gtex_gwas_metadata_df)