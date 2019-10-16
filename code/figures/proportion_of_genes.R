suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(ggplot2))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_plot_style.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

DATA<-"data/summaries"
dir.create(DATA, showWarnings = FALSE, recursive=TRUE)
dp_ <- function(p) file.path(DATA, p)

FOLDER <-"output"
dir.create(FOLDER, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(FOLDER, p)

options(gargle_oauth_email = TRUE)

###############################################################################

gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as phenotype, Deflation as deflation, new_abbreviation as abbreviation, Category as category
                       FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}
                       WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
pheno_whitelist <- sprintf("(%s)", toString(sprintf("'%s'", gwas_metadata$phenotype)))

gencode <- (function(){
  query <- glue::glue("SELECT *
                       FROM {gencode_all_annotation_tbl$dataset_name}.{gencode_all_annotation_tbl$table_name}
                       WHERE gene_type in ('lincRNA', 'protein_coding', 'pseudogene')")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

introns <- (function(){
  query <- glue::glue("SELECT * FROM {intron_annotation_tbl$dataset_name}.{intron_annotation_tbl$table_name}")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

###############################################################################

get_spredixcan_count <- function(sp_tbl = predixcan_mashr_tbl_eqtl, e_tbl = enloc_tbl_eqtl_eur,  s_threshold = NULL, e_threshold = NULL){
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

get_sp_count_for_thresholds <- function(thresholds, e_threshold=NULL) {
  spg <- list()
  for (i in 1:length(thresholds)) {
    spg[[i]] <- get_spredixcan_count(s_threshold=thresholds[i], e_threshold=e_threshold)
  }
  unlist(spg)
}

###############################################################################

n_sp_eqtl_entries <- get_spredixcan_count()
thresholds <- c(0.05, 0.05/49) %>% c(10^(-1*(0:30))) %>% sort(decreasing=TRUE)

data <- read_or_get(dp_("spredixcan_proportions.txt"), function(){
spg_e0 <- get_sp_count_for_thresholds(thresholds, NULL)
spg_e0.5 <- get_sp_count_for_thresholds(thresholds, 0.5)

data.frame(threshold = thresholds, n = spg_e0, prop = spg_e0/spg_e0[1], enloc = 0) %>%
  rbind(data.frame(threshold = thresholds, n = spg_e0.5, prop = spg_e0.5/spg_e0.5[1], enloc = 0.5))
})

build_decor_ <- function(data) {
  k_ <- data %>% filter(enloc == 0)
  p_at_0.05 <- approx(k_$prop, k_$threshold, 0.05)$y
  data.frame(x=c(0,    p_at_0.05, p_at_0.05),
                       y=c(0.05, 0.05,      0))
}
decor_ <- build_decor_(data)

p_ <- data %>% mutate(t = ifelse(enloc == 0.5, "genes", "enloc genes")) %>%
  ggplot(aes(x=-log10(threshold), y=prop, color=t)) + theme_bw(base_size = 25) + paper_base_theme_ +
  geom_point() + geom_line() +
  geom_line(aes(x=x,y=y), decor_, show.legend = FALSE, linetype=2, color="gray") +
  xlab(expression(paste(-log[10],'[',italic(P),']',sep=""))) +
  ylab('proportion of genes') +
  scale_x_continuous(breaks=c(0,5,10,20,-log10(p_at_0.05),30),labels=c(0,5,10,20,round(-log10(p_at_0.05),digits=2),30),expand=c(0.01,0))+
  scale_y_continuous(limits = c(NA,1.05), breaks=c(0,0.05,0.25,0.5,0.75,1.00),labels=c(0,0.05,0.25,0.5,0.75,1.00),expand=c(0,0.01))

save_plot(p_, fp_("spe_proportion_genes.png"), 600, 600)


