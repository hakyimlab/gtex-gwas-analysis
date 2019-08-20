suppressWarnings(library(dplyr))
suppressWarnings(library(readr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(cowplot))
suppressWarnings(library(ggrepel))
#suppressWarnings(library(RCurl))
#suppressWarnings(library(bigrquery))
suppressWarnings(library(RSQLite))

suppressWarnings(source("_helpers.R"))
suppressWarnings(source("_helpers_data.R"))
suppressWarnings(source("_helpers_qq.R"))
suppressWarnings(source("_data_adapter.R"))

RESULT<-"data/stats/enloc_multi"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

RESULT_SQTL<-"data/stats/enloc_multi_sqtl"
dir.create(RESULT_SQTL, showWarnings = FALSE, recursive=TRUE)
fps_ <- function(p) file.path(RESULT_SQTL, p)

###############################################################################

gtex_metadata <- "data/gtex_metadata.csv" %>% r_csv_
gwas_metadata <- "data/gwas_metadata.txt" %>% r_tsv_

multixcan_eqtl <- g_results_logic("data/smp_v8", "smultixcan_imputed_gwas_gtexelv8_(.*)_ccn30.txt")
multixcan_sqtl <- g_results_logic("data/smp_v8_splicing", "smultixcan_imputed_gwas_gtexelv8_splicing_(.*)_ccn30.txt")

enloc_eqtl <- get_enloc_results_logic()
enloc_sqtl <- get_enloc_results_logic("v8_splicing")

regions <- "data/eur_ld.annot.gz" %>% r_tsv_(col_types = cols_only(chromosome="c", start_location="i", end_location="i", gene_id="c")) %>%
  rename(region=gene_id) %>% mutate(chromosome = as.integer(gsub("chr", "", chromosome))) %>%
  filter(!is.na(chromosome), !is.na(start_location), !is.na(end_location))

keep_gwas <- gwas_metadata %>% filter(Deflation == 0) %>% .$Tag
gi <- g_results_logic("data/gwas", "imputed_(.*).txt.gz") %>% filter(trait %in% keep_gwas)

gencode_annotation <- "data/gencode_v26_stranded.txt" %>% r_tsv_ %>%
  select(id=gene_id, chromosome, position=start_location) %>%
  mutate(chromosome = as.integer(gsub("chr(.*)", "\\1", chromosome))) %>%
  filter(!is.na(chromosome))

intron_annotation <- "data/intron_annotation.txt.gz" %>% r_tsv_ %>%
  select(id=gene_id, chromosome, position=start_location) %>%
  mutate(chromosome = as.integer(gsub("chr(.*)", "\\1", chromosome))) %>%
  filter(!is.na(chromosome))

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) { 
  row <- as.integer(args[1])
  cat ("Kept row ", row, "\n")
  gi <- gi[row,]
} else { 
  row <- NULL 
}

###############################################################################

regionize_ <- function(regions, d) {
  r <- list()
  for (j in 1:nrow(regions)) {
    region_ <- regions[j,]
    r[[j]] <- d %>% filter(d$chromosome == region_$chromosome & 
                             region_$start_location < d$position &
                             d$position <= region_$end_location) %>%
      mutate(region = region_$region)
    
  }
  return(do.call(rbind, r))
}

complete_r_ <- function(r_, result, trait) {
  if (!(result %in% unique(r_$results))) {
    r_ <- rbind(r_, data.frame(stringsAsFactors = FALSE, 
                               results=result, n=0, trait=trait))
  }
  r_
}

regionize_ <- function(regions, d) {
  r <- list()
  for (j in 1:nrow(regions)) {
    region_ <- regions[j,]
    r[[j]] <- d %>% filter(d$chromosome == region_$chromosome & 
                             region_$start_location < d$position &
                             d$position <= region_$end_location) %>%
      mutate(region = region_$region)
    
  }
  return(do.call(rbind, r))
}


do_summary <- function(gwas, e, m , regions, annotation, pheno) {
  e_ <- e %>% group_by(gene)  %>% arrange(-rcp) %>% slice(1) %>% ungroup %>% select(gene, rcp) %>% mutate(gene = basename(gene))
  es_ <- e_ %>% filter(rcp>0.5) %>% inner_join(annotation, by=c("gene"="id"))
  er_ <- regionize_(regions, es_) %>% group_by(region) %>% summarise(n=n())
  
  ms_ <- m %>% filter(pvalue < 0.05/nrow(m)) %>% inner_join(annotation, by=c("gene"="id"))
  mr_ <- regionize_(regions, ms_) %>% group_by(region) %>% summarise(n=n())
  
  gs_ <- gwas %>% filter(pvalue < 0.05/nrow(gwas)) %>% 
    mutate(chromosome=as.integer(gsub("chr(.*)_(\\d+)_(.*)_(.*)_b38", "\\1", panel_variant_id))) %>%
    mutate(position=as.integer(gsub("chr(.*)_(\\d+)_(.*)_(.*)_b38", "\\2", panel_variant_id)))
  gr_ <- regionize_(regions, gs_) %>% group_by(region) %>% summarise(n=n())
  
  data.frame(trait=pheno,
             regions = regions %>% nrow,
             enlocalized_regions = er_ %>% nrow,
             multixcan_s_regions = mr_ %>% nrow,
             gwas_s_regions = gr_ %>% nrow,
             gwas_s_regions_enlocalized = gr_ %>% inner_join(er_, by="region") %>% nrow,
             gwas_s_regions_multixcan_s = gr_ %>% inner_join(mr_, by="region") %>% nrow,
             gwas_s_regions_enlocalized_multixcan_s = gr_ %>% inner_join(mr_, by="region") %>% inner_join(er_, by="region") %>% nrow,
             multixcan_genes = m %>% nrow,
             multixcan_s_genes = ms_ %>% nrow,
             enloc_genes = e_ %>% nrow,
             enlocalized_genes = es_ %>% nrow,
             enlocalized_multixcan_s_genes = es_ %>% inner_join(ms_, by="gene") %>% nrow,
             variants = gwas %>% nrow,
             s_variants = gs_ %>% nrow)
}

r <- list()
s <- list()
for (i in 1:nrow(gi)) {
  pheno <- gi$trait[i]
  message(pheno)
  message("gwas")
  gwas <- gi$path[i] %>% r_tsv_(col_types = cols_only(panel_variant_id="c", pvalue="d"))  %>% filter(!is.na(pvalue))
  message("enloc eqtl")
  e_e <- get_enloc_result(enloc_eqtl, pheno)
  message("multixcan eqtl")
  m_e <- multixcan_eqtl %>% filter(trait == pheno) %>% .$path %>% r_tsv_(col_types=cols_only(gene="c", pvalue="d")) %>% filter(!is.na(pvalue))
  
  r[[i]] <- do_summary(gwas, e_e, m_e, regions, gencode_annotation, pheno)
  rm(e_e, m_e)
  
  message("enloc sqtl")
  e_s <- get_enloc_result(enloc_sqtl, pheno, type="vanilla")
  message("multixcan sqtl")
  m_s <- multixcan_sqtl %>% filter(trait == pheno) %>% .$path %>% r_tsv_(col_types=cols_only(gene="c", pvalue="d")) %>% filter(!is.na(pvalue))
  
  s[[i]] <- do_summary(gwas, e_s, m_s, regions, intron_annotation, pheno)
  rm(e_s, m_s)
}
r <- do.call(rbind, r)
s <- do.call(rbind, s)

if (!is.null(row)) {
  save_delim(r, fp_(paste0("summary_", row, ".txt")))
  save_delim(s, fps_(paste0("summary_", row, ".txt")))
} else {
  save_delim(r, fp_("summary.txt"))
  save_delim(s, fps_("summary.txt"))
}


message("done")