suppressWarnings(library(dplyr))
suppressWarnings(library(readr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(cowplot))
suppressWarnings(library(ggrepel))
#suppressWarnings(library(RCurl))
#suppressWarnings(library(bigrquery))
suppressWarnings(library(RSQLite))

#suppressWarnings(source("_helpers_table_info.R"))
#suppressWarnings(source("_helpers_big_query.R"))
suppressWarnings(source("_helpers.R"))
suppressWarnings(source("_helpers_data.R"))
suppressWarnings(source("_helpers_qq.R"))
suppressWarnings(source("_data_adapter.R"))

RESULT<-"results/count_hits"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

###############################################################################

gtex_metadata <- "data/gtex_metadata.csv" %>% r_csv_
gwas_metadata <- "data/gwas_metadata.txt" %>% r_tsv_
keep_gwas <- gwas_metadata %>% filter(Deflation == 0) %>% .$Tag

sqtl_predixcan <- get_predixcan_results_logic("v8_en_splicing")
sqtl_multixcan <- g_results_logic("data/smp_v8_splicing", "smultixcan_imputed_gwas_gtexelv8_splicing_(.*)_ccn30.txt")
sqtl_enloc <- get_enloc_results_logic("v8_splicing")
eqtl_predixcan <- get_predixcan_results_logic()
eqtl_multixcan <- g_results_logic("data/smp_v8", "smultixcan_imputed_gwas_gtexelv8_(.*)_ccn30.txt")
eqtl_enloc <- get_enloc_results_logic()

gencode <- "data/gencode_v26_stranded.txt" %>% read_tsv %>%
  mutate(chromosome = as.integer(gsub("chr", "", chromosome)))

#predixcan_en_splicing <- get_predixcan_results_logic("v8_en_splicing")

regions <- "data/eur_ld.annot.gz" %>% r_tsv_(col_types = cols_only(chromosome="c", start_location="i", end_location="i", gene_id="c")) %>%
  rename(region=gene_id) %>% mutate(chromosome = as.integer(gsub("chr", "", chromosome)))

#gi <- g_results_logic("results/gwas", "imputed_(.*).txt") %>% filter(trait %in% keep_gwas)
gi <- g_results_logic("results/gwas_signif", "imputed_(.*).txt") %>% filter(trait %in% keep_gwas)

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) { 
  row <- as.integer(args[1])
  cat ("Kept row ", row, "\n")
  gi <- gi[row,]
} else { 
  row <- NULL 
}


d <- list()
r <- list()
j <- 1
k <- 1
for (i in 1:nrow(gi)) {
  message(gi$trait[i])

###############################################################################
  message("gwas")
  gwas <- gi$path[i] %>% r_tsv_(col_types = cols_only(panel_variant_id="c", pvalue="d"))  %>% filter(!is.na(pvalue))
  gwas_signif <- gwas %>% # filter(pvalue < 0.05/nrow(gwas)) %>% #filter(pvalue < 1e-6) %>%
    mutate(chromosome = as.integer(gsub("chr(\\d+)_(\\d+)_(.*)", "\\1", panel_variant_id))) %>%
    mutate(position = as.integer(gsub("chr(\\d+)_(\\d+)_(.*)", "\\2", panel_variant_id)))
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="gwas", n = nrow(gwas_signif))
  j <- j+1
  rm(gwas)
  rm(gwas_signif)

###############################################################################
  message("predixcan sqtl")
  sqtl_predixcan_ <- get_predixcan_result(sqtl_predixcan, gi$trait[i], col_types=cols_only(gene="c", pvalue="d")) %>% filter(!is.na(pvalue))
  sqtl_predixcan_signif <- sqtl_predixcan_ %>% filter(pvalue < 0.05/nrow(sqtl_predixcan_)) %>%
    mutate(chromosome = as.integer(gsub("intron_(\\d+)_(\\d+)_(.*)", "\\1", gene))) %>%
    mutate(position = as.integer(gsub("intron_(\\d+)_(\\d+)_(.*)", "\\2", gene)))
  
  d[[k]] <- sqtl_predixcan_signif %>% group_by(tissue) %>% summarise(n=n()) %>% mutate(trait = gi$trait[i], method="predixcan_sqtl")
  k <- k+1
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_sqtl_all_pairs", n = nrow(sqtl_predixcan_signif))
  j <- j+1
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_sqtl", n = length(unique(sqtl_predixcan_signif$gene)))
  j <- j+1
    
  message("enloc_sqtl")
  sqtl_enloc_ <- if ( nrow(sqtl_enloc %>% filter(trait == gi$trait[i])) != 0) {
    get_enloc_result(sqtl_enloc, gi$trait[i], type="vanilla")
  } else {
    data.frame(gene=character(0), rcp=numeric(0), experiment_rcp=numeric(0), trait=character(0), tissue=character(0))
  }
  
  enloc_introns_ <- sqtl_enloc_ %>% filter(rcp > 0.5) %>% select(gene, tissue)
  enloc_introns <- enloc_introns_ %>% .$gene %>% unique

  rm(sqtl_enloc_)
  
  sqtl_predixcan_signif_enloc <- sqtl_predixcan_signif %>% inner_join(enloc_introns_, by=c("gene", "tissue"))
  
  d[[k]] <- sqtl_predixcan_signif_enloc %>% group_by(tissue) %>% summarise(n=n()) %>% mutate(trait = gi$trait[i], method="predixcan_sqtl_enloc")
  k <- k+1

  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_sqtl_mrcp", n = nrow(sqtl_predixcan_signif %>% filter(gene %in% enloc_introns)))
  j <- j+1
    
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_sqtl_enloc_all_pairs", n = nrow(sqtl_predixcan_signif_enloc))
  j <- j+1
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_sqtl_enloc", n = length(unique(sqtl_predixcan_signif_enloc$gene)))
  j <- j+1
  
  rm(sqtl_predixcan_)
  rm(sqtl_predixcan_signif)
  rm(sqtl_predixcan_signif_enloc)
  
  ###############################################################################  
  message("multixcan_sqtl")
  sqtl_multixcan_signif <- if ( nrow(sqtl_multixcan %>% filter(trait == gi$trait[i])) != 0) {
    sqtl_multixcan_ <- sqtl_multixcan %>% filter(trait == gi$trait[i]) %>% .$path %>% r_tsv_(col_types=cols_only(gene="c", pvalue="d"))%>% filter(!is.na(pvalue))
    sqtl_multixcan_ %>% filter(pvalue < 0.05/nrow(sqtl_multixcan_)) %>%
      mutate(chromosome = as.integer(gsub("intron_(\\d+)_(\\d+)_(.*)", "\\1", gene))) %>%
      mutate(position = as.integer(gsub("intron_(\\d+)_(\\d+)_(.*)", "\\2", gene)))
  } else {
    data.frame(gene=character(0), pvalue=numeric(0), chromosome=integer(0), position=integer(0))
  }
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="multixcan_sqtl", n = nrow(sqtl_multixcan_signif))
  j <- j+1
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="multixcan_sqtl_enloc", 
      n = sqtl_multixcan_signif %>% filter(gene %in% enloc_introns) %>% nrow)
  j <- j+1
  
  rm(sqtl_multixcan_signif)

  ###############################################################################  
  message("predixcan eqtl")
  eqtl_predixcan_ <- get_predixcan_result(eqtl_predixcan, gi$trait[i], col_types=cols_only(gene="c", pvalue="d")) %>% filter(!is.na(pvalue))
  eqtl_predixcan_signif <- eqtl_predixcan_ %>% filter(pvalue < 0.05/nrow(eqtl_predixcan_)) %>%
    inner_join(gencode %>% select(chromosome=chromosome, position=start_location, gene=gene_id), by="gene")
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_eqtl_all_pairs", n = nrow(eqtl_predixcan_signif))
  j <- j+1
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_eqtl", n = length(unique(eqtl_predixcan_signif$gene)))
  j <- j+1
  
  d[[k]] <- eqtl_predixcan_signif %>% group_by(tissue) %>% summarise(n=n()) %>% mutate(trait = gi$trait[i], method="predixcan_eqtl")
  k <- k+1

  ###############################################################################  
  message("enloc_eqtl")
  eqtl_enloc_ <- if ( nrow(eqtl_enloc %>% filter(trait == gi$trait[i])) != 0) {
    get_enloc_result(eqtl_enloc, gi$trait[i])
  } else {
    data.frame(gene=character(0), rcp=numeric(0), experiment_rcp=numeric(0), trait=character(0), tissue=character(0))
  }
  enloc_genes_ <- eqtl_enloc_ %>% filter(rcp > 0.5) %>% select(gene, tissue)
  enloc_genes <- enloc_genes_ %>% .$gene %>% unique
  rm(eqtl_enloc_)
  eqtl_predixcan_signif_enloc <- eqtl_predixcan_signif %>% inner_join(enloc_genes_, by=c("gene", "tissue"))

  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_eqtl_enloc_all_pairs", n = nrow(eqtl_predixcan_signif_enloc))
  j <- j+1
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_eqtl_mrcp", n = nrow(eqtl_predixcan_signif %>% filter(gene %in% enloc_genes)))
  j <- j+1
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="predixcan_eqtl_enloc", n = length(unique(eqtl_predixcan_signif_enloc$gene)))
  j <- j+1
  
  d[[k]] <- eqtl_predixcan_signif_enloc %>% group_by(tissue) %>% summarise(n=n()) %>% mutate(trait = gi$trait[i], method="predixcan_eqtl_enloc")
  k <- k+1
  
  rm(eqtl_predixcan_)
  rm(eqtl_predixcan_signif)
  rm(eqtl_predixcan_signif_enloc)
  
  
  ###############################################################################  
  message("multixcan_eqtl")
  eqtl_multixcan_signif <- if ( nrow(eqtl_multixcan %>% filter(trait == gi$trait[i])) != 0) {
    eqtl_multixcan_ <- eqtl_multixcan %>% filter(trait == gi$trait[i]) %>% .$path %>% r_tsv_(col_types=cols_only(gene="c", pvalue="d"))%>% filter(!is.na(pvalue))
    eqtl_multixcan_ %>% filter(pvalue < 0.05/nrow(eqtl_multixcan_)) %>%
      inner_join(gencode %>% select(chromosome=chromosome, position=start_location, gene=gene_id), by="gene")
  } else {
    data.frame(gene=character(0), pvalue=numeric(0), chromosome=integer(0), position=integer(0))
  }

  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="multixcan_eqtl", n = nrow(eqtl_multixcan_signif))
  j <- j+1
  
  r[[j]] <- data.frame(stringsAsFactors = FALSE, trait = gi$trait[i], method="multixcan_eqtl_enloc", 
                       n = eqtl_multixcan_signif %>% filter(gene %in% enloc_genes) %>% nrow)
  j <- j+1
  
  rm(eqtl_multixcan_signif)

  if(i>1) break
}
d <- do.call(rbind, d)
r <- do.call(rbind, r)

name_ <- if(!is.null(row)) {
  paste0("pred_by_tissue_count_", row, ".txt")
} else {
  "pred_by_tissue_count.txt"
}
d %>% save_delim(fp_(name_))


name_ <- if(!is.null(row)) {
  paste0("count_", row, ".txt")
} else {
  "count.txt"
}
r %>% save_delim(fp_(name_))
message("done")
