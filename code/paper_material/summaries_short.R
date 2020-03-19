suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_plot_style.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

options(gargle_oauth_email = TRUE)
suppressWarnings(source("code/helpers/_helpers_big_query_getters.R"))

RESULT<-"output"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

DATA<-"data/summaries"
dir.create(DATA, showWarnings = FALSE, recursive=TRUE)
dp_ <- function(p) file.path(DATA, p)

k <- get_gtex_tissue_metadata()

gtex_metadata <- get_gtex_metadata()
order_ <- gtex_metadata %>% arrange(-v8_eur) %>% .$name

#########################################################################################
# How many more modes we had, compared to gtex v7
(function(){
  d_ <- suppressMessages((function(){
    v7 <- get_models_count(v7_prediction_en_models_extra_tbl_eqtl) %>%
      mutate(tissue = ifelse(tissue == "Cells_Transformed_fibroblasts", "Cells_Cultured_fibroblasts", tissue))
    v8 <- get_models_count(prediction_mashr_models_extra_tbl_eqtl)
    v8 %>% rename(v8 = models) %>% full_join(v7 %>% rename(v7 = models), by="tissue") %>% fill_with_zeros
  })())

  message("models gain: ", d_%>% filter(v7 != 0) %>%  mutate(f = v8/v7-1) %>% .$f %>% median)
})()

#########################################################################################
# how many of the tested genes are significant
(function(){
  sgenes <- suppressMessages(get_spredixcan_significant_genes(predixcan_mashr_tbl_eqtl, predixcan_mashr_tbl_count_eqtl))
  genes <- suppressMessages(get_spredixcan_genes(predixcan_mashr_tbl_eqtl))
  gene_tissue_pairs <- suppressMessages(get_spredixcan_gene_tissue_pairs(predixcan_mashr_tbl_eqtl))
  message("Unique genes in spredixcan mashr eQTL:", genes)
  message("Unique significant genes in spredixcan mashr eQTL:", sgenes)
  message("Proportion", sgenes/genes)
})()


(function(){
  mashr_gene_tissue_models <- suppressMessages(get_models_gene_tissue_pairs(prediction_mashr_models_extra_tbl_eqtl))
  mashr_intron_tissue_models <- suppressMessages(get_models_gene_tissue_pairs(prediction_mashr_models_extra_tbl_sqtl))
  en_gene_tissue_models <- suppressMessages(get_models_gene_tissue_pairs(prediction_en_models_extra_tbl_eqtl))
  en_intron_tissue_models <- suppressMessages(get_models_gene_tissue_pairs(prediction_en_models_extra_tbl_sqtl))
  data.frame(models=c("en", "en", "mashr", "mashr"), data=c("expression", "splicing", "expression", "splicing"),
             n=c(en_gene_tissue_models, en_intron_tissue_models, mashr_gene_tissue_models, mashr_intron_tissue_models)) %>%
    mutate(b = 0.05/n)

})()

########################################################################################
#numbers of models
(function(){
  model_count <- suppressMessages(get_models_count(prediction_mashr_models_extra_tbl_eqtl))
  message("models in testis:", model_count %>% filter(tissue == "Testis") %>% .$models)
})
