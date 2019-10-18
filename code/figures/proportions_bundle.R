suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(ggplot2))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_plot_style.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

options(gargle_oauth_email = TRUE)
suppressWarnings(source("code/helpers/_helpers_big_query_getters.R"))

DATA<-"data/summaries/proportions"
dir.create(DATA, showWarnings = FALSE, recursive=TRUE)
dp_ <- function(p) file.path(DATA, p)

FOLDER <-"output/proportions"
dir.create(FOLDER, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(FOLDER, p)

s_thresholds <- c(0.05, 0.05/49) %>% c(10^(-1*(0:30))) %>% sort(decreasing=TRUE)
c_thresholds <- c(0, 0.01, 0.05) %>% c(seq(0.1,0.9,0.1)) %>% sort(decreasing=FALSE)

palette_ <- c( "from_gene" =  rgb(42, 155, 204, maxColorValue = 255),
               "from_intron" = rgb(200, 90, 40, maxColorValue = 255))

label_for_type_ <- function(n, type_, x) {
  label_ <- if (type_ == "from_gene") {
    "# genes: {n}"
  } else {
    "# genes with\navailable introns: {n}"
  }
  label_ <- glue::glue(label_, n =n)
  annotate('text', x=x, y=0.9, label=paste0(label_), hjust=0,size=6)
}

###############################################################################

message("acquiring e data")
data_gti_e <- read_or_get(dp_("genes_with_enlocalized_proportions.txt"), function(){
  gt_ <- get_count_genes_with_e_intron_for_thresholds(gencode_all_annotation_tbl,  intron_gene_mapping_tbl, enloc_tbl_sqtl_eur, c_thresholds)
  g_ <- get_e_count_for_thresholds(enloc_tbl_eqtl_eur,  c_thresholds)
  data.frame(enloc_threshold = c_thresholds, n = gt_, prop = gt_/g_[1], type="from_intron") %>%
    rbind(data.frame(enloc_threshold = c_thresholds, n = g_, prop = g_/g_[1], type="from_gene"))
})

proportion_plot_e_ <- function(data, type_, label_="") {
  d_ <- data %>% filter(type == type_)
  c_ <- palette_[type_]

  p_ <- d_ %>%  ggplot(aes(x=enloc_threshold, y=prop)) +
    theme_bw(base_size = 25) + paper_base_theme_ +
    label_for_type_(n = d_ %>% filter(enloc_threshold == 0) %>% .$n, type_, 0.5) +
    xlab('rcp threshold') +
    ylab('proportion of genes') +
    geom_point(color=c_, size=2) + geom_line(color=c_, size=1) +
    scale_x_continuous(breaks= seq(0,1, 0.1))+
    scale_y_continuous(limits = c(NA,1.05), breaks=c(0,0.05,0.25,0.5,0.75,1.00),labels=c(0,0.05,0.25,0.5,0.75,1.00),expand=c(0,0.01))

  p_
}

message("plotting e")
(function(){
  p_ <- proportion_plot_e_(data_gti_e, "from_gene", label_) + ggtitle("Colocalized genes\nusing expression signal")
  save_plot(p_, fp_("SFIG_PROPORTIONS_GENES_ENLOC_FROM_EQTL.png"), 600, 600)

  p_ <- proportion_plot_e_(data_gti_e, "from_intron", label_)  + ggtitle("Colocalized genes\nusing splicing signal")
  save_plot(p_, fp_("SFIG_PROPORTIONS_GENES_ENLOC_FROM_SQTL.png"), 600, 600)
})()

###############################################################################

data_spe <- read_or_get(dp_("genes_with_spe_proportions.txt"), function(){
  sp_e0 <- get_spe_count_for_thresholds(sp_tbl, e_tbl, thresholds, NULL)
  sp_e0.5 <- get_spe_count_for_thresholds(sp_tbl, e_tbl, thresholds, 0.5)
  sp_e0.1 <- get_spe_count_for_thresholds(sp_tbl, e_tbl, thresholds, 0.1)

  sp_s0 <- get_count_genes_with_spe_intron_for_thresholds(gencode_all_annotation_tbl,  intron_gene_mapping_tbl, predixcan_mashr_tbl_eqtl, enloc_tbl_sqtl_eur, s_thresholds, NULL)
  sp_s0.5 <- get_count_genes_with_spe_intron_for_thresholds(gencode_all_annotation_tbl,  intron_gene_mapping_tbl, predixcan_mashr_tbl_eqtl, enloc_tbl_sqtl_eur, s_thresholds, 0.5)
  sp_s0.1 <- get_count_genes_with_spe_intron_for_thresholds(gencode_all_annotation_tbl,  intron_gene_mapping_tbl, predixcan_mashr_tbl_eqtl, enloc_tbl_sqtl_eur, s_thresholds, 0.1)

  data.frame(threshold = s_thresholds, n = sp_e0, prop_a = sp_e0/sp_e0[1], prop_r = sp_e0/sp_e0[1], enloc_threshold = 0) %>%
    rbind(data.frame(threshold = s_thresholds, n = sp_e0.5, prop_a = sp_e0.5/sp_e0[1], prop_r = sp_e0.5/sp_e0.5[1], enloc_threshold = 0.5)) %>%
    rbind(data.frame(threshold = s_thresholds, n = sp_e0.1, prop_a = sp_e0.1/sp_e0[1], prop_r = sp_e0.1/sp_e0.1[1], enloc_threshold = 0.1)) %>%
    rbind(data.frame(threshold = s_thresholds, n = sp_s0, prop_a = sp_s0/sp_e0[1], prop_r = sp_s0/sp_se0[1], enloc_threshold = 0)) %>%
    rbind(data.frame(threshold = s_thresholds, n = sp_s0.1, prop_a = sp_s0.1/sp_e0[1], prop_r = sp_s0.1/sp_s0.1[1], enloc_threshold = 0.1)) %>%
    rbind(data.frame(threshold = s_thresholds, n = sp_s0.5, prop_a = sp_s0.5/sp_e0[1], prop_r = sp_s0.5/sp_s0.5[1], enloc_threshold = 0.5))
})

b_eqtl <- get_bonferroni(predixcan_en_tbl_count_eqtl, pheno_whitelist_) %>% .$b %>% min

message("Plotting spe")

build_decor <- function(data) {
  data %>% filter(threshold == 1e-10) %>%
    mutate(text = ifelse(enloc_threshold == 0, "all genes",
                  ifelse(enloc_threshold == 0.1, "rcp > 0.1", "rcp > 0.5"))) %>%
    mutate(text = factor(text, levels = c("all genes", "rcp > 0.1", "rcp > 0.5"))) %>%
    arrange(text) %>%
    mutate(x = 1e-15, x_start =1e-13, x_end = 1e-10) %>%
    mutate(y = c(0.3, 0.25, 0.20), y_end = prop_a)
}

proportion_plot_spe_ <- function(data, type_, b = NULL) {
  d_ <- data %>% filter(type == type_)

  c_ <- palette_[type_]

  decor_ <- build_decor(d_)

  p_ <- d_ %>%
    ggplot(aes(x=-log10(threshold), y=prop_a)) +
    theme_bw(base_size = 25) + paper_base_theme_ +
    geom_point(aes(group = enloc_threshold), show.legend = FALSE, size = 2, color=c_) +
    geom_line(aes(group = enloc_threshold), show.legend = FALSE, size=1.5, color=c_) +
    geom_text(aes(-log10(x), y, label=text), decor_, size = 5) +
    geom_segment(aes(x=-log10(x_start), y=y, xend=-log10(x_end), yend=y_end), decor_, arrow=arrow(length=unit(0.02,'npc'))) +
    label_for_type_(d_ %>% filter(type ==type_) %>% filter(threshold == 1, enloc_threshold ==0) %>% .$n, type_, -log10(1e-15)) +
    xlab(expression(paste(-log[10],'[',italic(P),']',sep=""))) +
    ylab("Proportion of genes") +
    scale_x_continuous(breaks=c(0,5,10,20, 30),labels=c(0,5,10,20,30),expand=c(0.01,0))+
    scale_y_continuous(limits = c(NA,1.05), breaks=c(0,0.05,0.25,0.5,0.75,1.00),labels=c(0,0.05,0.25,0.5,0.75,1.00),expand=c(0,0.01))

    scale_color_manual(values=palette_)

  if (!is.null(b)) {
    p_ <- p_ + geom_vline(xintercept = -log10(b), linetype=2, color="gray")
  }

  p_
}

proportion_plot_spe_(data_spe, "from_gene", b_eqtl)

(function() {

  p_ <- proportion_plot_spe_(data_spe, "from_gene", b_eqtl) + ggtitle("Significant genes\nusing expression association")
  save_plot(p_, fp_("SFIG_PROPORTIONS_GENES_SPE_FROM_EQTL.png"), 600, 600)

  p_ <- proportion_plot_spe_(data_spe, "from_intron", b_eqtl) + ggtitle("Significant genes\nusing splicing association")
  save_plot(p_, fp_("SFIG_PROPORTIONS_GENES_SPE_FROM_SQTL.png"), 600, 600)
})()

message("colocalized prop eqtl: ", data_gti_e %>% filter(type == "from_gene", enloc_threshold == 0.5) %>% .$prop)
message("colocalized prop sqtl: ", data_gti_e %>% filter(type == "from_intron", enloc_threshold == 0.5) %>% .$

d_ <- data_spe %>% filter(type == "from_gene", enloc_threshold == 0)
message("gene e0 ", approx(d_$threshold, d_$prop_a, b_eqtl)$y)

d_ <- data_spe %>% filter(type == "from_gene", enloc_threshold == 0.5)
message("gene e0.5 ", approx(d_$threshold, d_$prop_a, b_eqtl)$y)

d_ <- data_spe %>% filter(type == "from_intron", enloc_threshold == 0)
message("intron e0 ", approx(d_$threshold, d_$prop_a, b_eqtl)$y)

d_ <- data_spe %>% filter(type == "from_intron", enloc_threshold == 0.5)
message("intron s0.5 ", approx(d_$threshold, d_$prop_a, b_eqtl)$y)
