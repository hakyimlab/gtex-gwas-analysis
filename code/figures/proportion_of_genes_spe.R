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


###############################################################################

thresholds <- c(0.05, 0.05/49) %>% c(10^(-1*(0:30))) %>% sort(decreasing=TRUE)

get_data_spe_ <- function(sp_tbl, e_tbl) {
  sp_e0 <- get_spe_count_for_thresholds(sp_tbl, e_tbl, thresholds, NULL)
  sp_e0.5 <- get_spe_count_for_thresholds(sp_tbl, e_tbl, thresholds, 0.5)
  sp_e0.1 <- get_spe_count_for_thresholds(sp_tbl, e_tbl, thresholds, 0.1)
  sp_e0.01 <- get_spe_count_for_thresholds(sp_tbl, e_tbl, thresholds, 0.01)

  data.frame(threshold = thresholds, n = sp_e0, prop_a = sp_e0/sp_e0[1], prop_r = sp_e0/sp_e0[1], enloc = 0) %>%
    rbind(data.frame(threshold = thresholds, n = sp_e0.5, prop_a = sp_e0.5/sp_e0[1], prop_r = sp_e0.5/sp_e0.5[1], enloc = 0.5)) %>%
    rbind(data.frame(threshold = thresholds, n = sp_e0.1, prop_a = sp_e0.1/sp_e0[1], prop_r = sp_e0.1/sp_e0.1[1], enloc = 0.1)) %>%
    rbind(data.frame(threshold = thresholds, n = sp_e0.01, prop_a = sp_e0.01/sp_e0[1], prop_r = sp_e0.01/sp_e0.01[1], enloc = 0.01))
}

data_spe_eqtl <- read_or_get(dp_("spredixcan_proportions_eqtl.txt"), function(){
  get_data_spe_(predixcan_mashr_tbl_eqtl, enloc_tbl_eqtl_eur)
})

data_spe_sqtl <- read_or_get(dp_("spredixcan_proportions_sqtl.txt"), function(){
  get_data_spe_(predixcan_mashr_tbl_sqtl, enloc_tbl_sqtl_eur)
})

b_eqtl <- get_bonferroni(predixcan_en_tbl_count_eqtl, pheno_whitelist_) %>% .$b %>% min
b_sqtl <- get_bonferroni(predixcan_en_tbl_count_sqtl, pheno_whitelist_) %>% .$b %>% min
################################################################################

build_decor_ <- function(data) {
  k_ <- data %>% filter(enloc == 0)
  p_at_0.05 <- approx(k_$prop_a, k_$threshold, 0.05)$y
  data.frame(x=c(1,    p_at_0.05, p_at_0.05),
             y=c(0.05,  0.05,      0))
}


proportion_plot_spe_ <- function(data, enloc_threshold, label, draw_0.05 = FALSE, b = NULL) {
  d_ <- data %>% filter(enloc == 0 | enloc == enloc_threshold)

  n_ <- data %>% filter(enloc == 0, threshold == 1) %>% .$n

  l1_ <- paste0("significant ", label)
  l2_ <- paste0("significant ", label, " with rcp>", enloc_threshold)

  palette_ <- c()
  palette_[l1_] = "#6209B5"
  palette_[l2_] <- "#23B509"

  p_ <- d_ %>% mutate(t = ifelse(enloc == 0.0, l1_, l2_)) %>%
    mutate(t = factor(t, levels=c(l1_,l2_))) %>%
    ggplot(aes(x=-log10(threshold), y=prop_a, color=t)) +
    theme_bw(base_size = 25) + paper_base_theme_ +
    theme(legend.position = c(1, 0.6), legend.justification = c("right", "top"), legend.title = element_blank()) +
    annotate('text',x=-log10(1e-20),y=1, label=paste0(n_, ' ', label),hjust=0,size=8)+
    geom_point() + geom_line() +
    xlab(expression(paste(-log[10],'[',italic(P),']',sep=""))) +
    scale_color_manual(values=palette_)

  if (!is.null(b)) {
    p_ <- p_ + geom_vline(xintercept = -log10(b), linetype=2, color="gray")
  }

  if (draw_0.05) {
    decor_ <- build_decor_(d_)
    p_at_0.05 <- decor_$x[2]
    p_ <- p_ +
      scale_x_continuous(breaks=c(0,5,10,20,-log10(p_at_0.05),30),labels=c(0,5,10,20,round(-log10(p_at_0.05),digits=2),30),expand=c(0.01,0))+
      scale_y_continuous(limits = c(NA,1.05), breaks=c(0,0.05,0.25,0.5,0.75,1.00),labels=c(0,0.05,0.25,0.5,0.75,1.00),expand=c(0,0.01)) +
    geom_line(aes(x=-log10(x),y=y), decor_, show.legend = FALSE, linetype=2, color="gray")
  } else {
    p_ <- p_ +
      scale_x_continuous(breaks=c(0,5,10,20, 30),labels=c(0,5,10,20,30),expand=c(0.01,0))+
      scale_y_continuous(limits = c(NA,1.05), breaks=c(0,0.05,0.25,0.5,0.75,1.00),labels=c(0,0.05,0.25,0.5,0.75,1.00),expand=c(0,0.01))

  }
  p_
}

proportion_plot_spe_eqtl <- function(enloc_threshold) {
  proportion_plot_spe_(data_spe_eqtl, enloc_threshold, "genes", b = b_eqtl) +
    ylab('proportion of genes')
}

proportion_plot_spe_sqtl <- function(enloc_threshold) {
  proportion_plot_spe_(data_spe_sqtl, enloc_threshold, "introns", b = b_sqtl) +
    ylab('proportion of introns')
}

(function(){
  p_ <- proportion_plot_spe_eqtl(0.5)
  save_plot(p_, fp_("proportion_spe_genes_0.5.png"), 600, 600)

  p_ <- proportion_plot_spe_eqtl(0.1)
  save_plot(p_, fp_("proportion_spe_genes_0.1.png"), 600, 600)

  p_ <- proportion_plot_spe_eqtl(0.01)
  save_plot(p_, fp_("proportion_spe_genes_0.01.png"), 600, 600)

  p_ <- proportion_plot_spe_sqtl(0.5)
  save_plot(p_, fp_("proportion_spe_introns_0.5.png"), 600, 600)

  p_ <- proportion_plot_spe_sqtl(0.1)
  save_plot(p_, fp_("proportion_spe_introns_0.1.png"), 600, 600)

  p_ <- proportion_plot_spe_sqtl(0.01)
  save_plot(p_, fp_("proportion_spe_introns_0.01.png"), 600, 600)
})()


