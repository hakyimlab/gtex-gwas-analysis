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

message("acquiring data")

data_gti_sp <- read_or_get(dp_("genes_with_sp_introns.txt"), function(){
  gt_ <- get_count_genes_with_sp_intron_for_thresholds(gencode_all_annotation_tbl,  intron_gene_mapping_tbl, predixcan_mashr_tbl_sqtl, thresholds)
  g_ <- get_sp_count_for_thresholds(gencode_all_annotation_tbl, predixcan_mashr_tbl_eqtl,  thresholds)
  data.frame(threshold = thresholds, n = gt_, prop = gt_/g_[1], type="from_intron") %>%
    rbind(data.frame(threshold = thresholds, n = g_, prop = g_/g_[1], type="from_gene"))
})

b_eqtl <- get_bonferroni(predixcan_en_tbl_count_eqtl, pheno_whitelist_) %>% .$b %>% min

proportion_plot_gti_sp_ <- function(data, b = NULL) {
  i_ <- data %>% filter(type == "from_intron", threshold == 1) %>% .$n
  g_ <- data %>% filter(type == "from_gene", threshold == 1) %>% .$n

  l1_ <- paste0("genes with significant associations")
  l2_ <- paste0("genes with a significant intron association")

  palette_ <- c()
  palette_[l1_] <-  rgb(42, 155, 204, maxColorValue = 255)
  palette_[l2_] <- rgb(200, 90, 40, maxColorValue = 255)

  p_ <- data %>% mutate(t = ifelse(type == "from_gene", l1_, l2_)) %>%
    mutate(t = factor(t, levels=c(l1_,l2_))) %>%
    ggplot(aes(x=-log10(threshold), y=prop, color=t)) +
    theme_bw(base_size = 25) + paper_base_theme_ +
    theme(legend.position = c(1, 0.8), legend.justification = c("right", "top"), legend.title = element_blank()) +
    annotate('text', x=-log10(1e-15),y=1,
             label=paste0("# genes: ", data %>% filter(type == "from_gene", prop == 1) %>% .$n),
             hjust=0,size=6) +
    # annotate('text', x=-log10(1e-15),y=0.9,
    #          label=paste0("# genes with introns: ", data %>% filter(type == "from_intron", prop == 1) %>% .$n),
    #          hjust=0,size=6) +
    geom_point() + geom_line() +
    xlab(expression(paste(-log[10],'[',italic(P),']',sep=""))) +
    scale_color_manual(values=palette_) +
    scale_x_continuous(breaks=c(0,5,10,20, 30),labels=c(0,5,10,20,30),expand=c(0.01,0))+
    scale_y_continuous(limits = c(NA,1.05), breaks=c(0,0.05,0.25,0.5,0.75,1.00),labels=c(0,0.05,0.25,0.5,0.75,1.00),expand=c(0,0.01))

  if (!is.null(b)) {
    p_ <- p_ + geom_vline(xintercept = -log10(b), linetype=2, color="gray")
  }

  p_
}

message("plotting")

(function(){
  p_ <- proportion_plot_gti_sp_(data_gti_sp, b_eqtl) + ylab('proportion of genes')
  save_plot(p_, fp_("proportion_gti_sp.png"), 600, 600)
})()


message("done")
