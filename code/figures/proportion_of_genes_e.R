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

DATA<-"data/summaries"
dir.create(DATA, showWarnings = FALSE, recursive=TRUE)
dp_ <- function(p) file.path(DATA, p)

FOLDER <-"output/proportions"
dir.create(FOLDER, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(FOLDER, p)


###############################################################################

thresholds <- c(0, 0.01, 0.05) %>% c(seq(0.1,0.9,0.1)) %>% sort(decreasing=FALSE)

message("acquiring data")

get_data_e_ <- function(e_tbl) {
  e_ <- get_e_count_for_thresholds(e_tbl, thresholds)
  data.frame(threshold = thresholds, n = e_, prop = e_/e_[1])
}

data_e_eqtl <- read_or_get(dp_("enloc_proportions_eqtl.txt"), function(){
  get_data_e_(enloc_tbl_eqtl_eur)
})

data_e_sqtl <- read_or_get(dp_("enloc_proportions_sqtl.txt"), function(){
  get_data_e_(enloc_tbl_sqtl_eur)
})


proportion_plot_e_ <- function(data, label) {

  p_ <- data %>%
    ggplot(aes(x=threshold, y=prop)) +
    theme_bw(base_size = 25) + paper_base_theme_ +
    geom_point(color = "#23B509") + geom_line(color = "#23B509") +
    xlab("colocalization threshold") +
    scale_x_continuous(breaks=seq(0,1,0.1))

  p_
}

proportion_plot_e_eqtl <- function() {
  proportion_plot_e_(data_e_eqtl, "colocalized genes") +
    ylab('proportion of genes')
}

proportion_plot_e_sqtl <- function() {
  proportion_plot_e_(data_e_sqtl, "colocalized introns") +
    ylab('proportion of introns')
}

message("plotting")

(function(){
  p_ <- proportion_plot_e_eqtl()
  save_plot(p_, fp_("proportion_e_genes.png"), 600, 600)

  p_ <- proportion_plot_e_sqtl()
  save_plot(p_, fp_("proportion_e_introns.png"), 600, 600)

})()

message("done")
