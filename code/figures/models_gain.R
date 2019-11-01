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


d_ <- read_or_get(dp_("models_en_vs_mashr.txt"), function(){
  v7 <- get_models_count(v7_prediction_en_models_extra_tbl_eqtl) %>%
    mutate(tissue = ifelse(tissue == "Cells_Transformed_fibroblasts", "Cells_Cultured_fibroblasts", tissue))
  v8 <- get_models_count(prediction_en_models_extra_tbl_eqtl)
  v8 %>% rename(v8 = models) %>% full_join(v7 %>% rename(v7 = models), by="tissue") %>% fill_with_zeros
})



(function(){
  d <- d_ %>%
    gather("release", "models", -tissue) %>%
    inner_join(gtex_metadata %>% select(tissue, name), by="tissue") %>%
    mutate(name = factor(name, order_)) %>%
    mutate(release = factor(release, levels = c("v8", "v7")))

  p_ <- ggplot(d) + theme_bw(base_size = 12) + paper_base_theme_ +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=rel(1.2)),
          legend.title=element_blank(), legend.position = "bottom",
          axis.text.y=element_text(size=rel(1.7)), axis.title.y = element_blank(),
          axis.text.x=element_text(size=rel(1.7))) +
    geom_col(aes(name, models, fill = release), position="dodge") +
    scale_fill_viridis_d() +
    coord_flip()
  save_plot(p_, fp_("SFIG_MODELS_GAIN_B.png"), 1200, 650)
})()

(function(){
  p_ <- d_ %>% inner_join(gtex_metadata %>% select(tissue, color), by="tissue") %>%
    filter(v7 != 0) %>%
    ggplot() + theme_bw(base_size = 12) + paper_base_theme_ +
    theme(plot.title = element_text(hjust = 0.5, face="bold", size=rel(2)),
          legend.title=element_blank(), legend.position = "bottom",
          axis.text.y=element_text(size=rel(1.7)), axis.title.y = element_text(size=rel(1.5), face = "bold"),
          axis.text.x=element_text(size=rel(1.7)), axis.title.x = element_text(size=rel(1.5), face = "bold")) +
    geom_abline(slope = 1, intercept =0) +
    #geom_point(aes(v7, v8, color=color), size = 3, show.legend = FALSE) +
    geom_point(aes(v7, v8), size = 3, show.legend = FALSE) +
    xlab("v7 - Elastic Net") + ylab("v8 - MASHR") +
    scale_x_continuous(limits = c(0,11000), breaks = seq(0, 10000, by = 2000)) +
    scale_y_continuous(limits = c(0,11000), breaks = seq(0, 10000, by = 2000)) +
    ggtitle("Number of models per tissue")
  save_plot(p_, fp_("SFIG_MODELS_GAIN.png"), 600, 600)
})()


d_%>% filter(v7 != 0) %>%  mutate(f = v8/v7-1) %>% .$f %>% median
