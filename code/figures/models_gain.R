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

k <- get_gtex_tissue_metadata()

gtex_metadata <- get_gtex_metadata()
order_ <- gtex_metadata %>% arrange(-v8_eur) %>% .$tissue

v7 <- get_models_count(v7_prediction_en_models_extra_tbl_eqtl) %>%
  mutate(tissue = ifelse(tissue == "Cells_Transformed_fibroblasts", "Cells_Cultured_fibroblasts", tissue))
v8 <- get_models_count(prediction_en_models_extra_tbl_eqtl)


d <- v8 %>% rename(v8 = models) %>% full_join(v7 %>% rename(v7 = models), by="tissue") %>% fill_with_zeros %>%
  gather("release", "models", -tissue) %>% mutate(tissue = factor(tissue, order_)) %>%
  mutate(release = factor(release, levels = c("v8", "v7")))

p_ <- ggplot(d) + theme_bw(base_size = 12) + paper_base_theme_ +
  geom_col(aes(tissue, models, fill = release), position="dodge") +
  scale_fill_viridis_d() +
  coord_flip()

save_plot(p_, fp_("SFIG_MODELS_GAIN.png"), 1200, 600)

d %>% spread(key="release", value="models") %>% mutate(f = v8/v7-1) %>% .$f %>% median
