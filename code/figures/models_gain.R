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
order_ <- gtex_metadata %>% arrange(-v8_eur) %>% .$name

v7 <- get_models_count(v7_prediction_en_models_extra_tbl_eqtl) %>%
  mutate(tissue = ifelse(tissue == "Cells_Transformed_fibroblasts", "Cells_Cultured_fibroblasts", tissue))
v8 <- get_models_count(prediction_en_models_extra_tbl_eqtl)


d <- v8 %>% rename(v8 = models) %>% full_join(v7 %>% rename(v7 = models), by="tissue") %>% fill_with_zeros %>%
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

save_plot(p_, fp_("SFIG_MODELS_GAIN.png"), 1200, 650)

d %>% spread(key="release", value="models") %>% filter(v7 != 0) %>%  mutate(f = v8/v7-1) %>% .$f %>% median
