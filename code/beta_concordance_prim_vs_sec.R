suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(argparse)
})

dapg_file_pattern <- "../data/dapg_selected_variants/expression/gwas_and_eqtl/DAPG_with_mashr__{tissue}.rds"
metadata <- read_tsv("../data/gwas_metadata.txt") %>% rename(phenotype=Tag)

# enloc_file_pattern <- "../data/enloc/{tissue}_w_{short_phenotype}_enloc_output.txt.gz"

tissue <- "Whole_Blood"
phenotype <- "Astle_et_al_2016_Eosinophil_counts"
short_phenotype <- as.character(metadata[metadata$phenotype==phenotype, "new_Phenotype"])

df <- readRDS(glue::glue(dapg_file_pattern))[[phenotype]]

df_2_eqtl <- df %>%
  group_by(tissue, phenotype, gene_id) %>%
  filter(n() == 2) %>%
  mutate(rank = order(order(abs(eqtl_effect_size), decreasing=TRUE))) %>%
  mutate(rank = paste0("rank", rank)) %>% # rank1 = primary, rank2 = secondary
  mutate(beta_gene=gwas_effect_size/eqtl_effect_size) %>%
  select(tissue, phenotype, gene_id, rank, beta_gene) %>%
  spread(key=rank, value=beta_gene)

# Discard outliers
df_2_eqtl <- df_2_eqtl[df_2_eqtl$rank1 > quantile(df_2_eqtl$rank1, 0.05, na.rm = TRUE) & df_2_eqtl$rank1 < quantile(df_2_eqtl$rank1, 0.95, na.rm = TRUE),]
df_2_eqtl <- df_2_eqtl[df_2_eqtl$rank2 > quantile(df_2_eqtl$rank2, 0.05, na.rm = TRUE) & df_2_eqtl$rank2 < quantile(df_2_eqtl$rank2, 0.95, na.rm = TRUE),]

# df_enloc <- read.table(enloc_file_pattern %>% glue::glue(), header=TRUE) %>% mutate(gene_id = substr(gene_id, 1, 15))

plot_title <- "{gsub(pattern='_', replacement=' ', x=tissue)} - {gsub(pattern='_', replacement=' ', x=short_phenotype)}\nPrimary vs. secondary eQTL" %>% glue::glue()

pp <- ggplot(df_2_eqtl, aes(rank1, rank2))
pp <- pp + geom_point(alpha=.2, size=4)
pp <- pp + theme_bw(base_size=20)
pp <- pp + geom_abline(slope=1,intercept=0)
# pp <- pp + stat_density2d(color = "white", contour = TRUE)
pp <- pp + stat_density2d(color = "black", contour = TRUE)
pp <- pp + ggtitle(plot_title)
pp <- pp + xlab(expression(beta[prim])) + ylab(expression(beta[sec]))
pp <- pp + xlim(c(-0.06,0.06)) + ylim(c(-0.06,0.06))
pp