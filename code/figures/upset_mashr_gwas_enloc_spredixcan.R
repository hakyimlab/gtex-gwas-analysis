suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(viridis))
suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_plot_style.R"))

SUMMARY <- "data/summaries/mashr_regions.txt"

gwas_metadata <- "data/gwas_metadata.txt" %>% read_tsv

RESULT<-"output/g"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

ldblockinfo = read_tsv(SUMMARY) %>% rename(trait=phenotype, n=count, results=method)

## spread data
tempo <- ldblockinfo %>% mutate(flag = n >= 1) %>% rename(results_type=results) %>%
  select(region, trait, results_type, flag) %>% spread(key=unique(results_type),value=flag)
## set NAs to 0
tempo[is.na(tempo)] <- 0
tempo$gwas <- 1

tempo <- tempo %>%
  mutate(spredixcan_eqtl_enloc = ifelse(spredixcan_eqtl & enloc_eqtl, 1, 0)) %>%
  mutate(spredixcan_sqtl_enloc = ifelse(spredixcan_sqtl & enloc_sqtl, 1, 0)) %>%
  mutate(smultixcan_eqtl_enloc = ifelse(smultixcan_eqtl & enloc_eqtl, 1, 0)) %>%
  mutate(smultixcan_sqtl_enloc = ifelse(smultixcan_sqtl & enloc_sqtl, 1, 0))

#color_pal <- viridis(3)
color_pal <-c('splicing' = rgb(200, 90, 40, maxColorValue = 255),
              'expression' = rgb(42, 155, 204, maxColorValue = 255))
eqtl_ <- (function(){
  d_ <- ldblockinfo %>% select(region,results,trait) %>% unique %>%
    filter(results %in% c("gwas", "spredixcan_eqtl", "enloc_eqtl")) %>% mutate(flag = 1) %>%
    spread(key=unique(results),value=flag)
  d_[is.na(d_)] <- 0
  d_ %>% rename(enloc=enloc_eqtl, spredixcan=spredixcan_eqtl) %>% select(gwas, enloc, spredixcan)
})()

svg(fp_("upset_mashr_spe_expression.svg"))
upset(eqtl_ %>% as.data.frame(),
      sets= c("gwas", "spredixcan", "enloc"), keep.order = TRUE, text.scale = 3,
      main.bar.color = color_pal['expression'], sets.x.label = "Loci", mainbar.y.label = "Shared Loci", mainbar.y.max = 3000)
dev.off()

png(fp_("SFIG_UPSET_MASHR_SPE_EXPRESSION.png"), 800, 600)
upset(eqtl_ %>% as.data.frame(),
      sets= c("gwas", "spredixcan", "enloc"), keep.order = TRUE, text.scale = 3,
      main.bar.color = color_pal['expression'], sets.x.label = "Loci", mainbar.y.label = "Shared Loci", mainbar.y.max = 3000)
dev.off()

sqtl_ <- (function(){
  d_ <- ldblockinfo %>% select(region,results,trait) %>% unique %>%
    filter(results %in% c("gwas", "spredixcan_sqtl", "enloc_sqtl")) %>% mutate(flag = 1) %>%
    spread(key=unique(results),value=flag)
  d_[is.na(d_)] <- 0
  d_ %>% rename(enloc=enloc_sqtl, spredixcan=spredixcan_sqtl) %>% select(gwas, enloc, spredixcan)
})()
svg(fp_("upset_mashr_spe_splicing.svg"))
upset(sqtl_ %>% as.data.frame(),
      sets= c("gwas", "spredixcan", "enloc"), keep.order = TRUE, text.scale = 3,
      main.bar.color = color_pal['splicing'], sets.x.label = "Loci", mainbar.y.label = "Shared Loci", mainbar.y.max = 3000)
dev.off()

png(fp_("SFIG_UPSET_MASHR_SPE_SPLICING.png"), 800, 600)
upset(sqtl_ %>% as.data.frame(),
      sets= c("gwas", "spredixcan", "enloc"), keep.order = TRUE, text.scale = 3,
      main.bar.color = color_pal['splicing'], sets.x.label = "Loci", mainbar.y.label = "Shared Loci", mainbar.y.max = 3000)
dev.off()
