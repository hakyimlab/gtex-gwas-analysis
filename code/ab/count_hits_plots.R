suppressWarnings(library(dplyr))
suppressWarnings(library(tidyr))
suppressWarnings(library(readr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(cowplot))
suppressWarnings(library(ggrepel))
#suppressWarnings(library(RCurl))
#suppressWarnings(library(bigrquery))
suppressWarnings(library(RSQLite))

#suppressWarnings(source("_helpers_table_info.R"))
#suppressWarnings(source("_helpers_big_query.R"))
suppressWarnings(source("_helpers.R"))
suppressWarnings(source("_helpers_data.R"))
suppressWarnings(source("_data_adapter.R"))

###############################################################################

RESULT<-"results/count_hits"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

theme_ <-function() {
  theme_bw() +
    theme(plot.title = element_text(hjust=0.5, face="bold", size=27),
          plot.subtitle = element_text(hjust=0.5, face="italic", size=25),
          axis.title = element_text(size=25),
          axis.text = element_text(size=20),
          axis.text.x = element_text(size=15, angle=90, hjust=1),
          legend.text = element_text(size = 15),
          legend.position="bottom")
}

###############################################################################

gwas_metadata <- "data/gwas_metadata.txt" %>% r_tsv_ %>% filter(Deflation == 0)
trait_key <- gwas_metadata %>% select(trait=Tag, abbreviation=new_abbreviation, color=color, category=Category)

counts <- if (TRUE) {
    k <- (function() {
      FOLDER <- RESULT
      files <- file.path(FOLDER, list.files(FOLDER))
      files <- files[grepl("/count_(\\d+).txt", files, perl = TRUE)]
      r <- list()
      for (i in 1:length(files)) {
        r[[i]] <- r_tsv_(files[i])
      }
      do.call(rbind, r)
    })() %>% inner_join(trait_key, by="trait") 
    k %>% save_delim(fp_("count.txt"))
    k
  } else {
    fp_("count.txt") %>% r_tsv_ %>% inner_join(trait_key, by="trait")
  }
counts <- counts %>% mutate(method = factor(method, levels = 
    c("gwas", "predixcan_sqtl", "predixcan_eqtl", "predixcan_sqtl_enloc", "predixcan_eqtl_enloc",
        "predixcan_sqtl_all_pairs", "predixcan_eqtl_all_pairs",  "predixcan_sqtl_enloc_all_pairs", "predixcan_eqtl_enloc_all_pairs",
        "predixcan_sqtl_mrcp", "predixcan_eqtl_mrcp",
        "multixcan_sqtl", "multixcan_eqtl", "multixcan_sqtl_enloc", "multixcan_eqtl_enloc")))
order_ <- counts %>% filter(method == "gwas") %>% arrange(n) %>% .$abbreviation
counts <- counts %>% mutate(abbreviation = factor(abbreviation, levels=order_))

pred_by_tissue_counts <- if (TRUE) {
  k <- (function() {
    FOLDER <- RESULT
    files <- file.path(FOLDER, list.files(FOLDER))
    files <- files[grepl("/pred_by_tissue_count_(\\d+).txt", files, perl = TRUE)]
    r <- list()
    for (i in 1:length(files)) {
      r[[i]] <- r_tsv_(files[i])
    }
    do.call(rbind, r)
  })() %>% inner_join(trait_key, by="trait") 
  k %>% save_delim(fp_("pred_by_tissue_count.txt"))
  k
} else {
  fp_("pred_by_tissue_count.txt") %>% r_tsv_ %>% inner_join(trait_key, by="trait")
}
pred_by_tissue_counts <- pred_by_tissue_counts %>% mutate(abbreviation = factor(abbreviation, levels = order_))

###############################################################################

palette_ <- trait_key %>% select(category, color) %>% unique %>% .$color
categories_ <- trait_key %>% select(category, color) %>% unique %>% .$category
names(palette_) <- categories_

###############################################################################

(function(){
  d_ <- counts %>% filter(method == "gwas") %>% select(trait=trait, gwas=n) %>%
    inner_join(counts %>% filter(method !="gwas"), by="trait")
  
  p_ <- d_ %>% filter(method %in% c("predixcan_eqtl_all_pairs", "predixcan_eqtl_enloc_all_pairs", "predixcan_eqtl_mrcp")) %>% ggplot() + theme_() +
    geom_point(aes(x=gwas, y=n, color=category), size=4) + facet_wrap(~method, ncol=3, scales = "free_y") +
    scale_color_manual(values = palette_) + xlab("gwas detections") + ylab("predixcan detections") +
    ggtitle("Predixcan (gene,tissue) associations")
  
  p_ %>% save_plot(fp_("a_count_hits_predixcan_pairs.png"), 600, 1300)
  
  p_ <- d_ %>% filter(method %in% c("predixcan_eqtl", "predixcan_eqtl_enloc")) %>% ggplot() + theme_() +
    geom_point(aes(x=gwas, y=n, color=category), size=4) + facet_wrap(~method, ncol=2, scales = "free_y") +
    scale_color_manual(values = palette_) + xlab("gwas detections") + ylab("predixcan detections") +
    ggtitle("Predixcan associations")
  
  p_ %>% save_plot(fp_("a_count_hits_predixcan.png"), 600, 1000)
  
  p_ <- d_ %>% filter(grepl("multixcan", method) & grepl("eqtl", method)) %>% ggplot() + theme_() +
    geom_point(aes(x=gwas, y=n, color=category), size=4) + facet_wrap(~method, ncol=2, scales = "free_y") +
    scale_color_manual(values = palette_) + xlab("gwas detections") + ylab("multixcan detections")
  
  p_ %>% save_plot(fp_("a_count_hits_multixcan.png"), 600, 1000)
  
  d_ <- counts %>% filter(method == "gwas") %>% select(trait=trait, gwas=n) %>%
    inner_join(pred_by_tissue_counts, by="trait") %>% 
    mutate(method = factor(method, levels = c("predixcan_sqtl", "predixcan_sqtl_enloc", "predixcan_eqtl", "predixcan_eqtl_enloc")))
  
  p_ <- d_ %>% filter(grepl("enloc", method)) %>% ggplot() + theme_() + 
    geom_boxplot(aes(abbreviation, n, color=method)) + coord_flip() + 
    ylab("predixcan associations (across tissues)") + xlab("trait")
  p_ %>% save_plot(fp_("a_count_hits_predixcan_by_tissue.png"), 1600, 800)
  
  p_ <- d_ %>% filter(grepl("predixcan_eqtl_enloc", method)) %>% mutate(gwas = ifelse(gwas>30000, 30000, gwas)) %>%  ggplot() + theme_() + 
    geom_boxplot(aes(gwas, n, group = cut_width(gwas, 1000))) + xlab("gwas associations") + ylab("Predixcan enlocalized associations")
  p_ %>% save_plot(fp_("a_count_hits_predixcan_by_tissue_2.png"), 600, 800)
  
})()