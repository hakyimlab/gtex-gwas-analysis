suppressWarnings(library(dplyr))
suppressWarnings(library(readr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(cowplot))
suppressWarnings(library(ggrepel))
suppressWarnings(library(viridis))
suppressWarnings(library(tidyr))
suppressWarnings(library(RSQLite))

suppressWarnings(source("_helpers.R"))
suppressWarnings(source("_helpers_data.R"))
suppressWarnings(source("_helpers_qq.R"))
suppressWarnings(source("_data_adapter.R"))

RESULT<-"results/prelimin"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)

gwas_metadata <- "data/gwas_metadata.txt" %>% r_tsv_

###############################################################################

d <- if (file.exists("data/stats/enloc_multi/summary.txt")) {
  r_tsv_("data/stats/enloc_multi/summary.txt")
} else {
  load_folder("data/stats/enloc_multi", function(x){grepl("summary_(.*).txt",x)})
}

d <- d %>%  inner_join(gwas_metadata %>% select(trait=Tag, abbreviation=new_abbreviation, category=Category), by="trait")

###############################################################################

d_s <- if (file.exists("data/stats/enloc_multi/summary.txt")) {
  r_tsv_("data/stats/enloc_multi_sqtl/summary.txt")
} else {
  load_folder("data/stats/enloc_multi_sqtl", function(x){grepl("summary_(.*).txt",x)})
}

d_s <- d_s %>%  inner_join(gwas_metadata %>% select(trait=Tag, abbreviation=new_abbreviation, category=Category), by="trait")

###############################################################################


v_ <- viridis(3)
scale_enloc <- c("gwas"="darkgray", "enloc"=v_[2])
scale_smultixcan <- c("gwas"="darkgray", "smultixcan"=v_[3], "enlocalized smultixcan"=v_[2])

###############################################################################

(function(){
  order_ <- d %>% arrange(s_variants) %>% .$abbreviation
  d <- d %>% mutate(abbreviation=factor(abbreviation, levels=order_))
  d_s <- d_s %>% mutate(abbreviation=factor(abbreviation, levels=order_))
  
  #############################################################################
  loci_n <- function(d) {
    d %>% select(abbreviation, enloc=gwas_s_regions_enlocalized, gwas=gwas_s_regions) %>%
      ggplot() + theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(data = d %>% select(abbreviation, gwas=gwas_s_regions) %>% mutate(method="gwas"), 
               mapping=aes(abbreviation, gwas, fill=method)) + 
      geom_col(data = d %>% select(abbreviation, enloc=gwas_s_regions_enlocalized) %>% mutate(method="enloc"), 
               mapping=aes(abbreviation, enloc, fill=method)) + 
      scale_fill_manual(values=scale_enloc, breaks=c("gwas", "enloc")) + coord_flip() +
      ylab("loci")
  }
  
  enloc_loci_n <- loci_n(d) + ggtitle("Number of loci with enlocalized genes")
  #enloc_loci_n %>% save_plot(fp_("enloc_loci_n.png"), 1200, 800)

  enloc_loci_sqtl_n <- loci_n(d_s) + ggtitle("Number of loci with enlocalized introns")
  #enloc_loci_sqtl_n %>% save_plot(fp_("sqtl_enloc_loci_n.png"), 1200, 800)
  
  #############################################################################  
  loci_prop <- function(d) {
    d %>% mutate(enloc=ifelse(gwas_s_regions>0, gwas_s_regions_enlocalized/gwas_s_regions, 0)) %>%
      mutate(gwas=ifelse(gwas_s_regions>0,1-enloc,0)) %>%
      select(abbreviation, enloc, gwas) %>% 
      gather("method", "loci", -abbreviation) %>%
      mutate(method=factor(method, levels=c("gwas", "enloc"))) %>% ggplot() + 
      theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, loci, fill=method)) + 
      scale_fill_manual(values=scale_enloc, breaks=c("gwas", "enloc")) + coord_flip() + ylab("loci proportion")
  }
  enloc_loci_prop <- loci_prop(d) + ggtitle("Proportion of loci with enlocalized genes")
  #enloc_loci_prop %>% save_plot(fp_("enloc_loci_prop.png"), 1200, 800)

  enloc_loci_sqtl_prop <- loci_prop(d_s) + ggtitle("Proportion of loci with enlocalized introns")
  #enloc_loci_sqtl_prop %>% save_plot(fp_("sqtl_enloc_loci_prop.png"), 1200, 800)
  
  #############################################################################
  genes_n <- function(d) {
    d %>% ggplot() + 
      theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, gwas, fill=method), d %>% mutate(gwas=0, method="gwas"), show.legend = FALSE) +
      geom_col(aes(abbreviation, enlocalized_genes, fill=method), d%>% mutate(method="enloc"), show.legend = FALSE) +
      scale_fill_manual(values=scale_enloc, breaks=c("gwas", "enloc")) +
      coord_flip() 
  }
  enloc_genes_n <- genes_n(d) + ggtitle("Number of enlocalized genes")
  #enloc_genes_n %>% save_plot(fp_("enloc_genes.png"), 1200, 800)

  enloc_introns_n <- genes_n(d_s) + ggtitle("Number of enlocalized introns")
  #enloc_introns_n %>% save_plot(fp_("enloc_introns.png"), 1200, 800)
  
  #############################################################################
  p_ <- plot_grid(enloc_genes_n, enloc_loci_n, enloc_loci_prop, labels = c("A", "B", "C"), ncol=3)
  cowplot::save_plot(fp_("ENLOC.png"), p_, base_width =16, base_height=9)
  
  p_ <- plot_grid(enloc_introns_n, enloc_loci_sqtl_n, enloc_loci_sqtl_prop, labels = c("A", "B", "C"), ncol=3)
  cowplot::save_plot(fp_("ENLOC_SQTL.png"), p_, base_width =16, base_height=9)
})()

enloc_prop_loci <- function(d) {
  d %>% mutate(enloc=ifelse(gwas_s_regions>0, gwas_s_regions_enlocalized/gwas_s_regions, 0)) %>%
    mutate(gwas=ifelse(gwas_s_regions>0,1-enloc,0)) %>%
    select(abbreviation, enloc, gwas) %>% 
    gather("method", "loci", -abbreviation) %>%
    mutate(method=factor(method, levels=c("gwas", "enloc")))
}

(function(){
  order_ <- d %>% arrange(s_variants) %>% .$abbreviation
  d <- d %>% mutate(abbreviation=factor(abbreviation, levels=order_))
  d_s <- d_s %>% mutate(abbreviation=factor(abbreviation, levels=order_))
  
  #############################################################################
  loci_n <- function(d) {
    d %>% 
      select(abbreviation, smultixcan=gwas_s_regions_multixcan_s, gwas=gwas_s_regions) %>% 
      gather("method", "loci", -abbreviation) %>% ggplot() + 
      theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, gwas_s_regions, fill=method), d %>% mutate(method="gwas")) + 
      geom_col(aes(abbreviation, gwas_s_regions_multixcan_s, fill=method), d %>% mutate(method="smultixcan")) + 
      geom_col(aes(abbreviation, gwas_s_regions_enlocalized_multixcan_s, fill=method), d %>% mutate(method="enlocalized smultixcan")) + 
      scale_fill_manual(values=scale_smultixcan, breaks=c("gwas", "smultixcan", "enlocalized smultixcan")) + coord_flip() + ylab("loci")
  }
  
  smultixcan_loci_n <- loci_n(d) + ggtitle("Number of loci with\n S-MultiXcan significant genes")
  #smultixcan_loci_n %>% save_plot(fp_("smultixcan_loci_n.png"), 1200, 800)
  
  smultixcan_loci_sqtl_n <- loci_n(d_s) + ggtitle("Number of loci with\n S-MultiXcan significant introns")
  #smultixcan_loci_sqtl_n %>% save_plot(fp_("sqtl_smultixcan_loci_n.png"), 1200, 800)
  
  #############################################################################
  loci_prop <- function(d) {
    d %>% 
      mutate(smultixcan=ifelse(gwas_s_regions>0, gwas_s_regions_multixcan_s/gwas_s_regions, 0)) %>%
      mutate(gwas=ifelse(gwas_s_regions>0,1-smultixcan,0)) %>%
      select(abbreviation, smultixcan, gwas) %>% 
      gather("method", "loci", -abbreviation) %>% 
      mutate(method=factor(method, levels=c("gwas", "enlocalized smultixcan", "smultixcan"))) %>%
      ggplot() + theme_bw() +  theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, loci, fill=method)) + #, position="fill") +
      geom_col(aes(abbreviation, e, fill=method), 
               d %>% mutate(e = ifelse(gwas_s_regions>0,gwas_s_regions_enlocalized_multixcan_s/gwas_s_regions, 0), method="enlocalized smultixcan")) +
      scale_fill_manual(values=scale_smultixcan, breaks=c("gwas", "smultixcan", "enlocalized smultixcan")) +
      coord_flip() + ylab("loci proportion")
  }
  
  smultixcan_loci_prop <- loci_prop(d) + ggtitle("Proportion of loci with\nS-MultiXcan significant genes")
  #smultixcan_loci_prop %>% save_plot(fp_("smultixcan_loci_prop.png"), 1200, 800)
  
  smultixcan_loci_sqtl_prop <- loci_prop(d_s) + ggtitle("Proportion of loci with\nS-MultiXcan significant introns")
  #smultixcan_loci_sqtl_prop %>% save_plot(fp_("sqtl_smultixcan_loci_prop.png"), 1200, 800)
  
  #############################################################################
  genes_n <- function(d) {
    d %>% mutate(gwas=0) %>% ggplot() + 
      theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, gwas, fill=method), d %>% mutate(gwas=0, method="gwas"), show.legend = FALSE) +
      geom_col(aes(abbreviation, multixcan_s_genes, fill=method), d %>% mutate(method="smultixcan"), show.legend = FALSE) +
      geom_col(aes(abbreviation, enlocalized_multixcan_s_genes, fill=method), d %>% mutate(method="enlocalized smultixcan"), show.legend = FALSE) +
      scale_fill_manual(values=scale_smultixcan, breaks=c("gwas", "smultixcan", "enlocalized smultixcan")) + coord_flip()
  }
  
  smultixcan_genes_n <- genes_n(d) + xlab("genes") + ggtitle("Number of\n S-MultiXcan significant genes")
  #smultixcan_genes_n %>% save_plot(fp_("smultixcan_genes.png"), 1200, 800)
  
  smultixcan_introns_n <- genes_n(d_s) + xlab("introns") + ggtitle("Number of\n S-MultiXcan significant introns")
  #smultixcan_introns_n %>% save_plot(fp_("smultixcan_introns.png"), 1200, 800)
  
  #############################################################################
  
  p_ <- plot_grid(smultixcan_genes_n, smultixcan_loci_n, smultixcan_loci_prop, labels = c("A", "B", "C"), ncol=3)#, rel_widths=c(1.4,1,1))
  cowplot::save_plot(fp_("MULTIXCAN.png"), p_, base_width =16, base_height=9)
  
  p_ <- plot_grid(smultixcan_introns_n, smultixcan_loci_sqtl_n, smultixcan_loci_sqtl_prop, labels = c("A", "B", "C"), ncol=3)#, rel_widths=c(1.4,1,1))
  cowplot::save_plot(fp_("MULTIXCAN_SQTL.png"), p_, base_width =16, base_height=9)
  
})()

smultixcan_prop_loci <- function(d) {
  d %>% 
    mutate(smultixcan=ifelse(gwas_s_regions>0, gwas_s_regions_multixcan_s/gwas_s_regions, 0)) %>%
    mutate(enlocalized_smultixcan = ifelse(gwas_s_regions>0,gwas_s_regions_enlocalized_multixcan_s/gwas_s_regions, 0)) %>%
    mutate(gwas=ifelse(gwas_s_regions>0,1-smultixcan,0)) %>%
    select(abbreviation, smultixcan, enlocalized_smultixcan, gwas) %>% 
    gather("method", "loci", -abbreviation) 
}

