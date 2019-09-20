###############################################################################
# USE MINICONDA
###############################################################################
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))


suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

################################################################################
# PRELIMINARIES
RESULT<-"output/enloc"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)


v_ <- viridis(3)
scale_enloc <- c("gwas"="darkgray", "enloc"=v_[2])
scale_smultixcan <- c("gwas"="darkgray", "smultixcan"=v_[3], "enlocalized smultixcan"=v_[2])

regions <- "data/summaries/regions.txt" %>% r_tsv_ %>% rename(trait=phenotype)
gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as trait, new_abbreviation as abbreviation, Category as category from
                      {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name} WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

pheno_whitelist <- sprintf("(%s)", toString(sprintf("'%s'", gwas_metadata$trait)))

d <- (function(){
  d <- regions
  d <- d %>% spread(key=method, value=count) %>%
    mutate_at(c("gwas", "enloc_eqtl", "smultixcan_eqtl", "enloc_sqtl", "smultixcan_sqtl"), function(x){ifelse(x>0,1,0)})
  d[is.na(d)] <- 0
  d <- d %>% group_by(trait) %>%
    summarise(gwas_s_regions=sum(gwas),
              gwas_s_regions_enlocalized_eqtl=sum(enloc_eqtl),
              gwas_s_regions_smultixcan_eqtl=sum(smultixcan_eqtl),
              gwas_s_regions_smultixcan_enlocalized_eqtl=sum(smultixcan_eqtl*enloc_eqtl),
              gwas_s_regions_enlocalized_sqtl=sum(enloc_sqtl),
              gwas_s_regions_smultixcan_sqtl=sum(smultixcan_sqtl),
              gwas_s_regions_smultixcan_enlocalized_sqtl=sum(smultixcan_sqtl*enloc_sqtl))
  left_ <- gwas_metadata %>% filter(!(trait %in% unique(d$trait))) %>% .$trait %>% sort %>% unlist
  left_ <- data.frame(trait=left_,
                      gwas_s_regions=0,
                      gwas_s_regions_enlocalized_eqtl=0,
                      gwas_s_regions_smultixcan_eqtl=0,
                      gwas_s_regions_smultixcan_enlocalized_eqtl=0,
                      gwas_s_regions_enlocalized_sqtl=0,
                      gwas_s_regions_smultixcan_sqtl=0,
                      gwas_s_regions_smultixcan_enlocalized_sqtl=0)
  rbind(d, left_) %>%
    inner_join(gwas_metadata, by="trait")
})()

order_ <- d %>% arrange(gwas_s_regions) %>% .$abbreviation

d <- d %>% mutate(abbreviation = factor(abbreviation, levels=order_))


################################################################################
#additional data
enloc_genes_ <- (function() {
  query_ <- glue::glue(
    "SELECT molecular_qtl_trait as gene_id, locus_rcp as rcp, phenotype as trait, tissue
             FROM {enloc_tbl_eqtl_eur$dataset_name}.{enloc_tbl_eqtl_eur$table_name} as e
             WHERE locus_rcp>0.5 and e.phenotype in {pheno_whitelist}"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
enloc_genes <- enloc_genes_ %>% group_by(trait, gene_id) %>% arrange(-rcp) %>% slice(1) %>% ungroup %>%
  group_by(trait) %>% summarise(n=n())  %>%  inner_join(gwas_metadata, by="trait") %>%
  mutate(abbreviation = factor(abbreviation, levels=order_))

enloc_introns_ <- (function() {
  query_ <- glue::glue(
    "SELECT molecular_qtl_trait as intron_id, locus_rcp as rcp, phenotype as trait, tissue
             FROM {enloc_tbl_sqtl_eur$dataset_name}.{enloc_tbl_sqtl_eur$table_name} as s
             WHERE locus_rcp>0.5 and s.phenotype in {pheno_whitelist}"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
enloc_introns <- enloc_introns_ %>% group_by(trait, intron_id) %>% arrange(-rcp) %>% slice(1) %>% ungroup %>%
  group_by(trait) %>% summarise(n=n())  %>%  inner_join(gwas_metadata, by="trait") %>%
  mutate(abbreviation = factor(abbreviation, levels=order_))

###############################################################################
# ENLOC Plots
(function(){
  #############################################################################
  genes_n <- function(d) {
    d %>% ggplot() +
      theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, n), show.legend = FALSE, fill=scale_enloc["enloc"]) +
      coord_flip()
  }

  enloc_genes_n <- genes_n(enloc_genes) + ggtitle("Number of colocalized genes")
  enloc_genes_n %>% save_plot(fp_("enloc_genes.png"), 1200, 800)

  enloc_introns_n <- genes_n(enloc_introns) + ggtitle("Number of colocalized introns")
  enloc_introns_n %>% save_plot(fp_("enloc_introns.png"), 1200, 800)

  #############################################################################
  loci_n <- function(d) {
    d %>% ggplot() + theme_bw(base_size = 10) + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, n, fill=method),d %>% select(abbreviation, n=gwas) %>% mutate(method="gwas")) +
      geom_col(aes(abbreviation, n, fill=method), d %>% select(abbreviation, n=enloc) %>% mutate(method="enloc")) +
      scale_fill_manual(values=scale_enloc, breaks=c("gwas", "enloc")) + coord_flip() +
      ylab("loci")
  }

  enloc_loci_n_eqtl <- d %>% select(abbreviation, gwas=gwas_s_regions, enloc=gwas_s_regions_enlocalized_eqtl) %>% loci_n() +
    ggtitle("Number of loci with colocalized genes")#, subtitle="comparison to loci with gwas associations")
  enloc_loci_n_eqtl %>% save_plot(fp_("enloc_loci_eqtl_n.png"), 1200, 800)

  enloc_loci_n_sqtl <- d %>% select(abbreviation, gwas=gwas_s_regions, enloc=gwas_s_regions_enlocalized_sqtl) %>% loci_n() +
    ggtitle("Number of loci with colocalized introns")#, subtitle="comparison to loci with gwas associations")
  enloc_loci_n_sqtl %>% save_plot(fp_("enloc_loci_sqtl_n.png"), 1200, 800)

  #############################################################################
  loci_prop <- function(d) {
    d %>% mutate(enloc=ifelse(gwas>0, enloc/gwas, 0)) %>%
      mutate(gwas=ifelse(gwas>0,1-enloc,0)) %>%
      select(abbreviation, enloc, gwas) %>%
      gather("method", "loci", -abbreviation) %>%
      mutate(method=factor(method, levels=c("gwas", "enloc"))) %>% ggplot() +
      theme_bw(base_size = 10) + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, loci, fill=method)) +
      scale_fill_manual(values=scale_enloc, breaks=c("gwas", "enloc")) + coord_flip() + ylab("loci proportion")
  }
  enloc_loci_prop_eqtl <- d %>% select(abbreviation, gwas=gwas_s_regions, enloc=gwas_s_regions_enlocalized_eqtl) %>% loci_prop() +
    ggtitle("Proportion of loci with colocalized genes")#, subtitle="comparison to loci with gwas associations")
  enloc_loci_prop_eqtl %>% save_plot(fp_("enloc_loci_eqtl_prop.png"), 1200, 800)

  enloc_loci_prop_sqtl <- d %>% select(abbreviation, gwas=gwas_s_regions, enloc=gwas_s_regions_enlocalized_sqtl) %>% loci_prop() +
    ggtitle("Proportion of loci with colocalized introns")#, subtitle="comparison to loci with gwas associations")
  enloc_loci_prop_eqtl %>% save_plot(fp_("enloc_loci_sqtl_prop.png"), 1200, 800)

  #############################################################################
  p_ <- cowplot::plot_grid(enloc_genes_n, enloc_loci_n_eqtl, enloc_loci_prop_eqtl, labels = c("A", "B", "C"), ncol=3)
  cowplot::save_plot(fp_("ENLOC_EQTL.png"), p_, base_width =16, base_height=9)

  p_ <- cowplot::plot_grid(enloc_introns_n, enloc_loci_n_sqtl, enloc_loci_prop_sqtl, labels = c("A", "B", "C"), ncol=3)
  cowplot::save_plot(fp_("ENLOC_SQTL.png"), p_, base_width =16, base_height=9)
})()



################################################################################
#additional data
smultixcan_genes <- (function() {
  query_ <- glue::glue("
SELECT m.phenotype as trait, m.gene, m.pvalue
FROM {multixcan_tbl_eqtl$dataset_name}.{multixcan_tbl_eqtl$table_name} as m
JOIN {multixcan_tbl_count_eqtl$dataset_name}.{multixcan_tbl_count_eqtl$table_name} as m_count
ON m_count.phenotype=m.phenotype
WHERE m.pvalue < m_count.b and m.phenotype in {pheno_whitelist}
")
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
smultixcan_genes <- smultixcan_genes %>%
  left_join(enloc_genes_ %>% select(trait, gene=gene_id) %>% mutate(enloc=1), by=c("trait", "gene")) %>%
  mutate(enloc=ifelse(is.na(enloc),0,1)) %>%
  group_by(trait) %>% summarise(smultixcan=n(), enloc=sum(enloc)) %>%
  inner_join(gwas_metadata, by="trait") %>%
  mutate(abbreviation = factor(abbreviation, levels=order_))

smultixcan_introns <- (function() {
  query_ <- glue::glue("
SELECT m.trait, m.gene as intron, m.pvalue
FROM {multixcan_tbl_sqtl$dataset_name}.{multixcan_tbl_sqtl$table_name} as m
JOIN {multixcan_tbl_count_sqtl$dataset_name}.{multixcan_tbl_count_sqtl$table_name} as m_count
ON m_count.trait=m.trait
WHERE m.pvalue < m_count.b and m.trait in {pheno_whitelist}
")
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
smultixcan_introns <- smultixcan_introns %>%
  left_join(enloc_introns_ %>% select(trait, intron=intron_id) %>% mutate(enloc=1), by=c("trait", "intron")) %>%
  mutate(enloc=ifelse(is.na(enloc),0,1)) %>%
  group_by(trait) %>% summarise(smultixcan=n(), enloc=sum(enloc)) %>%
  inner_join(gwas_metadata, by="trait") %>%
  mutate(abbreviation = factor(abbreviation, levels=order_))

###############################################################################
# SMultiXcan
(function(){
  #############################################################################
  loci_n <- function(d) {
    d %>%  gather("method", "loci", -abbreviation) %>% ggplot() +
      theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, gwas, fill=method), d %>% mutate(method="gwas")) +
      geom_col(aes(abbreviation, smultixcan, fill=method), d %>% mutate(method="smultixcan")) +
      geom_col(aes(abbreviation, smultixcan_enloc, fill=method), d %>% mutate(method="enlocalized smultixcan")) +
      scale_fill_manual(values=scale_smultixcan, breaks=c("gwas", "smultixcan", "enlocalized smultixcan")) + coord_flip() + ylab("loci")
  }

  smultixcan_loci_n <- d %>%
    select(abbreviation, smultixcan=gwas_s_regions_smultixcan_eqtl, gwas=gwas_s_regions, smultixcan_enloc=gwas_s_regions_smultixcan_enlocalized_eqtl)  %>%
    loci_n() + ggtitle("Number of loci with\n S-MultiXcan significant genes") #, subtitle="")
  smultixcan_loci_n %>% save_plot(fp_("smultixcan_loci_n_eqtl.png"), 1200, 800)

  smultixcan_loci_sqtl_n <- d %>%
    select(abbreviation, smultixcan=gwas_s_regions_smultixcan_sqtl, gwas=gwas_s_regions, smultixcan_enloc=gwas_s_regions_smultixcan_enlocalized_sqtl)  %>%
    loci_n() + ggtitle("Number of loci with\n S-MultiXcan significant introns")
  smultixcan_loci_sqtl_n %>% save_plot(fp_("smultixcan_loci_n_sqtl.png"), 1200, 800)

  #############################################################################
  loci_prop <- function(d) {
    d %>%
      mutate(smultixcan=ifelse(gwas>0, smultixcan/gwas, 0)) %>%
      mutate(gwas=ifelse(gwas>0,1-smultixcan,0)) %>%
      select(abbreviation, smultixcan, gwas) %>%
      gather("method", "loci", -abbreviation) %>%
      mutate(method=factor(method, levels=c("gwas", "enlocalized smultixcan", "smultixcan"))) %>%
      ggplot() + theme_bw() +  theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, loci, fill=method)) + #, position="fill") +
      geom_col(aes(abbreviation, e, fill=method),
               d %>% mutate(e = ifelse(gwas>0,smultixcan_enloc/gwas, 0), method="enlocalized smultixcan")) +
      scale_fill_manual(values=scale_smultixcan, breaks=c("gwas", "smultixcan", "enlocalized smultixcan")) +
      coord_flip() + ylab("loci proportion")
  }

  smultixcan_loci_prop <- d %>%
    select(abbreviation, smultixcan=gwas_s_regions_smultixcan_eqtl, gwas=gwas_s_regions, smultixcan_enloc=gwas_s_regions_smultixcan_enlocalized_eqtl) %>%
    loci_prop() + ggtitle("Proportion of loci with\nS-MultiXcan significant genes")
  smultixcan_loci_prop %>% save_plot(fp_("smultixcan_loci_prop_eqtl.png"), 1200, 800)

  smultixcan_loci_sqtl_prop <- d %>%
    select(abbreviation, smultixcan=gwas_s_regions_smultixcan_sqtl, gwas=gwas_s_regions, smultixcan_enloc=gwas_s_regions_smultixcan_enlocalized_sqtl) %>%
    loci_prop() + ggtitle("Proportion of loci with\nS-MultiXcan significant introns")
  smultixcan_loci_sqtl_prop %>% save_plot(fp_("smultixcan_loci_prop_sqtl.png"), 1200, 800)

  #############################################################################
  genes_n <- function(d) {
    d %>% mutate(gwas=0) %>% ggplot() +
      theme_bw() + theme(legend.title=element_blank(), legend.position = "bottom") +
      geom_col(aes(abbreviation, smultixcan, fill=method), d %>% mutate(method="smultixcan")) +
      geom_col(aes(abbreviation, enloc, fill=method), d %>% mutate(method="enlocalized smultixcan")) +
      scale_fill_manual(values=scale_smultixcan, breaks=c("gwas", "smultixcan", "enlocalized smultixcan")) + coord_flip()
  }

  smultixcan_genes_n <- smultixcan_genes %>% genes_n() + xlab("genes") +
    ggtitle("Number of\n S-MultiXcan significant genes")#, subtitle="comparison to colocalized genes")
  smultixcan_genes_n %>% save_plot(fp_("smultixcan_genes.png"), 1200, 800)

  smultixcan_introns_n <- smultixcan_introns %>% genes_n() + xlab("introns") +
    ggtitle("Number of \nS-MultiXcan significant introns")#, subtitle="comparison to colocalized introns")
  smultixcan_introns_n %>% save_plot(fp_("smultixcan_introns.png"), 1200, 800)


  p_ <- cowplot::plot_grid(smultixcan_genes_n, smultixcan_loci_n, smultixcan_loci_prop, labels = c("A", "B", "C"), ncol=3)#, rel_widths=c(1.4,1,1))
  cowplot::save_plot(fp_("MULTIXCAN_EQTL.png"), p_, base_width =16, base_height=9)

  p_ <- cowplot::plot_grid(smultixcan_introns_n, smultixcan_loci_sqtl_n, smultixcan_loci_sqtl_prop, labels = c("A", "B", "C"), ncol=3)#, rel_widths=c(1.4,1,1))
  cowplot::save_plot(fp_("MULTIXCAN_SQTL.png"), p_, base_width =16, base_height=9)
})()

cat("median proportion of gwas loci with enloc genes:\n")
d %>% filter(gwas_s_regions>0) %>% mutate(f = gwas_s_regions_enlocalized_eqtl/gwas_s_regions) %>% .$f %>% median() %>% cat("\n")

cat("median proportion of gwas loci with enloc introns:\n")
d %>% filter(gwas_s_regions>0) %>% mutate(f = gwas_s_regions_enlocalized_sqtl/gwas_s_regions) %>% .$f %>% median() %>% cat("\n")

cat("median proportion of gwas loci with smultixcan genes:\n")
d %>% filter(gwas_s_regions>0) %>% mutate(f = gwas_s_regions_smultixcan_eqtl/gwas_s_regions) %>% .$f %>% median() %>% cat("\n")

cat("median proportion of gwas loci with smultixcan enlocalized genes:\n")
d %>% filter(gwas_s_regions>0) %>% mutate(f = gwas_s_regions_smultixcan_enlocalized_eqtl/gwas_s_regions) %>% .$f %>% median() %>% cat("\n")

cat("median proportion of gwas loci with smultixcan introns:\n")
d %>% filter(gwas_s_regions>0) %>% mutate(f = gwas_s_regions_smultixcan_sqtl/gwas_s_regions) %>% .$f %>% median() %>% cat("\n")

cat("median proportion of gwas loci with smultixcan enlocalized introns:\n")
d %>% filter(gwas_s_regions>0) %>% mutate(f = gwas_s_regions_smultixcan_enlocalized_sqtl/gwas_s_regions) %>% .$f %>% median() %>% cat("\n")
