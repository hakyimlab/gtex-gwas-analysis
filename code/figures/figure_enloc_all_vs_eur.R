suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

FOLDER <-"output"
dir.create(FOLDER, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(FOLDER, p)

gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as phenotype, Deflation as deflation, new_abbreviation as abbreviation
                       FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}
                       WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

gtex_metadata <- (function(){
  query <- glue::glue("SELECT tissue, tissue_name, tissue_abbrv, tissue_color_hex, v8_eur
                       FROM {gtex_tissue_metadata_tbl$dataset_name}.{gtex_tissue_metadata_tbl$table_name}
                       WHERE v8_eur>50")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

gene_metadata <- (function(){
  query <- glue::glue("SELECT gene_id, gene_name, gene_type
                       FROM {gencode_all_annotation_tbl$dataset_name}.{gencode_all_annotation_tbl$table_name}
                       WHERE gene_type in ('pseudogene', 'protein_coding', 'lincRNA')")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

trait_whitelist <- "('GIANT_HEIGHT', 'CARDIoGRAM_C4D_CAD_ADDITIVE', 'UKB_20002_1111_self_reported_asthma')"
tissue_whitelist <- "('Kidney_Cortex', 'Brain_Cerebellum', 'Muscle_Skeletal')"

all <- (function() {
  query_ <- glue::glue(
    "SELECT gene_id, rcp, phenotype, tissue
             FROM {enloc_tbl_eqtl$dataset_name}.{enloc_tbl_eqtl$table_name} as e
             WHERE phenotype in {trait_whitelist} AND tissue in {tissue_whitelist}"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

eur <- (function() {
  query_ <- glue::glue(
    "SELECT molecular_qtl_trait as gene_id, locus_rcp as rcp, phenotype, tissue
             FROM {enloc_tbl_eqtl_eur$dataset_name}.{enloc_tbl_eqtl_eur$table_name} as e
             WHERE phenotype in {trait_whitelist} AND tissue in {tissue_whitelist}"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()


p_ <- (function(){
  d <- all %>% full_join(eur, by=c("gene_id", "phenotype", "tissue")) %>%
    rename(all=rcp.x, eur=rcp.y) %>%
    mutate(all=ifelse(is.na(all), 0, all)) %>%
    mutate(eur=ifelse(is.na(eur), 0, eur)) %>%
    mutate(complete=ifelse(all>0&eur>0, "complete", "incomplete")) %>%
    filter(gene_id %in% gene_metadata$gene_id) %>%
    inner_join(gwas_metadata %>% select(phenotype, abbreviation), by="phenotype") %>%
    select(-phenotype) %>% rename(phenotype=abbreviation) %>%
    inner_join(gtex_metadata %>% select(tissue, tissue_abbrv), by="tissue") %>%
    select(-tissue) %>% rename(tissue=tissue_abbrv) %>%
    mutate(phenotype = factor(phenotype, levels=c("CAD", "ATH_UKBS", "HEIGHT"))) %>%
    mutate(tissue = factor(tissue, levels=c("KDNCTX", "BRNCHA", "MSCLSK")))

  ggplot(d) + theme_bw(base_size=14) +
    theme(axis.text.x = element_text(angle = 45)) +
    geom_abline(intercept=0, slope=1) +
    geom_point(aes(x=all,y=eur,color=complete), show.legend = FALSE) +
    scale_color_manual(values = c("complete"="black", "incomplete"="darkgray")) +
    facet_grid(tissue ~phenotype) +
    xlab("All individuals") + ylab("Europeans") +
    ggtitle("ENLOC results agreement", subtitle="all individuals vs Europeans-only")
})()

#save_plot(p_, fp_("ENLOC_ALL_EUR.png"), 600, 600)


all_ <- (function() {
  query_ <- glue::glue(
    "SELECT gene_id, rcp, phenotype, tissue
             FROM {enloc_tbl_eqtl$dataset_name}.{enloc_tbl_eqtl$table_name} where rcp>0.5"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

eur_ <- (function() {
  query_ <- glue::glue(
    "SELECT molecular_qtl_trait as gene_id, locus_rcp as rcp, phenotype, tissue
             FROM {enloc_tbl_eqtl_eur$dataset_name}.{enloc_tbl_eqtl_eur$table_name} WHERE locus_rcp>0.5"
  )
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

f_ <- (function(){
  a_ <- all_ %>% filter(gene_id %in% gene_metadata$gene_id) %>% group_by(phenotype, tissue) %>% arrange(-rcp) %>% slice(1) %>% ungroup
  e_ <- eur_ %>% group_by(phenotype, tissue) %>% arrange(-rcp) %>% slice(1) %>% ungroup
  d_ <- a_ %>% select(gene_id, phenotype) %>% mutate(all=1) %>%
    full_join(e_ %>% select(gene_id, phenotype) %>% mutate(eur=1), by=c("phenotype", "gene_id"))
  d_[is.na(d_)] <- 0

  d_ <- d_ %>% mutate(type=ifelse(all & eur, "both\npopulations",
                     ifelse(all, "all\nindividuals", "Europeans\nonly")))
  ggplot(d_) + theme_bw(base_size=14) + geom_bar(aes(type, fill=type), show.legend = FALSE) + scale_color_viridis_d() +
    xlab("population") +
    ggtitle("ENLOC detections", subtitle = "gene-phenotype pairs with rcp>0.5 across tissues")
})()

plot_ <- plot_grid(p_, f_, labels = c("A", "B"), ncol=2)
cowplot::save_plot(fp_("ENLOC_ALL_EUR.png"), plot_, base_width =12, base_height=6)
