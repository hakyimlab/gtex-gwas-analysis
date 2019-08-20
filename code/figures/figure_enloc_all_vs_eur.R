suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(bigrquery))

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

p_ <- ggplot(d) + theme_bw(base_size=18) +
  geom_abline(intercept=0, slope=1) +
  geom_point(aes(x=all,y=eur,color=complete), show.legend = FALSE) +
  scale_color_manual(values = c("complete"="black", "incomplete"="darkgray")) +
  facet_grid(tissue ~phenotype)

save_plot(p_, fp_("ENLOC_ALL_EUR.png"), 600, 600)
