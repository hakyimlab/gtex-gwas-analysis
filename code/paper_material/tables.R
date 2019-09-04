suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(xtable))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

RESULT<-"output/tables"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)


gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as phenotype, new_abbreviation as abbreviation, Category as category, Sample_Size as sample_size
                       FROM {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}
                       WHERE Deflation=0")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

order_ <- gwas_metadata %>% filter(!grepl("UKB", phenotype)) %>% group_by(category) %>% arrange(-sample_size) %>% slice(1) %>% ungroup %>% arrange(-sample_size) %>% .$category
all_ <- unique(gwas_metadata$category)
order_ <- order_ %>% c(all_[!(all_ %in% order_)])

t <- gwas_metadata %>% mutate(category=factor(category, levels=order_)) %>% arrange(category, phenotype) %>%
  mutate(category=gsub("_", "-", category), phenotype=gsub("_", " ", phenotype), abbreviation=gsub("_", "\\_", abbreviation), sample_size=as.integer(sample_size))
t %>% select(Category=category, Phenotype=phenotype, Abbreviation=abbreviation, `Sample Size`=sample_size) %>% xtable %>% print.xtable(file = fp_("gwas_table.tex"), include.rownames=FALSE)
