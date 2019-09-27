###############################################################################
# USE MINICONDA
###############################################################################
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(bigrquery))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(viridis))

suppressWarnings(source("code/helpers/_helpers.R"))
suppressWarnings(source("code/helpers/_helpers_big_query.R"))
suppressWarnings(source("code/helpers/_helpers_big_query_tables.R"))

FOLDER <-"output"
dir.create(FOLDER, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(FOLDER, p)


scale_ <- viridis(4)
###############################################################################

gwas_metadata <- (function(){
  query <- glue::glue("SELECT Tag as phenotype, Deflation as deflation, new_abbreviation as abbreviation, Category as category from
                      {gwas_metadata_tbl$dataset_name}.{gwas_metadata_tbl$table_name}")
  query_exec(query, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()
pheno_whitelist <- sprintf("(%s)", toString(sprintf("'%s'", gwas_metadata$phenotype)))

gwas_count <- (function(){
  query_ <- glue::glue(
    "select phenotype, count(*) as variants
    FROM (
      SELECT * FROM
      {gwas_formatted_tbl$dataset_name}.{gwas_formatted_tbl$table_name}
      where pvalue is not NULL
    ) group by phenotype

") %>% glue::glue()
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf)
})()

gwas_decor_ <- gwas_count %>%
  inner_join(
    gwas_metadata %>% select(phenotype, abbreviation, deflation) %>% mutate(deflation=ifelse(deflation==0, "undeflated", "deflated")),
    by="phenotype") %>%
  arrange(variants)
gwas_decor <- gwas_decor_ %>% filter(deflation == "deflated") %>%
  rbind(gwas_decor_ %>% filter(phenotype == "Astle_et_al_2016_Red_blood_cell_count"))
order_ <- gwas_decor$abbreviation

###############################################################################
# I forgive you, BigQuery. I hope you can forget all the bad shit I said about you.
get_gwas_q_ <- function(tbl, pheno, the_type) {
  query_ <- glue::glue(
"with quantiles as (
  SELECT phenotype, APPROX_QUANTILES(pvalue, 4) AS values
  FROM {tbl$dataset_name}.{tbl$table_name}
  GROUP BY phenotype
) select phenotype, values, q
FROM quantiles
CROSS JOIN UNNEST(quantiles.values) as values
WITH OFFSET as q
") %>% glue::glue()
  query_exec(query_, project = "gtex-awg-im", use_legacy_sql = FALSE, max_pages = Inf) %>%
    mutate(type = the_type) %>%
    mutate(q=paste0("q",q))
}
d <- rbind(get_gwas_q_(gwas_formatted_tbl, pheno, "formatted"),
           get_gwas_q_(gwas_tbl, pheno, "imputed")) %>%
    mutate(type=ifelse(type == "formatted", "intersection", "intersection+\nimputed"))

d_ <- d %>% mutate(values=-log10(values)) %>%
  mutate(values = pmin(30, values)) %>%
  spread(key="q", value="values") %>%
  rename(y0 = q4, y25 = q3, y50 = q2, y75=q1, y100=q0) %>%
  mutate(ymin = y25 - 1.5*(y75-y25)) %>% mutate(ymin = ifelse(ymin <y0, y0, ymin)) %>%
  mutate(ymax = y75 + 1.5*(y75-y25)) %>% mutate(ymax = ifelse(ymax >y100, y100, ymax)) %>%
  inner_join(
    gwas_decor %>% select(phenotype, deflation, abbreviation),
    by="phenotype") %>%
  mutate(abbreviation=factor(abbreviation, levels=order_))

p_ <- ggplot(d_) + theme_bw() +
  theme(legend.position="bottom") + #, axis.text.x = element_text(angle = 90)) +
  coord_flip() +
  geom_boxplot(aes(abbreviation, ymin=ymin, lower=y25, middle=y50, upper=y75, ymax=ymax, fill=deflation, color=type), stat="identity") +
  scale_color_manual(values=c("intersection"=scale_[3], "intersection+\nimputed"=scale_[1])) +
  scale_fill_manual(values=c("undeflated"="white", "deflated"="lightgray"))
save_plot(p_, fp_("GWAS_DEFLATION.png"), 1200, 500)

# phenos <- c("ADIPOGen_Adiponectin")
# for (pheno in phenos) {
#   d <- get_gwas(gwas_formatted_tbl, pheno, "imputed")
# }
