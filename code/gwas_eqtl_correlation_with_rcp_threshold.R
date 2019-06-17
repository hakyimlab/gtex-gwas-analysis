suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--tissue", default="Whole_Blood")
parser$add_argument("--phenotype", default="Astle_et_al_2016_Eosinophil_counts")
parser$add_argument("--rcp_thresholds", default=c(0,0.01, 0.05, 0.1, 0.2, 0.5))
args <- parser$parse_args()

tissue <- args$tissue
phenotype <- args$phenotype

# metadata <- read.table("data/gwas_metadata.txt", header=TRUE, stringsAsFactors=FALSE)
metadata <- read_tsv("../data/gwas_metadata.txt")
dapg_file_pattern <- "../data/dapg_variants/with_mashr_results/DAPG_with_mashr__{tissue}.rds"
enloc_file_pattern <- "../data/enloc/{tissue}_w_{short_phenotype}_enloc_output.txt.gz"
ld_file_pattern <- "../data/dapg_variants/LD/dapg_filterd_LD_{tissue}.txt"
output_file_pattern <- "../output/eqtl_gwas_correlation/{phenotype}__{tissue}__correlation.txt"

short_phenotype <- as.character(metadata[metadata$Tag==phenotype, "new_Phenotype"])

# --------------------------------------------------------------------------------- #

dd <- readRDS(glue::glue(dapg_file_pattern))[[phenotype]]
rcp_df <- read.table(glue::glue(enloc_file_pattern), header=T)
ld_df <- read.table(glue::glue(ld_file_pattern), header=T)

# gtex_tissue_metadata <-
#   basicQuery(gtex_tissue_metadata_tbl) %>%
#   select(tissue, v8_all, tissue_abbrv) %>%
#   rename(sample_size=v8_all)


dd_gt_1 <- dd %>%
  group_by(tissue, phenotype, gene_id) %>%
  filter(n() > 1) %>%
  mutate(ranking=order(order(abs(eqtl_effect_size), decreasing=TRUE))) %>%
  ungroup() %>%
  mutate(beta_gene=gwas_effect_size/eqtl_effect_size) %>%
  select(gene_id, gtex_variant_id, ranking, beta_gene) %>%
  filter(ranking <= 2)


dd_wide <-
  inner_join(
    dd_gt_1 %>% filter(ranking==1),
    dd_gt_1 %>% filter(ranking==2),
    by="gene_id"
  ) %>%
  inner_join(
    ld_df,
    by=c("gtex_variant_id.x"="SNP1", "gtex_variant_id.y"="SNP2")
  ) %>%
  inner_join(
    rcp_df %>% mutate(gene_id=substr(gene_id, 1, 15)),
    by="gene_id"
  ) %>%
  filter(
    beta_gene.x > quantile(beta_gene.x, 0.05, na.rm = T) &
    beta_gene.y > quantile(beta_gene.y, 0.05, na.rm = T) &
    beta_gene.x < quantile(beta_gene.x, 0.95, na.rm = T) &
    beta_gene.y < quantile(beta_gene.y, 0.95, na.rm = T)
  )


# --------------------------------------------------------------------------------- #

calculate_baseline <- function(df, ld, binning=seq(0, 1, 0.1)) {
  df <- df[sample(1:nrow(df)),]
  df <- mutate(df, bin=cut(r2, binning)) %>% mutate(bin=as.character(bin))

  counts <- table(cut(ld, binning))
  ratio <- min(as.integer(min(table(df$bin)/counts)), 10)
  counts <- ratio * counts # use as many SNPs as possible, keeping the relationships
  
  cc <- vector(mode="integer", length = length(seq(0, 1, 0.1)))
  names(cc) <- names(counts)
  
  dd <- df[NULL,]
  for (i in 1:nrow(df)) {
    bin <- as.character(df[i, "bin"])
    if (cc[bin] < counts[bin]) {
      dd <- rbind(dd, df[i,])
      cc[bin] <- cc[bin] + 1
    }
  }
  dd
}

# ------------------------- Iterate through RCP thresholds ------------------------- #

rows <- list()
for (rcp_threshold in args$rcp_thresholds) {
  
  dd_rcp <- dd %>%
    inner_join(
      rcp_df %>% mutate(gene_id=substr(gene_id, 1, 15)),
      by="gene_id"
    ) %>% 
    filter(rcp > rcp_threshold)
  
  dd_wide_rcp <- dd_wide %>% 
    filter(rcp > rcp_threshold)
  
  if (nrow(dd_wide_rcp) <= 1)
    next
  
  dd_baseline <- calculate_baseline(
    dd_wide %>% filter(rcp < 0.01),
    dd_wide_rcp$r2
  )
  
  conc_baseline_pearson <- cor(
    dd_baseline$beta_gene.x,
    dd_baseline$beta_gene.y,
    use="complete.obs", method="pearson"
  )
  
  conc_baseline_spearman <- cor(
    dd_baseline$beta_gene.x,
    dd_baseline$beta_gene.y,
    use="complete.obs", method="spearman"
  )
  
  corr_eqtl_gwas_pearson <- cor(
    abs(dd_rcp$eqtl_effect_size),
    abs(dd_rcp$gwas_effect_size),
    use="complete.obs", method="pearson"
  )
  
  corr_eqtl_gwas_spearman <- cor(
    abs(dd_rcp$eqtl_effect_size),
    abs(dd_rcp$gwas_effect_size),
    use="complete.obs", method="spearman"
  )
  
  rows <- c(
    rows,
    list(
      data.frame(
        "tissue"=tissue,
        "phenotype"=phenotype,
        "rcp_threshold"=rcp_threshold,
        "corr_eqtl_gwas_pearson"=corr_eqtl_gwas_pearson,
        "corr_eqtl_gwas_spearman"=corr_eqtl_gwas_spearman,
        "concordance_pearson"=cor(dd_wide_rcp$beta_gene.x, dd_wide_rcp$beta_gene.y, use="complete.obs", method="pearson"),
        "concordance_spearman"=cor(dd_wide_rcp$beta_gene.x, dd_wide_rcp$beta_gene.y, use="complete.obs", method="spearman"),
        "conc_baseline_pearson"=conc_baseline_pearson,
        "conc_baseline_spearman"=conc_baseline_spearman,
        "n"=nrow(dd_wide_rcp),
        "ld_1st_quartile"=quantile(dd_wide_rcp$r2, 0.25, na.rm=T),
        "ld_median"=median(dd_wide_rcp$r2, na.rm=T),
        "ld_3rd_quartile"=quantile(dd_wide_rcp$r2, 0.75, na.rm=T)
      )
    )
  )
}
output_df <- bind_rows(rows)

output_filename <- glue::glue(output_file_pattern)
if (!dir.exists(dirname(output_filename))) {
  dir.create(dirname(output_filename), recursive = TRUE)
}

if (nrow(output_df) != 0) {
  write.table(
    output_df, 
    file = output_filename, 
    quote = FALSE, sep = "\t", row.names = FALSE
  )
}