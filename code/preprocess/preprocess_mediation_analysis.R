suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(argparse)
  library(glue)
})

file_pattern <- "/Volumes/im-lab/nas40t2/yanyul/GTExV8/mediation_analysis_new_setting/expression/pack_into_dataframe/packed_rds_by_tissue.{tissue}.rds"
phenotype_lst <- readLines("../data/phenotype_list.txt")
tissue_lst <- readLines("../data/tissue_list.txt")
metadata <- read_tsv("../data/gwas_metadata.txt") %>% rename(phenotype=Tag)

i=1
for (tissue in tissue_lst)
{
  df_tissue <- readRDS(glue::glue(file_pattern))
  print(tissue)
  for (phenotype in phenotype_lst)
  {
    df <- df_tissue[[phenotype]]
    enloc_genes <- df[df$rcp > 0.1,] %>% select(gene_id) %>% distinct()
    df$gwas_effect_size[is.na(df$gwas_effect_size)] <- df$gwas_effect_size_imputed[is.na(df$gwas_effect_size)] #replace NA's with imputed values
    df['eqtl_effect_size'] <- df['eqtl_effect_size']/sd(df$eqtl_effect_size, na.rm=TRUE)
    df['gwas_effect_size'] <- df['gwas_effect_size']/sd(df$gwas_effect_size, na.rm=TRUE)
    df['delta_rand'] <- rnorm(nrow(df),mean(df$gwas_effect_size, na.rm=TRUE),sd(df$gwas_effect_size, na.rm=TRUE))
    df['eqtl_zscore'] <- df$eqtl_effect_size/df$slope_se

    df_2 <- df %>%
      group_by(tissue, phenotype, gene_id) %>%
      filter(n() == 2)

    if (nrow(df_2) == 0) {next}

    df_2_eqtl <- df_2 %>%
      mutate(rank = order(order(abs(eqtl_effect_size), decreasing=TRUE))) %>%
      mutate(rank = paste0("rank", rank)) %>% # rank1 = primary, rank2 = secondary
      mutate(beta_gene=gwas_effect_size/eqtl_effect_size) %>%
      mutate(beta_gene_rand=delta_rand/eqtl_effect_size)
    df_2_eqtl_1 <- df_2_eqtl %>%
      select(tissue, phenotype, gene_id, rank, beta_gene) %>%
      spread(key=rank, value=beta_gene) %>% rename("effect_size_rank1_beta"="rank1","effect_size_rank2_beta"="rank2")
    df_2_eqtl_2 <- df_2_eqtl %>%
      select(tissue, phenotype, gene_id, rank, gtex_variant_id) %>%
      spread(key=rank, value=gtex_variant_id) %>% rename("effect_size_rank1_variant_id"="rank1","effect_size_rank2_variant_id"="rank2")
    df_2_eqtl_3 <- df_2_eqtl %>%
      select(tissue, phenotype, gene_id, rank, gwas_effect_size) %>%
      spread(key=rank, value=gwas_effect_size) %>% rename("effect_size_rank1_delta"="rank1","effect_size_rank2_delta"="rank2")
    df_2_eqtl_4 <- df_2_eqtl %>%
      select(tissue, phenotype, gene_id, rank, eqtl_effect_size) %>%
      spread(key=rank, value=eqtl_effect_size) %>% rename("effect_size_rank1_gamma"="rank1","effect_size_rank2_gamma"="rank2")
    df_2_eqtl_5 <- df_2_eqtl %>%
      select(tissue, phenotype, gene_id, rank, beta_gene_rand) %>%
      spread(key=rank, value=beta_gene_rand) %>% rename("effect_size_rank1_beta_rand"="rank1","effect_size_rank2_beta_rand"="rank2")

    df_2_eqtl <- merge(x=df_2_eqtl_1,y=df_2_eqtl_2,by=c("tissue","phenotype","gene_id"))
    df_2_eqtl <- merge(x=df_2_eqtl,y=df_2_eqtl_3,by=c("tissue","phenotype","gene_id"))
    df_2_eqtl <- merge(x=df_2_eqtl,y=df_2_eqtl_4,by=c("tissue","phenotype","gene_id"))
    df_2_eqtl <- merge(x=df_2_eqtl,y=df_2_eqtl_5,by=c("tissue","phenotype","gene_id"))

    tmp <- df_2 %>%
      mutate(rank = order(order(abs(eqtl_zscore), decreasing=TRUE))) %>%
      mutate(rank = paste0("rank", rank)) %>% # rank1 = primary, rank2 = secondary
      mutate(beta_gene=gwas_effect_size/eqtl_effect_size) %>%
      mutate(beta_gene_rand=delta_rand/eqtl_effect_size)
    tmp_1 <- tmp %>%
      select(tissue, phenotype, gene_id, rank, beta_gene) %>%
      spread(key=rank, value=beta_gene) %>% rename("zscore_rank1_beta"="rank1","zscore_rank2_beta"="rank2")
    tmp_2 <- tmp %>%
      select(tissue, phenotype, gene_id, rank, gtex_variant_id) %>%
      spread(key=rank, value=gtex_variant_id) %>% rename("zscore_rank1_variant_id"="rank1","zscore_rank2_variant_id"="rank2")
    tmp_3 <- tmp %>%
      select(tissue, phenotype, gene_id, rank, gwas_effect_size) %>%
      spread(key=rank, value=gwas_effect_size) %>% rename("zscore_rank1_delta"="rank1","zscore_rank2_delta"="rank2")
    tmp_4 <- tmp %>%
      select(tissue, phenotype, gene_id, rank, eqtl_effect_size) %>%
      spread(key=rank, value=eqtl_effect_size) %>% rename("zscore_rank1_gamma"="rank1","zscore_rank2_gamma"="rank2")
    tmp_5 <- tmp %>%
      select(tissue, phenotype, gene_id, rank, beta_gene_rand) %>%
      spread(key=rank, value=beta_gene_rand) %>% rename("zscore_rank1_beta_rand"="rank1","zscore_rank2_beta_rand"="rank2")

    tmp <- merge(x=tmp_1,y=tmp_2,by=c("tissue","phenotype","gene_id"))
    tmp <- merge(x=tmp,y=tmp_3,by=c("tissue","phenotype","gene_id"))
    tmp <- merge(x=tmp,y=tmp_4,by=c("tissue","phenotype","gene_id"))
    tmp <- merge(x=tmp,y=tmp_5,by=c("tissue","phenotype","gene_id"))

    df_2_eqtl <- merge(x=df_2_eqtl,y=tmp,by=c("tissue","phenotype","gene_id"))

    df_2_eqtl <- df_2_eqtl[df_2_eqtl$effect_size_rank1_beta > quantile(df_2_eqtl$effect_size_rank1_beta, 0.05, na.rm = TRUE) & df_2_eqtl$effect_size_rank1_beta < quantile(df_2_eqtl$effect_size_rank1_beta, 0.95, na.rm = TRUE),]
    df_2_eqtl <- df_2_eqtl[df_2_eqtl$effect_size_rank2_beta > quantile(df_2_eqtl$effect_size_rank2_beta, 0.05, na.rm = TRUE) & df_2_eqtl$effect_size_rank2_beta < quantile(df_2_eqtl$effect_size_rank2_beta, 0.95, na.rm = TRUE),]

    tmp <- inner_join(enloc_genes, df_2_eqtl, by="gene_id")

    if (i == 1) {
      write.table(tmp, file = "../data/eur_pvs_analysis/eur_genes_rcp0.1.txt",row.names=FALSE,quote=FALSE,sep="\t")
    } else {
      write.table(tmp, file = "../data/eur_pvs_analysis/eur_genes_rcp0.1.txt",row.names=FALSE,quote=FALSE,sep="\t",append=TRUE,col.names=FALSE)
    }
    i = i + 1
  }
}


