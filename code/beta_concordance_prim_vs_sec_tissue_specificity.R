suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(argparse)
  library(glue)
})

enloc <- read_tsv("../data/eur_pvs_analysis/eur_genes_all.txt")

enloc['index'] <- seq.int(nrow(enloc))
tissue_sharing <- read_tsv("../data/tissue_sharing_dapg_eqtl_for_Josh/DAPG-tissue-sharing-eQTL.txt", col_names=FALSE)
tissue_specific <- read_tsv("../data/tissue_sharing_dapg_eqtl_for_Josh/DAPG-tissue-specific-eQTL.txt", col_names=FALSE)
d_ <- rbind(tissue_sharing,tissue_specific)
d_ <- d_[!names(d_) %in% "X4"]
d_ <- d_ %>% rename("variant_id"="X1","gene_id"="X2", "tissue_sharing"="X3")

effect_size_rank1 <- enloc[names(enloc) %in% c("effect_size_rank1_variant_id","gene_id","index")]
effect_size_rank2 <- enloc[names(enloc) %in% c("effect_size_rank2_variant_id","gene_id","index")]
zscore_rank1 <- enloc[names(enloc) %in% c("zscore_rank1_variant_id","gene_id","index")]
zscore_rank2 <- enloc[names(enloc) %in% c("zscore_rank2_variant_id","gene_id","index")]

effect_size_rank1 <- merge(x=effect_size_rank1, y=d_, by.x=c("gene_id","effect_size_rank1_variant_id"), by.y=c("gene_id","variant_id"), all.x=TRUE)
effect_size_rank2 <- merge(x=effect_size_rank2, y=d_, by.x=c("gene_id","effect_size_rank2_variant_id"), by.y=c("gene_id","variant_id"), all.x=TRUE)
zscore_rank1 <- merge(x=zscore_rank1, y=d_, by.x=c("gene_id","zscore_rank1_variant_id"), by.y=c("gene_id","variant_id"), all.x=TRUE)
zscore_rank2 <- merge(x=zscore_rank2, y=d_, by.x=c("gene_id","zscore_rank2_variant_id"), by.y=c("gene_id","variant_id"), all.x=TRUE)

effect_size_rank1 <- effect_size_rank1[names(effect_size_rank1) %in% c("index", "tissue_sharing")] %>% rename("effect_size_rank1_tissue_sharing" = "tissue_sharing")
effect_size_rank2 <- effect_size_rank2[names(effect_size_rank2) %in% c("index", "tissue_sharing")] %>% rename("effect_size_rank2_tissue_sharing" = "tissue_sharing")
zscore_rank1 <- zscore_rank1[names(zscore_rank1) %in% c("index", "tissue_sharing")] %>% rename("zscore_rank1_tissue_sharing" = "tissue_sharing")
zscore_rank2 <- zscore_rank2[names(zscore_rank2) %in% c("index", "tissue_sharing")] %>% rename("zscore_rank2_tissue_sharing" = "tissue_sharing")

enloc <- merge(x=enloc, y=effect_size_rank1, by="index")
enloc <- merge(x=enloc, y=effect_size_rank2, by="index")
enloc <- merge(x=enloc, y=zscore_rank1, by="index")
enloc <- merge(x=enloc, y=zscore_rank2, by="index")

enloc <- enloc[!names(enloc) %in% c("index","PP_H4_abf","rcp","nindep")]
enloc['effect_size_share_group'] <- ifelse(enloc$effect_size_rank1_tissue_sharing==TRUE & enloc$effect_size_rank2_tissue_sharing==TRUE, "Both",
                                    ifelse(enloc$effect_size_rank1_tissue_sharing==TRUE & enloc$effect_size_rank2_tissue_sharing==FALSE, "Prim only",
                                    ifelse(enloc$effect_size_rank1_tissue_sharing==FALSE & enloc$effect_size_rank2_tissue_sharing==TRUE, "Sec only",
                                    ifelse(enloc$effect_size_rank1_tissue_sharing==FALSE & enloc$effect_size_rank2_tissue_sharing==FALSE, "Neither",
                                    NA))))
enloc['zscore_share_group'] <- ifelse(enloc$zscore_rank1_tissue_sharing==TRUE & enloc$zscore_rank2_tissue_sharing==TRUE, "Both",
                               ifelse(enloc$zscore_rank1_tissue_sharing==TRUE & enloc$zscore_rank2_tissue_sharing==FALSE, "Prim only",
                               ifelse(enloc$zscore_rank1_tissue_sharing==FALSE & enloc$zscore_rank2_tissue_sharing==TRUE, "Sec only",
                               ifelse(enloc$zscore_rank1_tissue_sharing==FALSE & enloc$zscore_rank2_tissue_sharing==FALSE, "Neither",
                               NA))))
enloc <- enloc %>% drop_na(zscore_share_group)

#CONCORDANCE BASED ON DISTANCE FROM IDENTITY
#distance <- function(x,y){
#  return (abs((-x+y)/sqrt(2)))
#}
#enloc['zscore_concordance'] <- ifelse(distance(enloc$zscore_rank1_beta,enloc$zscore_rank2_beta) <= 1, "concordant",
#                                      "discordant")
#enloc['zscore_rand_concordance'] <- ifelse(distance(enloc$zscore_rank1_beta_rand,enloc$zscore_rank2_beta_rand) <= 1, "concordant",
#                                           "discordant")

#CONE SHAPED CONCORDANCE
#enloc['zscore_concordance'] <- ifelse((enloc$zscore_rank2_beta <= 2*enloc$zscore_rank1_beta & enloc$zscore_rank2_beta >= 0.5*enloc$zscore_rank1_beta) | 
#                                        (enloc$zscore_rank2_beta >= 2*enloc$zscore_rank1_beta & enloc$zscore_rank2_beta <= 0.5*enloc$zscore_rank1_beta), "concordant",
#                                      "discordant")
#enloc['zscore_rand_concordance'] <- ifelse((enloc$zscore_rank2_beta_rand <= 2*enloc$zscore_rank1_beta_rand & enloc$zscore_rank2_beta_rand >= 0.5*enloc$zscore_rank1_beta_rand) |
#                                             (enloc$zscore_rank2_beta_rand >= 2*enloc$zscore_rank1_beta_rand & enloc$zscore_rank2_beta_rand <= 0.5*enloc$zscore_rank1_beta_rand) , "concordant",
#                                      "discordant")
#enloc['effect_size_concordance'] <- ifelse((enloc$effect_size_rank2_beta <= 2*enloc$effect_size_rank1_beta & enloc$effect_size_rank2_beta >= 0.5*enloc$effect_size_rank1_beta) | 
#                                        (enloc$effect_size_rank2_beta >= 2*enloc$effect_size_rank1_beta & enloc$effect_size_rank2_beta <= 0.5*enloc$effect_size_rank1_beta), "concordant",
#                                      "discordant")
#enloc['effect_size_rand_concordance'] <- ifelse((enloc$effect_size_rank2_beta_rand <= 2*enloc$effect_size_rank1_beta_rand & enloc$effect_size_rank2_beta_rand >= 0.5*enloc$effect_size_rank1_beta_rand) |
#                                             (enloc$effect_size_rank2_beta_rand >= 2*enloc$effect_size_rank1_beta_rand & enloc$effect_size_rank2_beta_rand <= 0.5*enloc$effect_size_rank1_beta_rand) , "concordant",
#                                           "discordant")

#ON QUADRANTS
enloc['zscore_concordance'] <- ifelse((enloc$zscore_rank1_beta*enloc$zscore_rank2_beta >= 0), "concordant",
                                           "discordant")
enloc['zscore_rand_concordance'] <- ifelse((enloc$zscore_rank1_beta_rand*enloc$zscore_rank2_beta_rand >= 0), "concordant",
                                      "discordant")

write.table(enloc, file="../data/eur_pvs_analysis/eur_genes_tissue_spec_all_quadrants.txt", quote=FALSE, row.names=FALSE, sep="\t")  
  
  
  
  