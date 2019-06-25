
DATA <-"local"

WRAP_ <- function(x,y) {
  if (DATA == "local") {
    x
  } else {
    stop("not implemented")
  } 
}

metadata_ <- WRAP_((function(){
  r_tsv_("data/gwas_metadata.txt", col_types=cols_only(Tag="c", Phenotype="c", abbreviation="c", new_Phenotype="c", new_abbreviation="c")) %>%
    rename(tag=Tag, phenotype=Phenotype, new_phenotype=new_Phenotype)
})(),NULL)

###############################################################################

get_gwas <- WRAP_(function(file_logic, trait_, fields=cols_only(gtex_variant_id="c", pvalue="d", imputation_status="c")){
  file_logic %>% filter(trait == trait_) %>% .$path %>%
    r_tsv_(col_types = fields) %>% mutate(trait = trait_)
},NULL)

process_gwas <- WRAP_(function(gwas_results_logic, callback, row=NULL){
  if (!is.null(row)) {
    gwas_results_logic <- gwas_results_logic[row,]
  }
  for (pheno_ in gwas_results_logic$trait) {
    cat(pheno_, "\n")
    d <- get_gwas(gwas_results_logic, pheno_)
    callback(d)
  }
}, NULL)


###############################################################################
#model file logic, tissue
get_models_logic <- WRAP_(function(identifier="v6p_en"){
  m_ <- c("v6p_en"="data/models/gtex_v6p", 
          "v6p_en_np"="data/models/gtex_v6p_np",
          "v7_en"="data/models/gtex_v7",
          "v8_en"="data/models/gtex_v8",
          "v8_en_splicing"="data/models/gtex_v8_splicing",
          "v8_ols_dapgs"="data/models/gtex_v8_ols_dapgs",
          "v8_en_dapgw"="data/models/gtex_v8_dapgw",
          "v8_en_dapgw_hapmap"="data/models/gtex_v8_dapgw_hapmap",
          "v8_dapg_gene_snp"="data/models/gene_snp_dapg",
          "v8_mashr"="data/models/mashr",
          "v8_ctimp"="data/models/ctimp",
          "v8_cgp"="data/models/conditional",
          "v8_cgs"="data/models/conditional",
          "v8_cmp"="data/models/marginal",
          "v8_cdp"="data/models/dap_g",
          "geuvadis"="data/models/geuvadis")
  p_ <- c("v6p_en"="TW_(.*)_0.5.db",
          "v6p_en_np" = "gtex_v6p_(.*)_tw_0.5.db",
          "v7_en"="gtex_v7_(.*)_imputed_europeans_tw_0.5_signif.db",
          "v8_en"="gtex_v8_(.*)_itm_signif.db",
          "v8_en_splicing"="gtex_splicing_v8_eur_(.*)_signif.db",
          "v8_ols_dapgs"="dapgw_(.*).db",
          "v8_en_dapgw"="dapgw_(.*).db",
          "v8_en_dapgw_hapmap"="dapgw_(.*).db",
          "v8_dapg_gene_snp"="gene_snps_(.*).db",
          "v8_mashr"="mashr_(.*).db",
          "v8_ctimp"="ctimp_(.*).db",
          "v8_cgp"="(.*)_primary_eQTLs.db",
          "v8_cgs"="(.*)_secondary_eQTLs.db",
          "v8_cmp"="(.*)_Conditional_analysis_primary_eQTL.db",
          "v8_cdp"="(.*)_DAPG_independent_eQTL_primary_by_beta.db",
          "geuvadis"="(.*)_HapMap_alpha0.5_window1e6_filtered.db")
  m_results_logic(m_[identifier][[1]], p_[identifier][[1]]) %>% mutate(family=identifier)
}, NULL)

#model file logic, selected tissue, fields
get_model_extra <- WRAP_(function(file_logic, tissue_=NULL, fields="gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`") {
  if (!is.null(tissue_)) {
    file_logic <- file_logic %>% filter(tissue == tissue_)
  }
  l <- list()
  for (i in 1:nrow(file_logic)) {
    path <- file_logic$path[i]
    l[[i]] <- db_(path, function(con){
      query <- paste0("SELECT ", fields, " FROM extra")
      dbGetQuery(con, query)
    }) %>% mutate(tissue = file_logic$tissue[i], family = file_logic$family[1])
  }
  do.call(rbind, l)
} ,NULL)

#model file logic,  fields
get_model_extra_2 <- WRAP_(function(file_logics, fields="gene, genename, `n.snps.in.model`, `pred.perf.R2`, `pred.perf.pval`"){
  l <- list()
  i <- 1
  for (file_logic in file_logics) {
    l[[i]] <- get_model_extra(file_logic, fields=fields)
    i <- i+1
  }
  do.call(rbind, l)
} ,NULL)

get_model_weight <- WRAP_( function(file_logic, tissue_, fields="gene, rsid, varID, ref_allele, eff_allele, weight"){
  path <- file_logic %>% filter(tissue == tissue_) %>% .$path
  db_(path, function(con){
    query <- paste0("SELECT ", fields, " FROM weights")
    dbGetQuery(con, query)
  })
},NULL)

#model file logic,  fields
get_model_weight_2 <- WRAP_(function(file_logics, fields="gene, rsid, varID, ref_allele, eff_allele, weight"){
  d <- list()
  i <- 1
  for (file_logic in file_logics) {
    family <- file_logic$family %>% unique
    cat("Family: ", family, "\n")
    for (tissue in file_logic$tissue)  {
      d[[i]] <- get_model_weight(file_logic, tissue) %>% 
        mutate(family = family, tissue=tissue)
      i <- i+1
    }
  }
  do.call(rbind, d)
} ,NULL)

###############################################################################

#no args
get_coloc_results_logic <- WRAP_(function(identifier="gtex_v8") {
      if (identifier == "gtex_v8") {
        cpt_results_logic("data/coloc_v8", "(.*)_w_(.*)_coloc_output.txt.gz", "_coloc_output.txt.gz") %>% 
          inner_join(metadata_ %>% select(tag, new_phenotype), by=c("pheno"="new_phenotype")) %>% mutate(pheno=tag) %>% select(-tag)
      } else if(identifier == "ep_v8") {
        cpt_results_logic("data/coloc_v8_ep", "(.*)__PM__(.*).txt.gz", ".txt.gz") %>% rename(pheno=tissue, tissue=pheno)
      } else if(identifier == "default_v8") {
        cpt_results_logic("data/coloc_v8_default", "(.*)__PM__(.*).txt.gz", ".txt.gz") %>% rename(pheno=tissue, tissue=pheno)
      } else if(identifier == "ep_v8_gp_eqtlp") {
        cpt_results_logic("data/coloc_v8_ep_gp_etlp", "(.*)__PM__(.*).txt.gz", ".txt.gz") %>% rename(pheno=tissue, tissue=pheno)
      } else  {
        #implement
        stop("unimplemented")
      }
    }, NULL)

read_gtex_coloc_ <- function(path) {
  path %>% r_tsv_(col_types=cols_only(gene_id="c",  PP.H0.abf="d", PP.H1.abf="d", PP.H2.abf="d", PP.H3.abf="d", PP.H4.abf="d")) %>%
    rename(p0=PP.H0.abf, p1=PP.H1.abf, p2=PP.H2.abf, p3=PP.H3.abf, p4=PP.H4.abf, gene=gene_id) %>% 
    mutate(p0=as.double(p0), p1=as.double(p1), p2=as.double(p2), p3=as.double(p3),p4=as.double(p4)) %>%
    mutate(p012 = p0+p1+p2, pc=p4)
}

read_coloc_ <- function(path) {
  path %>% r_tsv_() %>% mutate(p012 = p0+p1+p2, pc=p4)
}

# result_logic, pheno, tissue
get_coloc_result_ <- WRAP_(function(file_logic, pheno_=NULL, tissue_=NULL, mode="gtex_run") {
  path <- file_logic %>% filter(pheno == pheno_, tissue == tissue_) %>% .$path 
  if (length(path) == 0) 
    return(data.frame(gene=character(0), pheno=character(0), tissue=character(0), p0=double(0), p1=double(0), p2=double(0), p3=double(0), p4=double(0), p012=double(0), pc=double(0)))
  
  if (mode == "gtex_run") {
    path %>% read_gtex_coloc_ %>% mutate(pheno=pheno_, tissue=tissue_)
  } else {
    path %>% read_coloc_ %>% mutate(pheno=pheno_, tissue=tissue_)
  }
}, NULL)

get_coloc_result <- WRAP_(function(file_logic, pheno_=NULL, tissue_=NULL, type="gtex") {
  if (!is.null(pheno_)) {
    file_logic <- file_logic %>% filter(pheno %in% pheno_)
  }
  if (!is.null(tissue_)) {
    file_logic <- file_logic %>% filter(tissue %in% tissue_)
  }
  r <- list()
  i <- 1
  for (trait_ in unique(file_logic$pheno)) {
    for (tiss_ in unique(file_logic$tissue)) {
      r[[i]] <- get_coloc_result_(file_logic, trait_, tiss_, type) %>% mutate(trait = trait_, tissue = tiss_)
      i <- i+1
    }
  }
  do.call(rbind, r)
}, NULL)

###############################################################################

get_enloc_results_logic <- WRAP_(function(identifier="v8_expression"){
  if (identifier == "v8_expression") {
    cpt_results_logic("data/enloc_v8", "(.*)_w_(.*)_enloc_output.txt.gz", "_enloc_output.txt.gz") %>% 
      inner_join(metadata_ %>% select(tag, new_phenotype), by=c("pheno"="new_phenotype")) %>% mutate(trait=tag) %>% select(-tag)
  } else if (identifier == "v8_splicing"){
    d_results_logic_("data/enloc_sqtl", "(.*)__PM__(.*).enloc.rst.gz", ".enloc.rst.gz") %>% rename(trait = name_1, tissue = name_2)
  } else {
    stop("unsupported identifier")
  }
}, NULL)

read_gtex_enloc_ <- function(path) {
  #One gene might be tested against multiple variants. Keep only the best.
  path %>% r_tsv_(col_types=cols_only(gene_id="c",  rcp="d", experiment_rcp="d")) %>% rename( gene=gene_id) %>%
    group_by(gene) %>% arrange(-rcp) %>% slice(1) %>% ungroup
}

read_vanilla_enloc_ <- function(path) {
  #One gene might be tested against multiple variants. Keep only the best.
  path %>% r_tsv_%>%  select(gene=molecular_qtl_trait, rcp=locus_rcp) %>% mutate(experiment_rcp = NA) %>%
    group_by(gene) %>% arrange(-rcp) %>% slice(1) %>% ungroup
}


get_enloc_result_ <- WRAP_(function(path, type="gtex") {
  if (type == "gtex")
    path %>% read_gtex_enloc_
  else if (type =="vanilla")
    path %>% read_vanilla_enloc_
  else 
    stop("unsupported type")
}, NULL)

get_enloc_result <- WRAP_(function(file_logic, phenos=NULL, tissues=NULL, type="gtex") {
  if (!is.null(phenos)) {
    file_logic <- file_logic %>% filter(trait %in% phenos)
  }
  if (!is.null(tissues)) {
    file_logic <- file_logic %>% filter(tissue %in% tissues)
  }
  r <- list()
  i <- 1
  for (trait_ in unique(file_logic$trait)) {
    for (tissue_ in unique(file_logic$tissue)) {
      r[[i]] <- file_logic %>% filter(trait == trait_, tissue == tissue_) %>% .$path %>% get_enloc_result_(type) %>% mutate(trait = trait_, tissue = tissue_)
      i <- i+1
    }
  }
  do.call(rbind, r)
}, NULL)


###############################################################################

#identifier
get_predixcan_results_logic <- WRAP_(function(identifier="v8_en"){
    m_ <- c("v6p_en"="data/spredixcan/sp_imputed_gwas_gtexv6p_en_old", 
            "v7_en"="data/spredixcan/sp_imputed_gwas_gtexv7_en", 
            "v8_en"="data/spredixcan/sp_imputed_gwas_gtexv8_en", 
            "v8_en_splicing"="data/spredixcan/sp_imputed_gwas_gtexv8_en_splicing",
            "v8_ols_dapgs"="data/spredixcan/sp_imputed_gwas_gtexv8_ols_dapgs",
            "v8_en_dapgw"="data/spredixcan/sp_imputed_gwas_gtexv8_en_dapgw",
            "v8_en_dapgw_hapmap"="data/spredixcan/sp_imputed_gwas_gtexv8_en_dapgw_hapmap",
            "v8_dapg_gene_snp"="data/spredixcan/sp_imputed_gwas_gtex_dapg_gene_snp",
            "v8_mashr"="data/spredixcan/sp_imputed_gwas_gtexv8_mashr",
            "v8_ctimp"="data/spredixcan/sp_imputed_gwas_gtexv8_ctimp",
            "v8_cgp"="data/spredixcan/sp_imputed_gwas_conditional_gtex_primary",
            "v8_cmp"="data/spredixcan/sp_imputed_gwas_conditional_marginal_gtex_primary",
            "v8_cdp"="data/spredixcan/sp_imputed_gwas_dapg_gtex_primary",
            "geuvadis"="data/spredixcan/sp_imputed_gwas_geuvintron_en_l",
            "streamed"="data/spredixcan/sp_stream")
    p_ <- c("v6p_en"="spredixcan_igwas_gtexenv6pold_(.*)__PM__(.*).csv",
            "v7_en"="spredixcan_igwas_gtexenv7_(.*)__PM__(.*).csv",
            "v8_en"="spredixcan_igwas_gtexenv8_(.*)__PM__(.*).csv",
            "v8_en_splicing"="spredixcan_igwas_gtexenv8_splicing_(.*)__PM__(.*).csv",
            "v8_ols_dapgs"="spredixcan_igwas_gtexolsdapgsv8_(.*)__PM__(.*).csv",
            "v8_en_dapgw"="spredixcan_igwas_gtexendapgwv8_(.*)__PM__(.*).csv",
            "v8_en_dapgw_hapmap"="spredixcan_igwas_gtexenv8dapgwhm_(.*)__PM__(.*).csv",
            "v8_dapg_gene_snp"="spredixcan_igwas_dgs_(.*)__PM__(.*).csv",
            "v8_mashr"="spredixcan_igwas_gtexmashrv8_(.*)__PM__(.*).csv",
            "v8_ctimp"="spredixcan_igwas_gtexctimpv8_(.*)__PM__(.*).csv",
            "v8_cgp"="spredixcan_igwas_cgtexp_(.*)__PM__(.*).csv",
            "v8_cmp"="spredixcan_igwas_cmgtexp_(.*)__PM__(.*).csv",
            "v8_cdp"="spredixcan_igwas_dgtexp_(.*)__PM__(.*).csv",
            "geuvadis"="spredixcan_igwas_geuvin_(.*)__PM__(.*).csv",
            "streamed"="spredixcan_igwas_gtexenv8_(.*)__PM__(.*).csv")
    p_results_logic(m_[identifier][[1]], p_[identifier][[1]], ".csv") %>% mutate(family=identifier)
  }, NULL)

#result_logic, phenoss, tissues
get_predixcan_result <- WRAP_(function(file_logic, phenos=NULL, tissues=NULL, col_types=cols_only(gene="c", zscore="d", pvalue="d", pred_perf_pval="d")) {
  if (!is.null(phenos)) {
    file_logic <- file_logic %>% filter(pheno %in% phenos)
  }
  if (!is.null(tissues)) {
    file_logic <- file_logic %>% filter(tissue %in% tissues)
  }
  d <- data.frame()
  
  r <- list()
  i <- 1
  for (pheno_ in file_logic$pheno %>% unique ){
    for (tissue_ in file_logic$tissue %>% unique) {
      #cat(pheno_, ":", tissue_, "\n")
      f_ <- file_logic %>% filter(pheno == pheno_, tissue == tissue_)
      d_ <- f_$path %>%  r_csv_(col_types=col_types) %>%
        mutate(pheno=pheno_, tissue=tissue_, family=f_$family)
      d_ <- if ("zscore" %in% colnames(d_)) {
        d_ %>% filter(!is.na(zscore))    
      } else if ("pvalue" %in% colnames(d_)) {
        d_ %>% filter(!is.na(pvalue))    
      } else {
        d_
      }

      r[[i]] <- d_
      i <- i+1
    }
  }
  do.call(rbind, r)
}, NULL)

get_predixcan_result_pheno <- function(file_logic, pheno_, col_types=cols_only(gene="c", zscore="d", pvalue="d", pred_perf_pval="d")) {
  f <- file_logic %>% filter(pheno == pheno_)
  get_predixcan_result(f, tissues=f$tissue %>% unique, col_types=col_types)
}

get_predixcan_result_pheno_2 <- function(file_logics, pheno_, col_types=cols_only(gene="c", zscore="d", pvalue="d", pred_perf_pval="d")) {
  d <- data.frame()
  r <- list()
  for (l_ in file_logics) {
    d_ <- get_predixcan_result_pheno(l_, pheno_, col_types=col_types) %>% mutate(family = unique(l_$family))
    d <- rbind(d,d_)
  }
  d
}

# result_logics (list of logic), pheno, tissue
get_predixcan_result_2 <- WRAP_(function(file_logics, pheno_=NULL, tissue_=NULL, col_types=cols_only(gene="c", zscore="d", pvalue="d", pred_perf_pval="d")) {
  d <- data.frame()
  for (l_ in file_logics) {
    d_ <- get_predixcan_result(l_, pheno_, tissue_, col_types=col_types) %>% mutate(family = unique(l_$family))
    d <- rbind(d,d_)
  }
  d
}, NULL)


#result_logics (list of logic), tissue, callback
process_predixcan_trait_2 <- WRAP_(function(predixcan_results_logics, callback, selected_tissue=NULL){
  phenos_ <- predixcan_results_logics[[1]]$pheno %>% unique()
  for (pheno_ in phenos_) {
    cat(pheno_, "\n")
    if (!is.null(selected_tissue)) {
      d <- get_predixcan_result_2(predixcan_results_logics, pheno_, selected_tissue)  
    } else {
      d <- get_predixcan_result_pheno_2(predixcan_results_logics, pheno_)
    }
    
    callback(d)
  }
}, NULL)


###############################################################################

get_predixcan_coloc_result <- WRAP_(function(predixcan_results_logic, coloc_results_logic, pheno, tissue){
  p_ <- get_predixcan_result(predixcan_results_logic, pheno, tissue)
  c_ <- get_coloc_result(coloc_results_logic, pheno, tissue)
  p_ %>% full_join(c_, by=c("gene", "pheno", "tissue"))
},NULL)

get_predixcan_coloc_result_pheno <- WRAP_(function(predixcan_results_logic, coloc_results_logic, pheno_, tissue_filter = NULL){
  tissues_ <- predixcan_results_logic$tissue %>% unique
  if (!is.null(tissue_filter)) {
    tissues_ <- tissues_[tissues_ %in% tissue_filter]
  }
  d <- data.frame()
  for (tissue_ in tissues_) {
    cat("Processing ", tissue_, "\n")
    p_ <- get_predixcan_result(predixcan_results_logic, pheno_, tissue_)
    c_ <- get_coloc_result(coloc_results_logic, pheno_, tissue_)
    d_ <- p_ %>% full_join(c_, by=c("gene", "pheno", "tissue"))
    d <- rbind(d, d_)
  }
  d
},NULL)

process_predixcan_coloc_pheno <- WRAP_(function(predixcan_results_logic, coloc_results_logic, callback, tissue_filter=NULL, pheno_filter=NULL){
  traits_ <- predixcan_results_logic$pheno %>% unique
  if (!is.null(pheno_filter)) {
    traits_ <- traits_[traits_ %in% pheno_filter]
  }
  for (trait_ in traits_) {
    cat("Processing ",trait_, "\n")
    d <- get_predixcan_coloc_result_pheno(predixcan_results_logic, coloc_results_logic, trait_, tissue_filter)
    callback(d)
  }
},NULL)

process_predixcan_coloc_pheno_tissue <- WRAP_(function(predixcan_results_logic, coloc_results_logic, callback, pheno_filter=NULL, tissue_filter=NULL){
  traits_ <- predixcan_results_logic$pheno %>% unique
  if (!is.null(pheno_filter)) {
    traits_ <- traits_[traits_ %in% pheno_filter]
  }
  for (trait_ in traits_) {
    cat("Processing ",trait_, "\n")
    pl <- predixcan_results_logic %>% filter(pheno == trait_)
    tissues_ <- pl$tissue %>% unique
    if (!is.null(tissue_filter)) {
      tissues_ <- tissues_[tissues_ %in% tissue_filter]
    }
    for (tissue_ in tissues_) {
      cat("Processing ",tissue_, "\n")
      d <- get_predixcan_coloc_result(predixcan_results_logic, coloc_results_logic, trait_, tissue_) 
      callback(d)
    }
    
  }
},NULL)

###############################################################################

get_predixcan_enloc_result <- WRAP_(function(predixcan_results_logic, enloc_results_logic, pheno, tissue){
  p_ <- get_predixcan_result(predixcan_results_logic, pheno, tissue)
  c_ <- get_enloc_result(enloc_results_logic, pheno, tissue)
  p_ %>% full_join(c_, by=c("gene", "pheno", "tissue"))
},NULL)

get_predixcan_enloc_result_pheno <- WRAP_(function(predixcan_results_logic, enloc_results_logic, pheno_, tissue_filter = NULL){
  tissues_ <- predixcan_results_logic$tissue %>% unique
  if (!is.null(tissue_filter)) {
    tissues_ <- tissues_[tissues_ %in% tissue_filter]
  }
  d <- data.frame()
  for (tissue_ in tissues_) {
    cat("Processing ", tissue_, "\n")
    p_ <- get_predixcan_result(predixcan_results_logic, pheno_, tissue_)
    c_ <- get_enloc_result(enloc_results_logic, pheno_, tissue_)
    d_ <- p_ %>% full_join(c_, by=c("gene", "pheno", "tissue"))
    d <- rbind(d, d_)
  }
  d
},NULL)

process_predixcan_enloc_pheno <- WRAP_(function(predixcan_results_logic, enloc_results_logic, callback, tissue_filter=NULL, pheno_filter=NULL){
  traits_ <- predixcan_results_logic$pheno %>% unique
  if (!is.null(pheno_filter)) {
    traits_ <- traits_[traits_ %in% pheno_filter]
  }
  for (trait_ in traits_) {
    cat("Processing ",trait_, "\n")
    d <- get_predixcan_enloc_result_pheno(predixcan_results_logic, enloc_results_logic, trait_, tissue_filter)
    callback(d)
  }
},NULL)

process_predixcan_enloc_pheno_tissue <- WRAP_(function(predixcan_results_logic, enloc_results_logic, callback, pheno_filter=NULL, tissue_filter=NULL){
  traits_ <- predixcan_results_logic$pheno %>% unique
  if (!is.null(pheno_filter)) {
    traits_ <- traits_[traits_ %in% pheno_filter]
  }
  for (trait_ in traits_) {
    cat("Processing ",trait_, "\n")
    pl <- predixcan_results_logic %>% filter(pheno == trait_)
    tissues_ <- pl$tissue %>% unique
    if (!is.null(tissue_filter)) {
      tissues_ <- tissues_[tissues_ %in% tissue_filter]
    }
    for (tissue_ in tissues_) {
      cat("Processing ",tissue_, "\n")
      d <- get_predixcan_enloc_result(predixcan_results_logic, enloc_results_logic, trait_, tissue_) 
      callback(d)
    }
    
  }
},NULL)


###############################################################################

get_coloc_enloc_result <- WRAP_(function(coloc_results_logic, enloc_results_logic, pheno, tissue){
  p_ <- get_coloc_result(coloc_results_logic, pheno, tissue)
  c_ <- get_enloc_result(enloc_results_logic, pheno, tissue)
  p_ %>% full_join(c_, by=c("gene", "pheno", "tissue"))
},NULL)

get_coloc_enloc_result_pheno <- WRAP_(function(coloc_results_logic, enloc_results_logic, pheno_, tissue_filter = NULL){
  tissues_ <- coloc_results_logic$tissue %>% unique
  if (!is.null(tissue_filter)) {
    tissues_ <- tissues_[tissues_ %in% tissue_filter]
  }
  d <- data.frame()
  for (tissue_ in tissues_) {
    cat("Processing ", tissue_, "\n")
    p_ <- get_coloc_result(coloc_results_logic, pheno_, tissue_)
    c_ <- get_enloc_result(enloc_results_logic, pheno_, tissue_)
    d_ <- p_ %>% full_join(c_, by=c("gene", "pheno", "tissue"))
    d <- rbind(d, d_)
  }
  d
},NULL)

process_coloc_enloc_pheno <- WRAP_(function(coloc_results_logic, enloc_results_logic, callback, tissue_filter=NULL, pheno_filter=NULL){
  traits_ <- coloc_results_logic$pheno %>% unique
  if (!is.null(pheno_filter)) {
    traits_ <- traits_[traits_ %in% pheno_filter]
  }
  for (trait_ in traits_) {
    cat("Processing ",trait_, "\n")
    d <- get_coloc_enloc_result_pheno(coloc_results_logic, enloc_results_logic, trait_, tissue_filter)
    callback(d)
  }
},NULL)

process_coloc_enloc_pheno_tissue <- WRAP_(function(coloc_results_logic, enloc_results_logic, callback, pheno_filter=NULL, tissue_filter=NULL){
  traits_ <- coloc_results_logic$pheno %>% unique
  if (!is.null(pheno_filter)) {
    traits_ <- traits_[traits_ %in% pheno_filter]
  }
  for (trait_ in traits_) {
    cat("Processing ",trait_, "\n")
    pl <- coloc_results_logic %>% filter(pheno == trait_)
    tissues_ <- pl$tissue %>% unique
    if (!is.null(tissue_filter)) {
      tissues_ <- tissues_[tissues_ %in% tissue_filter]
    }
    for (tissue_ in tissues_) {
      cat("Processing ",tissue_, "\n")
      d <- get_coloc_enloc_result(coloc_results_logic, enloc_results_logic, trait_, tissue_) 
      callback(d)
    }
    
  }
},NULL)

###############################################################################

process_predixcan_coloc_enloc_pheno_tissue <- WRAP_(function(predixcan_results_logic, coloc_results_logic, enloc_results_logic, callback, pheno_filter=NULL, tissue_filter=NULL){
  traits_ <- predixcan_results_logic$pheno %>% unique
  if (!is.null(pheno_filter)) {
    traits_ <- traits_[traits_ %in% pheno_filter]
  }
  for (trait_ in traits_) {
    cat("Processing ",trait_, "\n")
    pl <- predixcan_results_logic %>% filter(pheno == trait_)
    tissues_ <- pl$tissue %>% unique
    if (!is.null(tissue_filter)) {
      tissues_ <- tissues_[tissues_ %in% tissue_filter]
    }
    for (tissue_ in tissues_) {
      cat("Processing ",tissue_, "\n")
      d <- get_predixcan_coloc_enloc_result(predixcan_results_logic, coloc_results_logic, enloc_results_logic, trait_, tissue_) 
      callback(d)
    }
    
  }
},NULL)

process_predixcan_coloc_enloc_pheno <- WRAP_(function(predixcan_results_logic, coloc_results_logic, enloc_results_logic, callback, tissue_filter=NULL, pheno_filter=NULL){
  traits_ <- predixcan_results_logic$pheno %>% unique
  if (!is.null(pheno_filter)) {
    traits_ <- traits_[traits_ %in% pheno_filter]
  }
  for (trait_ in traits_) {
    cat("Processing ",trait_, "\n")
    d <- get_predixcan_coloc_enloc_result_pheno(predixcan_results_logic, coloc_results_logic, enloc_results_logic, trait_, tissue_filter)
    callback(d)
  }
},NULL)

get_predixcan_coloc_enloc_result_pheno <- WRAP_(function(predixcan_results_logic, coloc_results_logic, enloc_results_logic, pheno_, tissue_filter = NULL){
  tissues_ <- predixcan_results_logic$tissue %>% unique
  if (!is.null(tissue_filter)) {
    tissues_ <- tissues_[tissues_ %in% tissue_filter]
  }
  d <- data.frame()
  for (tissue_ in tissues_) {
    cat("Processing ", tissue_, "\n")
    p_ <- get_predixcan_result(predixcan_results_logic, pheno_, tissue_)
    c_ <- get_coloc_result(coloc_results_logic, pheno_, tissue_)
    e_ <- get_enloc_result(enloc_results_logic, pheno_, tissue_)
    d_ <- p_ %>% full_join(c_, by=c("gene", "pheno", "tissue")) %>% full_join(e_, by=c("gene", "pheno", "tissue"))
    d <- rbind(d, d_)
  }
  d
},NULL)

get_predixcan_coloc_enloc_result <- WRAP_(function(predixcan_results_logic, coloc_results_logic, enloc_results_logic, pheno, tissue){
  pc_ <- get_predixcan_coloc_result(predixcan_results_logic, coloc_results_logic, pheno, tissue)
  e_ <- get_enloc_result(enloc_results_logic, pheno, tissue)
  pc_ %>% full_join(e_, by=c("gene", "pheno", "tissue"))
},NULL)

get_intron_name_map <- function(folder) {
  f_ <- file_logic_(folder, "(.*)_key_list.txt.gz")
  d <- list()
  for (i in 1:nrow(f_))  {
      #cat(tissue, "\n")
      d_ <-f_$path[i] %>% r_tsv_() %>% mutate(tissue=f_$name[i])
      d[[i]] <- d_
  }
  do.call(rbind, d)
}
