# GTEx v8 GWAS analysis

Resources and analysis shared by ["Widespread dose-dependent effects of RNA expression and splicing on complex diseases and traits"](https://doi.org/10.1101/814350).

# Data

The results underlying these analyses can be found in **zenodo.org**:

* [Finemapping on eQTL and sQTL](https://zenodo.org/record/3517189#.XbMe6NF7m90): computed using [DAP-G](https://github.com/xqwen/dap)
* [GWAS and sQTL/sQTL integration](https://zenodo.org/record/3518299#.XbMgFNF7m90). This zenodo package contains the following:
    * [coloc](https://github.com/chr1swallace/coloc) results on eQTL, using priors computed from [enloc](https://github.com/xqwen/integrative) enrichment
    * [enloc](https://github.com/xqwen/integrative) results on eQTL/sQTL, using only individuals of European Ancestry and variants with MAF>0.01.
GWAS regions from [ldetect](https://bitbucket.org/nygcresearch/ldetect/src/master/), lifted over to hg38
    * expression and splicing prediction models using [MASHR](https://github.com/stephenslab/mashr) effect sizes
    * predicted expression and splicing associations (S-MultiXcan and S-Predixcan, [here](https://github.com/hakyimlab/MetaXcan)). 
    Model training, GWAS harmonization and imputation available [here](https://github.com/hakyimlab/summary-gwas-imputation)
    * [SMR](https://cnsgenomics.com/software/smr/#Overview) analysis on eQTL
* Addendum: SMR sQTL [results](https://zenodo.org/record/3525070#.XbxH79F7m90)
* [Elastic Net Prediction models](https://zenodo.org/record/3519321#.XbxAANF7m90): This provides robust, if less powerful, models than the new MASHR-based models.
Provided for compatibility.
* Harmonized/imputed GWAS: The underlying GWAS summary statistics harmonized and imputed to GTEx v8 (*you need to provide a google cloud account with billing information*) is available  [here](https://storage.googleapis.com/gtex-gwas-public/gwas/harmonized_imputed_gwas.tar).

# Reproducible analysis for manuscript

The code for the manuscript's analyses is available here as R scripts with the following dependencies:

- bigrquery
- tidyverse
- upsetR

To use BigQuery, it is helpful to install [Google Cloud SDK](https://cloud.google.com/sdk/) and have read access to the required Google Cloud tables.

# GTEx GWAS subgroup paper

The markdows can be built from an R session executing:

```R
wflow_build("analysis/miscellaneous_statistics_2.Rmd")
```

The R scripts can be run from a bash session:

```bash
Rscript code/figures/figure_enloc_all_vs_eur.R
```

## Manuscript material

Support scripts:

* `code/helpers` folder contains miscellaneous R functions and definitions used throughtout the analyses.
    * `code/helpers/_helpers_big_query_tables.R` contains a centralized definition of Bigquery  tables to be used by other scripts.
* `code/preprocess` folder contains scripts that were ran once to setup auxiliary data.
    * `code/preprocess/preprocess_gwas_regions.R` counts detections per loci (independent LD regions) that will be used as inputs for other analyses.
    * `code/preprocess/preprocess_setup_auxiliary_tables.R` builds auxiliary tables in big query
    * `code/preprocess/preprocess_mediation_analysis.R` builds data for primary vs secondary concordance analysis. Can download the data from `download_aux_data.sh`.
    
Main Paper material:

* `code/figures/figure_enloc_all_vs_eur.R` figure comparing ENLOC RCP when using all individuals vs using European only (Main Paper, suppl fig 24)
  
Companion paper material:

* `code/paper_material/tables.R` Generates latex tables to be included in the paper. At the moment: 
    * Supplementary Table S1: the list of 87 selected traits
    * Supplementary Table S2: expression and splicing models tally
* `code/figures/gwas_imputation_deflation.R` figure showing the deflation of GWAS' p-value distribution after imputation for 27 traits 
(Supplementary Figure S4)
* `code/figures/gwas_imputation_quality.R` scatter plot of original vs imputed GWAS zscores (Supplementary Figure S3)
* `code/figures/predixcan_enloc_eqtl_sqtl.R` plots summarising loci detection per emthod (Supplementary Figure S17, S18, S19, S20)
* `code/figures/models_gain.R` plots a comparison of numbers of models between Elastic Net models and MASHR-based models
* `code/figures/proportions_bundle.R` plots the proportion of enloc and s-predixcan detections (Supplementary Figure S14-A,B,C,D)
* `code/figures/upset_mashr_gwas_enloc_spredixcan.R` generates upset plots underlying Supplementary Figure S15 (S-PrediXcan/enloc loci with MASHR models)
* `code/figures/upset_mashr_gwas_enloc_spredixcan.R` generates upset plots underlying Supplementary Figure S16 (S-PrediXcan/enloc loci with Elastic Net models)
* For more details about figure please go to `code/`

Markdowns:

* `analysis/miscellaneous_statistics_2.Rmd`: This R markdown generates miscellaneous statistics and summaries from other methods results. 
i.e. This tallies summaries from S-PrediXcan results, enloc results; numbers of genomic loci with detections, etc.
* `analysis/gwas_enloc_predixcan_multixcan.Rmd`: analyzes GWAS, S-PrediXcan, S-MultiXcan, enloc results and builds upset plots.
This overlaps a bity with the previous markdown. Also Figure 5-e


A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr


# Deprecated

* `code/figures/FIG-DOSE-RESPONSE-CONCORDANCE-C.R` scatter plot of beta prim vs. sec for Whole Blood, Europeans, rcp>0.1 (Figure 3-C)
* `code/figures/SFIG-CONCORDANCE-MEDIATING-EFFECTS-RANK-BY-EFFECT-SIZE.R` scatter plot, residual plot, and p-value (Supplementary Figure S13)
