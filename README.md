# Resources shared by manuscript 
["Widespread dose-dependent effects of RNA expression and splicing on complex diseases and traits"](https://doi.org/10.1101/814350) data and resources


currently being uploaded to Zenodo

- colocalization results
- predicted gene association results
- fine-mapping 


check soon...





# gtex-gwas-analysis

In order to run these scripts, the following R packages are required:
- bigrquery
- tidyverse
- upsetR

To use BigQuery, you need to install [Google Cloud SDK](https://cloud.google.com/sdk/) and have read access to the required Google Cloud tables.

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

* `code/paper_material/tables.R` Generates latex tables to be included in the paper 
(at the moment, only one table: Supplementary Table 1: the list of 87 selected traits)
* `code/figures/gwas_imputation_deflation.R` figure showing the deflation of GWAS' p-value distribution after imputation for 27 traits 
(Supplementary Figure S3)
* `code/figures/gwas_imputation_quality.R` scatter plot of original vs imputed GWAS zscores (Supplementary Figure S2)
* `code/figures/FIG-DOSE-RESPONSE-CONCORDANCE-C.R` scatter plot of beta prim vs. sec for Whole Blood, Europeans, rcp>0.1 (Figure 3-C)
* `code/figures/SFIG-CONCORDANCE-MEDIATING-EFFECTS-RANK-BY-EFFECT-SIZE.R` scatter plot, residual plot, and p-value (Supplementary Figure S13)
* `code/figures/predixcan_enloc_eqtl_sqtl.R` plots summarising loci detection per emthod (Supplementary Figure S12, S13, S14, S15)

Markdowns:

* `analysis/miscellaneous_statistics_2.Rmd`: This R markdown generates miscellaneous statistics and summaries from other methods results. 
i.e. This computes summaries from S-PrediXcan results, enloc results; numbers of genomic loci with detections, etc.
* `analysis/gwas_enloc_predixcan_multixcan.Rmd`: analyzes GWAS, S-PrediXcan, S-MultiXcan, enloc results and builds upset plots.
This overlaps a bity with the previous markdown. Also Figure 5-e



A [workflowr][] project.

[workflowr]: https://github.com/jdblischak/workflowr
