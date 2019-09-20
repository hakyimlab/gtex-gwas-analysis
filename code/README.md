# Code

Save command-line scripts and shared R code here.

# Structure

`code/preprocess` contains scripts that were once for the project, to setup data and Google Cloud BigQuery tables.

`code/figures` contains scripts to generate specific figures for the paper.

`code/paper_material` contain miscellaneous scripts to build sumamry statistics, counts, tables, etc.

`code/checks` contains R scripts to explore and verify results; this is support material that didn't end in the paper(s)

`code/helpers` contains R functions reused across scripts.

## Preprocessing

`code/preprocess/preprocess_setup_auxiliary_tables.R` is a script left here for reference purposes only. It was used to create support tables in google cloud.

`code/preprocess_gwas_regions.R` is a script left here for reference purposes only. It was used to compute counts of detections from different methods in LD-independent regions.

## Figures

`code/figures/figure_enloc_all_vs_eur.R` yields a figure comparing ENLOC results when using all individuals, to ENLOC using only Europeans.

`code/figures/predixcan_enloc_eqtl_sqtl.R` yields many figures ummarising the number of loci with GWAS-significant detections, and enloc and S-MultiXcan results.
