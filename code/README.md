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

`code/figures/predixcan_enloc_eqtl_sqtl.R` yields many figures summarising the number of loci with GWAS-significant detections, and enloc and S-MultiXcan results.

`code/figures/mediation_analysis.R` yields two figures of mediation analysis: 1) correlation of |GWAS effect size| and |QTL effect size| (Fig2B); 2) gene-level effects estimated from mixed effect model (Fig2C). 

`code/figures/scatter_plot_script_whole_blood_rcp_0.1.R` yields two figures on plotting beta_gene's between primary and secondary QTL (Fig2E, SFig11)

`code/figures/PR-ROC-CURVES.R` yields figures, tables on silver standard analysis (ROC and PR curves, etc)

`code/figures/upset_top_gene_per_locus.R` yields upset figures for top gene per locus in silver standard analysis   

* Fig2B: `output/SFIG-COR-DELTA-GAMMA.png`
* Fig2C: `output/FIG-DOSE-RESPONSE-CONCORDANCE-B.png`
* Fig2E: `output/scatter_plot_script_whole_blood_rcp_0.1.png` 
* SFig11: `output/scatter_plot_script_whole_blood_rcp_0.1_by_effect_size.png`
* Fig3C/D: `output/ROC-*-OMIM.png`
* SFig25: `output/PR-*-OMIM.png`
* SFig26: `output/PR-*-EWAS.png`
* SFig24: `output/top_per_locus_upset_*.png`
* STab8: `output/AUC-and-ENRICH.tsv`

 
