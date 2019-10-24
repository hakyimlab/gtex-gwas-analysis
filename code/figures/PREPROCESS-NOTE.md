# `mediation_analysis.R` and `scatter_plot_script_whole_blood_rcp_0.1.R`

Preprocessing pipeline is at [here](https://bitbucket.org/yanyul/rotation-at-imlab/src/master/analysis/reproduce_mediation_analysis/) with job config files at `qsub`/

# `PR-ROC-CURVES.R`

Preprocessing pipeline is at [here](https://bitbucket.org/yanyul/rotation-at-imlab/src/master/analysis/fdr_power_specificity/)
with job config files at `qsub/config.ldblock-*.yaml`

# `causal_tissue.R`

TBA

# `eQTL_sQTL_fraction.R`

Initial contribution by Boxiang Liu. 
Original scripts and data are shared at [here](https://drive.google.com/open?id=1wNhfUgQ1MDHSZJjXDNmsEkWpfWlAPtfi)
The data is in `processed_data`, and the scripts are in scripts
For this analysis, set scripts as the home directory, and run `Rscript 006_eQTL_and_sQTL.R`