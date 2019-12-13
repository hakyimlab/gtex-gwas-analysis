# verifying harmonization+imputation on Rapid GWAS

The scripts in this folder were run to obtain a harmonized, summary-statistics-imputed GWAS
from the rapid GWAS project.
The overarching goal was the assessment of sign mismatches in ambiguous snps for this numerous family of GWAS.
The product is in this repository:
`data/new_ukb_gwas/imputed_50_standing_height.txt.gz`

It was built with the following in CRI gardner:

0) Standing Height GWAS was downloaded via the following:
```
wget https://www.dropbox.com/s/od6dr8kdrrornuz/50_raw.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 50_standing_height.txt.gz
```
1) `ukb_harmonization.sh` was run to harmonize the input GWAS
2) The scripts in `jobs_summary_imputation` were run to generate imputed summary statistics
3) `collect_imputed_50_standing_height.sh` was run to build the properly formmated, harmonized+imputed GWAS
