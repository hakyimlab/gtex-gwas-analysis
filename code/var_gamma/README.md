These scripts are meant to be run on CRI.
First execute the following command:

```bash
$bash compute_var_gamma.sh
```

It will first generate a set of eQTL files (one per tissue) containing only two columns, `slope` and `slope_se`.
Then, an R script will create a file per tissue, with a single row and two columns: 1)variance of the gamma coefficients and 2) the mean squared standard errors.
Then execute:

```bash
$Rscript join_files.R
```
