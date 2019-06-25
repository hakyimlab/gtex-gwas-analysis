---
title: "Readme.rmd"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Miscellaneous scripts

All scripts in this folder are meant to be ran from this folder. Note that they assume a local `data` folder.

## Counting significant associations

The following is included for documentation reasons.
`count_hits.sh` is a PBS queue array job that will gather summary statistics (one trait per job) of significant results.
`count_hits.R` is the actual script being run in the above.
`count_hits_plots.R` picks up the output of the above and performs some plotting. Most of this code was migrated to a markdown under `/analysis` folder in this workflowr repository.

## S-MultiXcan and ENLOC summary plots.

`gtex_paper_enloc_multi_stats.sh` submits jobs to a PBS queue, that summarise S-MultiXcan and ENLOC results per trait.
`gtex_paper_enloc_multi_stats.R` is the R code executed in the above.
`gtex_paper_enloc_multi_stats_plot.R` can be ran after the previous, to build the plots. 
