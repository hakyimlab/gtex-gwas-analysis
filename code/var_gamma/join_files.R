library(dplyr)
library(stringr)

rows <- vector(mode="list", length=0)
dir <- "var_gammas"

for (file in list.files(dir)) {
  df <- read.table(file.path(dir, file), header=TRUE)
  tissue <- str_match(pattern="(.*)_var_gamma.txt", string=file)[2]
  df$tissue <- tissue 
  rows <- c(rows, list(df))
}

df <- bind_rows(rows)

write.table(df, row.names=FALSE, quote=FALSE, "kk.txt", sep="\t")
