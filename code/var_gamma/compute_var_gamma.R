library(argparse)
library(stringr)

# READ ARGUMENTS
parser <- ArgumentParser()
parser$add_argument("--file")
parser$add_argument("--preprocessed_files_dir")
parser$add_argument("--output_dir")
args <- parser$parse_args()

tissue <- str_match(pattern="(.*).allpairs_beta_and_se.txt", string=args$file)[2]

if (!dir.exists(args$output_dir))
  dir.create(args$output_dir)

df <- read.table(file.path(args$preprocessed_files_dir, args$file), header=TRUE)

row <- data.frame("var_gamma"=var(df$slope), "var_se"=sum(df$slope_se**2, na.rm=TRUE)/nrow(df))
write.table(row, glue::glue("{output_dir}/{tissue}_var_gamma.txt"), quote=FALSE, row.names=FALSE, sep="\t")
