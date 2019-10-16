
remove_id_from_ensemble <- function(v) {
  k <- gsub("\\.(.*)", "", v)
  return(k)
}

fill_with_zeros <- function(d) {
  d %>% mutate_all(funs(ifelse(is.na(.), 0, .)))
}

file_logic_ <- function(folder, pattern) {
  names <- list.files(folder)
  names <- names[grepl(pattern, names)] %>% sort
  data.frame(name = gsub(pattern, "\\1", names),
             path = file.path(folder, names),
             stringsAsFactors = FALSE)
}

g_results_logic <- function(folder, pattern) {
  file_logic_(folder, pattern) %>% rename(trait = name)
}

m_results_logic <- function(folder, pattern) {
  file_logic_(folder, pattern) %>% rename(tissue = name)
}

d_results_logic_ <- function(folder, pattern, strip) {
  names <- list.files(folder)
  names <- names[grepl(pattern, names)] %>% sort
  data.frame(name = gsub(strip, "", names),
             path = file.path(folder,names),
             name_1 = gsub(pattern, "\\1", names),
             name_2 = gsub(pattern, "\\2", names),
             stringsAsFactors = FALSE)

}

cpt_results_logic <- function(folder, pattern, strip) {
  d_results_logic_(folder, pattern, strip) %>% rename(tissue = name_1, pheno=name_2)
}

p_results_logic <- function(folder, pattern, strip) {
  d_results_logic_(folder, pattern, strip) %>% rename(pheno = name_1, tissue=name_2)
}

load_from_logic_ <- function(logic, name_) {
  logic %>% filter(name == name_) %>% .$path %>% r_tsv_
}


load_from_logic_2_ <- function(logic, pheno_, tissue_) {
  logic %>% filter(pheno == pheno, tissue == tissue_) %>% .$path %>% r_tsv_
}

db_ <- function(path, callback) {
  con <- dbConnect(SQLite(), path)
  results <- tryCatch({
    callback(con)
  }, error = function(e) {
    print(e)
    NULL
  }, finally = {
    dbDisconnect(con)
  })
  results
}

r_tsv_ <- function(path, col_types=NULL, col_names=TRUE) { suppressMessages(read_tsv(path, col_types=col_types, col_names=col_names)) }

r_csv_ <- function(path, col_types=NULL, col_names=TRUE) { suppressMessages(read_csv(path, col_types=col_types, col_names=col_names)) }

save_plot <- function(plot, path, height, width, res=NA) {
  png(path, height=height, width=width, res=res)
  print(plot)
  dev.off()
}

save_delim <- function(x, path) {
  if (grepl(".gz$",path)) {
    gz1 <- gzfile(path, "w")
    write.table(x, gz1, sep="\t", row.names = FALSE, quote=FALSE)
    close(gz1)
  } else {
    write.table(x, path, sep="\t", row.names = FALSE, quote=FALSE)
  }
}

d_theme_ <- function() {
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5, face="bold", size=27),
        plot.subtitle = element_text(hjust=0.5, face="italic", size=25),
        axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size = 15)) +
  ggplot2::theme(legend.position="bottom",legend.direction="horizontal")
}

load_folder <- function(path, filter_=NULL, tag_=NULL, col_types=NULL) {
  files <- list.files(path)
  if (!is.null(filter_)) {
    files <- files[filter_(files)]
  }
  files <- sort(files)
  r <- list()
  for (i in 1:length(files)) {
    #message(files[i], " ", i)
    d <- file.path(path, files[i]) %>% r_tsv_(col_types=col_types)
    if (!is.null(tag_))
      d <- d %>% tag_(files[i])
    r[[i]] <- d
  }
  do.call(rbind, r)
}

read_or_get <- function(path, method) {
  if (!file.exists(path)) {
    d <- method()
    d %>% save_delim(path)
    d
  } else {
    path %>% r_tsv_
  }
}
