filter_palindromic <- function(d) {
  d %>% filter(is_palindromic(gtex_variant_id))
                 
}

is_palindromic <- function(x) {
  grepl("A_T", x) |
    grepl("T_A", x) |
    grepl("C_G", x) |
    grepl("G_C", x)
}

filter_not_biallelic <- function(d) {
  g <- d %>% group_by(chromosome, position) %>% summarise(n=n()) %>% filter(n>1) %>% select(chromosome, position)
  d %>% inner_join(g, by=c("chromosome", "position"))
}

round_any = function(x, accuracy, f=round){
  f(x/ accuracy) * accuracy
}