get_trait = function(str) {
  str_remove(str_remove(str, 'imputed_'), '.txt.gz')
}