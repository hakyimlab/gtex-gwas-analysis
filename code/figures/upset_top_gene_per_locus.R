exp = read.table('../../data/data_silver_standard/expression/logistic-based-test.datamatrix.OMIM-LD-block-PrediXcan-MASH-EUR.tsv', header = T, sep = '\t')
spl = read.table('../../data/data_silver_standard/splicing/logistic-based-test.datamatrix.OMIM-LD-block-PrediXcan-MASH-EUR.tsv', header = T, sep = '\t')

is_top = function(is_omim, myrank) {
  e = myrank[is_omim]
  (0 %in% e) * 1
}

library(dplyr)
library(UpSetR)
mat = exp %>% group_by(lead_var, trait) %>% summarize(total = 1, top_proximity = is_top(is_omim, rank_proximity), top_predixcan = is_top(is_omim, predixcan_mashr_eur_rank), top_enloc = is_top(is_omim, enloc_rank)) %>% ungroup()
 
mat2 = spl %>% group_by(lead_var, trait) %>% summarize(total = 1, top_proximity = is_top(is_omim, rank_proximity), top_predixcan = is_top(is_omim, predixcan_mashr_eur_rank), top_enloc = is_top(is_omim, enloc_rank)) %>% ungroup()



splicing_col = rgb(200, 90, 40, maxColorValue = 255)
expression_col = rgb(42, 155, 204, maxColorValue = 255)

png('../../output/top_per_locus_upset_expr.png', width = 6, height = 4, units = 'in', res = 300)
upset(as.data.frame(mat), order.by = 'freq', main.bar.color = expression_col, point.size = 5, line.size = 1, text.scale = 1.5)
dev.off()

png('../../output/top_per_locus_upset_spli.png', width = 6, height = 4, units = 'in', res = 300)
upset(as.data.frame(mat2), order.by = 'freq', main.bar.color = splicing_col, point.size = 5, line.size = 1, text.scale = 1.5)
dev.off()