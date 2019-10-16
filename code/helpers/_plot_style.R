
paper_base_theme_ <- theme(panel.border = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "#170a45", size = .5),
                axis.ticks = element_line(colour = "#170a45", size = .2),
                axis.text = element_text(color = '#170a45'))

method_palette_ <- c('enloc' = "#23B509",
                     'smr' = "#D3BE0D",
                     "PrediXcan-mashr" = "#6209B5",
                     'coloc' = '#FA12E8' )
data_palette_ <- c(
  'splicing' = rgb(200, 90, 40, maxColorValue = 255),
  'expression' = rgb(42, 155, 204, maxColorValue = 255)
)
