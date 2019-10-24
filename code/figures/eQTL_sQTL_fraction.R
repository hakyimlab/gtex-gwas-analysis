library(dplyr)
library(ggplot2)
library(data.table)
theme_set(theme_bw(base_size = 25))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/8151c6fe70e3d4ee43d9ce340ecc0eb65172e616/my_ggplot_theme.R')
mymixer = c(
  'v8-sQTL' = rgb(200, 90, 40, maxColorValue = 255),
	'v8-sQTL-GWAS' = rgb(200, 90, 40, maxColorValue = 255),
  'v8-eQTL' = rgb(42, 155, 204, maxColorValue = 255),
  'v8-eQTL-GWAS' = rgb(42, 155, 204, maxColorValue = 255)
)
myshape = c(
  'v8-sQTL' = 19,
	'v8-sQTL-GWAS' = 15,
  'v8-eQTL' = 19,
  'v8-eQTL-GWAS' = 15
)

fontsize = 6

# ugly load in minimal data
minimal_data = readRDS('../../data/data_qtl_fraction_curve_from_boxiang/minimal_data.rds')
pct_vs_pval_threshold_v8_eQTL = minimal_data$pct_vs_pval_threshold_v8_eQTL
pct_vs_pval_threshold_v8_gwas_eQTL = minimal_data$pct_vs_pval_threshold_v8_gwas_eQTL
pct_vs_pval_threshold_v8_sQTL = minimal_data$pct_vs_pval_threshold_v8_sQTL
pct_vs_pval_threshold_v8_gwas_sQTL = minimal_data$pct_vs_pval_threshold_v8_gwas_sQTL
neglogq05_v8_sQTL = minimal_data$neglogq05_v8_sQTL
pct_snp_with_p_lt_05_v8_sQTL = minimal_data$pct_snp_with_p_lt_05_v8_sQTL
neglogq05_v8_eQTL = minimal_data$neglogq05_v8_eQTL
pct_snp_with_p_lt_05_v8_eQTL = minimal_data$pct_snp_with_p_lt_05_v8_eQTL
bonf_threshold_v8 = minimal_data$bonf_threshold_v8
pct_snp_with_p_lt_bonf_threshold_v8_eQTL = minimal_data$pct_snp_with_p_lt_bonf_threshold_v8_eQTL
pct_snp_with_p_lt_bonf_threshold_v8_gwas_eQTL = minimal_data$pct_snp_with_p_lt_bonf_threshold_v8_gwas_eQTL
pct_snp_with_p_lt_bonf_threshold_v8_sQTL = minimal_data$pct_snp_with_p_lt_bonf_threshold_v8_sQTL
pct_snp_with_p_lt_bonf_threshold_v8_gwas_sQTL = minimal_data$pct_snp_with_p_lt_bonf_threshold_v8_gwas_sQTL
pct_snp_with_p_lt_05_v8_gwas_eQTL = minimal_data$pct_snp_with_p_lt_05_v8_gwas_eQTL
pct_snp_with_p_lt_05_v8_gwas_sQTL = minimal_data$pct_snp_with_p_lt_05_v8_gwas_sQTL
fig_dir = '../../'


# plot_for_companion
message('plotting...')
pct_vs_pval_threshold_v8_sQTL$Version = 'v8-sQTL'
pct_vs_pval_threshold_v8_gwas_sQTL$Version = 'v8-sQTL-GWAS'

pct_vs_pval_threshold_v8_eQTL$Version = 'v8-eQTL'
pct_vs_pval_threshold_v8_gwas_eQTL$Version = 'v8-eQTL-GWAS'


### eQTL
pct_vs_pval_threshold = rbind(pct_vs_pval_threshold_v8_eQTL,pct_vs_pval_threshold_v8_gwas_eQTL)
pct_vs_pval_threshold[,Version := factor(Version,levels=c('v8-eQTL','v8-eQTL-GWAS'))]

p = ggplot(pct_vs_pval_threshold,aes(x=neglog10ub,y=cum_pct,color = Version, shape = Version))+
	geom_point(size = 3) +
	geom_line() +
	xlab(expression(paste(-log[10],'[',italic(P),']',sep=""))) +
	ylab('proportion of SNPs') + 
	# Setting scales:
	scale_x_continuous(breaks=c(0,5,10,20,neglogq05_v8_eQTL,30),labels=c(0,5,10,20,round(neglogq05_v8_eQTL,digits=2),30),expand=c(0.01,0))+
	scale_y_continuous(limits = c(NA,1.05), breaks=c(0,0.05,0.25,0.5,0.75,1.00),labels=c(0,0.05,0.25,0.5,0.75,1.00),expand=c(0,0.01)) +
	# Setting v8-eQTL arrows:
	geom_segment(aes(x = neglogq05_v8_eQTL, y = -0.02, xend = neglogq05_v8_eQTL, yend = 0.05),color='gray',linetype=2) +
	geom_segment(aes(x = -0.5, y = 0.05, xend = neglogq05_v8_eQTL, yend = 0.05),color='gray',linetype=2) +
	geom_segment(aes(x = -log10(0.05)+1,y = pct_snp_with_p_lt_05_v8_eQTL-0.025,xend = -log10(0.05),yend = pct_snp_with_p_lt_05_v8_eQTL),arrow=arrow(length=unit(0.02,'npc')),color='black')+
	annotate('text',x=-log10(0.05)+1.5,y=pct_snp_with_p_lt_05_v8_eQTL-0.025,label=paste0(100*round(pct_snp_with_p_lt_05_v8_eQTL,4),'% SNPs with P < 0.05'),hjust=0,size=fontsize)+
	geom_segment(aes(x=-log10(bonf_threshold_v8)+1,y=pct_snp_with_p_lt_bonf_threshold_v8_eQTL,xend=-log10(bonf_threshold_v8),yend=pct_snp_with_p_lt_bonf_threshold_v8_eQTL),arrow=arrow(length=unit(0.02,'npc')),color='black')+
	annotate('text',x=-log10(bonf_threshold_v8)+1.5,y=pct_snp_with_p_lt_bonf_threshold_v8_eQTL,label=paste0(100*round(pct_snp_with_p_lt_bonf_threshold_v8_eQTL,4),'% SNPs with P < 0.05/49'),hjust=0,size=fontsize) +
	# Setting legend:
	theme(legend.position = c(0.99, 0.6),legend.justification = c("right", "top"))
p = p + th
p = p + scale_color_manual(values = mymixer) + scale_shape_manual(values = myshape)
ggsave(paste0(fig_dir, 'output/eQTL_fraction.png'), p, height = 7, width = 7)



### sQTL
pct_vs_pval_threshold = rbind(pct_vs_pval_threshold_v8_sQTL,pct_vs_pval_threshold_v8_gwas_sQTL)
pct_vs_pval_threshold[,Version := factor(Version,levels=c('v8-sQTL','v8-sQTL-GWAS'))]

p2 = ggplot(pct_vs_pval_threshold,aes(x=neglog10ub,y=cum_pct,color = Version, shape = Version))+
  geom_point(size = 3) +
  geom_line() +
  xlab(expression(paste(-log[10],'[',italic(P),']',sep=""))) +
  ylab('proportion of SNPs') + 
  # Setting scales:
  scale_x_continuous(breaks=c(0,5,10,neglogq05_v8_sQTL,20,30),labels=c(0,5,10,round(neglogq05_v8_sQTL,digits=2),20,30),expand=c(0.01,0))+
  scale_y_continuous(limits = c(NA,1.05), breaks=c(0,0.05,0.25,0.5,0.75,1.00),labels=c(0,0.05,0.25,0.5,0.75,1.00),expand=c(0,0.01)) +
  # Setting v8-eQTL arrows:
  geom_segment(aes(x = neglogq05_v8_sQTL, y = -0.02, xend = neglogq05_v8_sQTL, yend = 0.05),color='gray',linetype=2) +
  geom_segment(aes(x = -0.5, y = 0.05, xend = neglogq05_v8_sQTL, yend = 0.05),color='gray',linetype=2) +
  geom_segment(aes(x = -log10(0.05)+1,y = pct_snp_with_p_lt_05_v8_sQTL-0.025,xend = -log10(0.05),yend = pct_snp_with_p_lt_05_v8_sQTL),arrow=arrow(length=unit(0.02,'npc')),color='black')+
  annotate('text',x=-log10(0.05)+1.5,y=pct_snp_with_p_lt_05_v8_sQTL-0.025,label=paste0(100*round(pct_snp_with_p_lt_05_v8_sQTL,4),'% SNPs with P < 0.05'),hjust=0,size=fontsize)+
  geom_segment(aes(x=-log10(bonf_threshold_v8)+1,y=pct_snp_with_p_lt_bonf_threshold_v8_sQTL,xend=-log10(bonf_threshold_v8),yend=pct_snp_with_p_lt_bonf_threshold_v8_sQTL),arrow=arrow(length=unit(0.02,'npc')),color='black')+
  annotate('text',x=-log10(bonf_threshold_v8)+1.5,y=pct_snp_with_p_lt_bonf_threshold_v8_sQTL,label=paste0(100*round(pct_snp_with_p_lt_bonf_threshold_v8_sQTL,4),'% SNPs with P < 0.05/49'),hjust=0,size=fontsize) +
  # Setting legend:
  theme(legend.position = c(0.99, 0.6),legend.justification = c("right", "top"))

p2 = p2 + th
p2 = p2 + scale_color_manual(values = mymixer) + scale_shape_manual(values = myshape)
ggsave(paste0(fig_dir, 'output/sQTL_fraction.png'), p2, height = 7, width = 7)