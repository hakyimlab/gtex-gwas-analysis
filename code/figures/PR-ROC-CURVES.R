# This scripts plot figures going to main panel and also supplement for silver standard analysis
# They are (by TEX label)
# FIG-IDENTIFYING-TARGET-GENES panel B  : ROC for expression PrediXcan-mashr and enloc (OMIM) 
# FIG-IDENTIFYING-TARGET-GENES panel C  : ROC for splicing PrediXcan-mashr and enloc (OMIM)
# SFIG-IDENTIFYING-TARGET-GENES all panels  : PR & joint PR for expression & splicing PrediXcan-mashr and enloc (OMIM)
# SFIG-IDENTIFYING-TARGET-GENES-EWAS all panels : PR & joint PR for expression & splicing PrediXcan-mashr and enloc (EWAS)
# Additional figures (not in paper yet): log odds ratio of logistic regression based test for expression & splicing PrediXcan-mashr and enloc (OMIM & EWAS)

# r set up
library(ggplot2)
library(dplyr)
library(reshape2)
theme_set(theme_bw(base_size=30))
outdir = '../../output'
# set up color guide here
cbPalette = c('enloc' = "#23B509", 'smr' = "#D3BE0D", "PrediXcan-mashr" = "#6209B5", 'coloc' = '#FA12E8' )
# , 'multixcan' = "#000000", 'predixcan' = "#E69F00", 'coloc' = "#009E73", "enloc-predixcan" = "#0072B2", 'predixcan-mashr' = "#F0E442", "predixcan-en-dapgw" = "#0072B2", "random" = "#808080", "predixcan-en-dapgw-eur" = "#00FF00"
myline = c('enloc' = 'solid', "PrediXcan-mashr" = 'solid', 'smr' = 'solid', 'coloc' = 'solid')
# set up ylim for PR plots
ylim_list = list(EWAS = 0.6, OMIM = 0.4)

source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/8151c6fe70e3d4ee43d9ce340ecc0eb65172e616/my_ggplot_theme.R')

compute_auc = function(curve) {
  # f = curve %>% filter(nTP > 0) %>% group_by(nTP) %>% summarize(y = mean(tpr), x = fpr[1] - fpr[2])
  o1 = rbind(c(1, 1), curve %>% filter(nTP > 0) %>% select(tpr, fpr))
  o2 = rbind(curve %>% filter(nTP > 0) %>% select(tpr, fpr), c(0, 0))
  # message('e')
  y = data.frame(x1 = o1$fpr, x2 = o2$fpr, y1 = o1$tpr, y2 = o2$tpr) 
  e = y %>% mutate(s = (x1 - x2) * (y1 + y2) / 2)
  data.frame(roc_auc = sum(e$s))
}

th2 = th # + theme(axis.text = element_blank(), axis.ticks = element_blank())

# load data
data_folder = '../../data/data_silver_standard'
subfolders = c('expression', 'splicing')
datasets = c('OMIM', 'EWAS')
## load PR/ROC curves
curve_data_list = list()
for(i in subfolders) {
  for(j in datasets) {
    input_name = paste0(data_folder, '/', i, '/data-for-plot.', j, '-LD-block.rds')
    curve_data_list[[i]][[j]] = readRDS(input_name)
  }
}
## load logistic estimates
logistic_data_list = list()
for(i in subfolders) {
  for(j in datasets) {
    input_name = paste0(data_folder, '/', i, '/logistic-based-ad-hoc-plot-data.', j, '-LD-block-PrediXcan-MASH-EUR.rds')
    logistic_data_list[[i]][[j]] = readRDS(input_name)
  }
}
## load data to compute enrichment
enrich_data_list = list()
for(i in subfolders) {
  for(j in datasets) {
    input_name = paste0(data_folder, '/', i, '/logistic-based-test.datamatrix.', j, '-LD-block-PrediXcan-MASH-EUR.tsv')
    enrich_data_list[[i]][[j]] = read.table(input_name, sep = '\t', header = T)
  }
}


# plot ROC
df_auc = data.frame()
for(i in subfolders) {
  for(j in datasets) {
    # load pre-computed ROC curve
    margin_curve2 = curve_data_list[[i]][[j]]$margin_curve2
    margin_curve2$method[margin_curve2$method == 'predixcan-mashr-eur'] = 'PrediXcan-mashr'
    p = ggplot() + geom_path(data = margin_curve2 %>% filter(method != 'coloc'), aes(fpr, tpr, color = method, linetype = method), size = 2, alpha = .6)
    # p = p + geom_point(data = point_curve2, aes(fpr, tpr, shape = method))
    p = p + scale_color_manual(values = cbPalette) + coord_cartesian(ylim = c(0, 1), xlim = c(0, 1))
    p = p + geom_abline(slope = 1, intercept = 0, color = 'lightgray') + theme(legend.position = c(0.67, 0.2)) + scale_linetype_manual(values = myline)
    p = p + th2
    p = p + xlab('False positive rate') + ylab('True positive rate')
    tmp = margin_curve2 %>% group_by(method) %>% do(compute_auc(.)) %>% ungroup()
    df_auc = rbind(df_auc, cbind(
      data.frame(
        regulation = i,
        dataset = j
      ),
      tmp
    )
    )
    ggsave(paste0(outdir, '/', 'ROC-', i, '-', j, '.png'), p + theme(aspect.ratio = 1))
  }
}
write.table(df_auc, paste0(outdir, '/', 'ROC-AUC.tsv'), quo = F, sep = '\t', row = F, col = T)


# plot PR 
for(i in subfolders) {
  for(j in datasets) {
    # load pre-computed PR curve
    margin_curve = curve_data_list[[i]][[j]]$margin_curve
    margin_curve$method[margin_curve$method == 'predixcan-mashr-eur'] = 'PrediXcan-mashr'
    p = ggplot() + geom_path(data = margin_curve %>% filter(method != 'coloc'), aes(recall, precision, color = method, linetype = method), size = 2, alpha = .6)
    p = p + scale_color_manual(values = cbPalette) + coord_cartesian(ylim = c(0, ylim_list[[j]]), xlim = c(0, 1)) + theme(legend.position = c(0.7, 0.8)) + scale_linetype_manual(values = myline)
    p = p + th
    p = p + xlab('Recall') + ylab('Precision')
    ggsave(paste0(outdir, '/', 'PR-marginal-', i, '-', j, '.png'), p + theme(aspect.ratio = 1))
  }
}

# plot PR joint with enloc
for(i in subfolders) {
  for(j in datasets) {
    # load pre-computed PR and PR-joint curves
    margin_curve = curve_data_list[[i]][[j]]$margin_curve
    joint_curve_g = curve_data_list[[i]][[j]]$joint_curve_g
    margin_curve$method[margin_curve$method == 'predixcan-mashr-eur'] = 'PrediXcan-mashr'
    joint_curve_g$method[joint_curve_g$method == 'predixcan-mashr-eur'] = 'PrediXcan-mashr'
    p = ggplot() + geom_path(data = margin_curve %>% filter(method == 'enloc'), aes(recall, precision, color = method), size = 2, alpha = .6)
    p = p + geom_path(data = joint_curve_g %>% filter(method != 'coloc'), aes(recall, precision, color = method), size = 2, alpha = .6)
    p = p + scale_color_manual(values = cbPalette) + coord_cartesian(ylim = c(0, ylim_list[[j]]), xlim = c(0, 1)) + theme(legend.position = c(0.7, 0.8)) + scale_linetype_manual(values = myline)
    p = p + th
    p = p + xlab('Recall') + ylab('Precision')
    ggsave(paste0(outdir, '/', 'PR-joint-', i, '-', j, '.png'), p + theme(aspect.ratio = 1))
  }
}

# plot enloc vs coloc
for(i in subfolders) {
  for(j in datasets) {
    # load pre-computed PR and PR-joint curves
    margin_curve = curve_data_list[[i]][[j]]$margin_curve
    margin_curve$method[margin_curve$method == 'predixcan-mashr-eur'] = 'PrediXcan-mashr'
    if('coloc' %in% margin_curve$method) {
      p = ggplot() + geom_path(data = margin_curve %>% filter(method %in% c('enloc', 'coloc')), aes(recall, precision, color = method, linetype = method), size = 2, alpha = .6)
      p = p + scale_color_manual(values = cbPalette) + coord_cartesian(ylim = c(0, ylim_list[[j]]), xlim = c(0, 1)) + theme(legend.position = c(0.7, 0.8)) + scale_linetype_manual(values = myline)
      p = p + th
      p = p + xlab('Recall') + ylab('Precision')
      ggsave(paste0(outdir, '/', 'PR-coloc-', i, '-', j, '.png'), p + theme(aspect.ratio = 1)) 
    }
  }
}


# plot estimates from logistic test
for(i in subfolders) {
  for(j in datasets) {
    # load pre-computed estimates 
    df_out = logistic_data_list[[i]][[j]]
    p = df_out %>% filter(type == 'rank_based', per_locus_method == 'rank_based', scoring_method == 'without_score', pvalue_scale == '-log10', variable != '(Intercept)') %>% ggplot() + geom_point(aes(x = variable, y = Estimate), size = 5) + geom_errorbar(aes(x = variable, ymin = Estimate - 1.96 * `Std. Error`, ymax = Estimate + 1.96 * `Std. Error`), width = .2) + geom_hline(yintercept = 0, linetype = 2) + ylab('log odds ratio') + xlab('method')  + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
    p = p + th
    ggsave(paste0(outdir, '/', 'logistic-estimates-', i, '-', j, '.png'), p)
  }
}

# for enrichment
enloc_cutoff = 0.5
df_enrich = data.frame()
for(i in subfolders) {
  for(j in datasets) {
    score = enrich_data_list[[i]][[j]]
    pval_cutoff_nlog10 = -log10(0.05 / nrow(score) / 49) 
    score = score %>% mutate(is_predixcan_signif = predixcan_mashr_eur_score > pval_cutoff_nlog10, is_enloc_signif = enloc_score > enloc_cutoff)  # , is_smr_signif = smr_score > pval_cutoff_nlog10)
    tmp = score %>% filter(!duplicated(paste(gene, trait))) %>% group_by(is_omim) %>% summarize(`PrediXcan-mashr` = mean(is_predixcan_signif), enloc = mean(is_enloc_signif))  # , smr = mean(is_smr_score)) 
    df_enrich = rbind(df_enrich, 
      rbind(tmp, (tmp[ tmp$is_omim, ] / tmp[ !tmp$is_omim, ]) %>% as.data.frame %>% mutate(is_omim = 'fold')) %>% mutate(regulation = i, dataset = j)
    )
  }
}
write.table(df_enrich, paste0(outdir, '/', 'Top-per-locus-enrich.tsv'), quo = F, sep = '\t', row = F, col = T)


# combining enrichment and auc
tmp = df_enrich %>% filter(is_omim == 'fold')
tmp = tmp %>% select(-is_omim) %>% melt(id.var = c('regulation', 'dataset'))
tmp = tmp %>% rename(method = variable, enrich_fold = value)
dfm = left_join(df_auc, tmp, by = c('regulation', 'dataset', 'method'))
write.table(dfm, paste0(outdir, '/', 'AUC-and-ENRICH.tsv'), quo = F, sep = '\t', row = F, col = T)