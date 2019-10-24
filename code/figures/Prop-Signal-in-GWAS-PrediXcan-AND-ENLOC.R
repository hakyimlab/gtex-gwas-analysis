# generate minimal data for plot

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(ggplot2))
# suppressPackageStartupMessages(library(UpSetR))
suppressPackageStartupMessages(library(viridis))
suppressWarnings(source("https://raw.githubusercontent.com/hakyimlab/gtex-gwas-analysis/master/code/helpers/_helpers.R?token=AC7RPMML4BAW75VV6Q66IWC5WBGSO"))


SUMMARY <- "https://raw.githubusercontent.com/hakyimlab/gtex-gwas-analysis/master/data/summaries/mashr_regions.txt?token=AC7RPMJPIUICXGCQIMJHEL25WBGWQ"

gwas_metadata <- "https://raw.githubusercontent.com/hakyimlab/gtex-gwas-analysis/master/data/gwas_metadata.txt?token=AC7RPMP2LDH6ILLAUMITC3S5WBGYM" %>% read_tsv

RESULT<-"output/g"
dir.create(RESULT, showWarnings = FALSE, recursive=TRUE)
fp_ <- function(p) file.path(RESULT, p)


ldblockinfo = read_tsv(SUMMARY) %>% rename(trait=phenotype, n=count, results=method)

## spread data 
tempo <- ldblockinfo %>% mutate(flag = n >= 1) %>% rename(results_type=results) %>% 
  select(region, trait, results_type, flag) %>% spread(key=unique(results_type),value=flag)
## set NAs to 0
tempo[is.na(tempo)] <- 0
tempo$gwas <- 1

tempo <- tempo %>%
  mutate(spredixcan_eqtl_enloc = ifelse(spredixcan_eqtl & enloc_eqtl, 1, 0)) %>%
  mutate(spredixcan_sqtl_enloc = ifelse(spredixcan_sqtl & enloc_sqtl, 1, 0)) %>%
  mutate(smultixcan_eqtl_enloc = ifelse(smultixcan_eqtl & enloc_eqtl, 1, 0)) %>%
  mutate(smultixcan_sqtl_enloc = ifelse(smultixcan_sqtl & enloc_sqtl, 1, 0))

tempo <- tempo %>% mutate(enloc_eqtl_sqtl = ifelse(enloc_eqtl & enloc_sqtl, 1, 0))
summary <- apply(tempo %>% select(-region,-trait),2,sum)

d_ <- tempo%>% group_by(trait) %>%
  summarise(n=n(), 
            spredixcan_eqtl=sum(spredixcan_eqtl), spredixcan_sqtl=sum(spredixcan_sqtl),
            enloc_eqtl=sum(enloc_eqtl), enloc_sqtl=sum(enloc_sqtl),
            spredixcan_enloc_eqtl=sum(spredixcan_eqtl_enloc), spredixcan_enloc_sqtl=sum(spredixcan_sqtl_enloc)) %>%
  gather(key="method_", value="count", -trait, -n) %>%
  mutate(f = count/n) %>% 
  mutate(data = ifelse(grepl("eqtl", method_), "expression", "splicing")) %>%
  mutate(method = ifelse(grepl("^spredixcan_(\\w)qtl$", method_), "S-PrediXcan",
                         ifelse(grepl("^enloc_(\\w)qtl$", method_),"ENLOC", "S-PrediXcan\n&\nENLOC")))

# original plot
# p_ <- ggplot(d_) + theme_bw(base_size=18) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   geom_boxplot(aes(method, f, fill=data)) +
#   ggtitle("Proportion of GWAS loci\nwith colocalized/associated genes or introns")  +
#   ylab("Proportion of loci")
# p_

# new plot
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/8151c6fe70e3d4ee43d9ce340ecc0eb65172e616/my_ggplot_theme.R')
mymixer = c(
  'splicing' = rgb(200, 90, 40, maxColorValue = 255),
  'expression' = rgb(42, 155, 204, maxColorValue = 255)
)
p = ggplot(d_) + theme_bw(base_size=25) + 
  th +
  geom_violin(aes(method, f, fill=data), position = position_dodge(width = .8)) +
  geom_boxplot(aes(method, f, fill=data), position = position_dodge(width = .8), width = .1) +
  scale_fill_manual(values = mymixer) +
  ylab("Proportion of loci")
ggsave('../../output/prop_gwas_loci_coloc_assoc.png', p, width = 9, height = 7)
