library(tidyverse)

N_genes <- 1e4

beta <- rnorm(N_genes, 0, 1)
gamma <- rnorm(N_genes, 0, 1)
delta <- beta * gamma + rnorm(N_genes, 0, 1)

df <- data.frame("beta"=beta, "gamma"=gamma, "delta"=delta)

pp <- ggplot(df, aes(x=abs(gamma), y=abs(delta)))
pp <-  pp + geom_point(alpha=.2)
pp <-  pp + theme_bw(base_size=20)
pp <-  pp + xlab(expression(paste("|", gamma[gene],"|"))) + ylab(expression(paste("|", delta[GWAS], "|")))
pp <-  pp + ggtitle("Simulated GWAS and eQTL effect sizes")
pp
