library(tidyverse)
v6mcalls <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/chrX.methylation.tsv")


v6mcalls$pos <- rep("chrx", nrow(v6mcalls))
v6mcalls[v6mcalls$start >= 57828561 & v6mcalls$end  <= 60664792 , ][, "pos"] <- "cenx"
v6mcalls[v6mcalls$start >=  113868842 & v6mcalls$end <= 114116851,][, "pos"] <- "dxz4" 

ggplot(data = v6mcalls, aes(x = log_lik_ratio, fill=v6mcalls$pos)) +geom_histogram(binwidth=.5,alpha=.5)+facet_wrap(~v6mcalls$pos, ncol=1, scales = "free")+theme_bw()+scale_fill_manual("Region",values=c( "red", "green", "blue"))+xlim(25,-25)

v6_cen <- v6mcalls %>%
  filter(pos=="cenx") %>%
  ggplot(aes(x=log_lik_ratio))+geom_histogram(binwidth=.5)+xlim(25,-25)+theme_bw()+labs(title="chm13 v6 assembly cenx")

v6_dxz4 <-v6mcalls %>%
  filter(pos=="cenx") %>%
  ggplot(aes(x=log_lik_ratio))+geom_histogram(binwidth=.5)+xlim(25,-25)+theme_bw()+labs(title="chm13 v6 assembly dxz4")

v6_chrx <- v6mcalls %>%
  filter(pos=="cenx") %>%
  ggplot(aes(x=log_lik_ratio))+geom_histogram(binwidth=.5)+xlim(25,-25)+theme_bw()+labs(title="chm13 v6 assembly chrx")

  
v4mcalls <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v4_assembly/methylation.tsv")

v4mcalls %>%
  ggplot(aes(x=log_lik_ratio))+geom_histogram(binwidth=.5)+xlim(25,-25)+annotate("text", x = 10, y = 2000000, label = "average methylation : 29.9%")+annotate("text", x = 10, y = 2100000, label = "called sites vs all sites : 58%")+labs(title="chm13 v4 assembly chrX")+theme_bw()

gm128 <- read_tsv("/dilithium/Data/Nanopore/projects/gm12878/analysis/mcall-cpg/GM12878.46.cpg.meth.tsv")

gm128 %>%
  ggplot(aes(x=log_lik_ratio))+geom_histogram(binwidth=.5)+xlim(25,-25)+annotate("text", x = 10, y = 230000, label = "average methylation : 3.87e-05%")+annotate("text", x = 10, y = 250000, label = "called sites vs all sites : 100%")+labs(title="gm12878 wgs subset")+theme_bw()

nome <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/nanonome/GM12878_nanoNOMe_subset.cpg.meth.tsv")

nome %>%
  ggplot(aes(x=log_lik_ratio))+geom_histogram(binwidth=.5)+xlim(25,-25)+annotate("text", x = 10, y = 1000, label = "average methylation : 46%")+annotate("text", x = 10, y = 1100, label = "called sites vs all sites : 49%")+labs(title="gm12878 nanonome subset")+theme_bw()
