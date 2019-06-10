library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)

# read in gm12878 data
gm128_meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/data/nanopolish_calls/GM12878.chrX.cpg.meth.tsv", col_names = T) %>%
  add_column(strand="+", .before = 2) %>%
  rename(num_motifs=num_recsites)

# annotate based on grch38 coordinates -- likely not accurate for cenx 
gm128_meth$pos <- rep("chrx", nrow(gm128_meth))
gm128_meth[gm128_meth$start >= 58605580 & gm128_meth$end  <= 62412542 , ][, "pos"] <- "cenx"
gm128_meth[gm128_meth$start >=  115840395 & gm128_meth$end <= 115969110,][, "pos"] <- "dxz4" 
gm128_meth$sample <- "gm12878"

# read in meth frequency file
gm128_freq <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/data/nanopolish_calls/chrX.methylation_frequency.tsv")

gm128_freq$pos <- rep("chrx", nrow(gm128_freq))
gm128_freq[gm128_freq$start >= 58605580 & gm128_freq$end  <= 62412542 , ][, "pos"] <- "cenx"
gm128_freq[gm128_freq$start >=  115840395 & gm128_freq$end <= 115969110,][, "pos"] <- "dxz4" 
gm128_freq$sample <- "gm12878"

# read in chm13 data to compate meth info
v6mcalls <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/chrX.methylation.tsv")


v6mcalls$pos <- rep("chrx", nrow(v6mcalls))
v6mcalls[v6mcalls$start >= 57828561 & v6mcalls$end  <= 60664792 , ][, "pos"] <- "cenx"
v6mcalls[v6mcalls$start >=  113868842 & v6mcalls$end <= 114116851,][, "pos"] <- "dxz4"
v6mcalls$sample <- "chm13"

meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/chrX.methylation_frequency.tsv")
meth$pos <- rep("chrx", nrow(meth))
meth[meth$start >= 57828561 & meth$end  <= 60664792 , ][, "pos"] <- "cenx"
meth[meth$start >=  113868842 & meth$end <= 114116851,][, "pos"] <- "dxz4"
meth$sample <- "chm13"

# merge gm128 and chm13 data sets
all_freq <- rbind(meth, gm128_freq)
all_meth <- rbind(v6mcalls, gm128_meth)

ggplot(data = all_freq, aes(x = methylated_frequency, fill=all_freq$sample)) +geom_density(alpha=.5)+facet_wrap(~all_freq$pos, ncol=2)+theme_bw()+labs(x="methylated frequency (per base)", color="Position")+scale_fill_manual("Assembly",values=c( "red", "green", "blue"))

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/plots/frequency_perbase.pdf", height = 10, width = 10)

binned <- all_freq %>%
  mutate(bin = ntile(start, (length(meth$start)/50))) %>%
  group_by(bin, pos, sample) %>%
  summarise(avgmeth = mean(methylated_frequency)) %>%
  arrange(desc(avgmeth))

ggplot((data=binned),aes(x=avgmeth, fill = sample))+geom_density(alpha=.5)+theme_bw()+labs(x="Average methylation frequency (bin size 50bps)")+facet_wrap(~pos, ncol=2)+scale_fill_manual("Assembly",values=c( "red", "green", "blue"))

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/plots/frequency_binned.pdf", height = 10, width = 10)

ggplot(data = all_meth, aes(x=log_lik_ratio, fill = all_meth$sample))+geom_density(alpha=.5)+xlim(25,-25)+theme_bw()+scale_fill_manual("Assembly",values=c( "red", "green"))+geom_vline(xintercept= c(-2.5, 2.5), linetype="dashed", color="black")+facet_wrap(~pos, ncol=2)

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/plots/log_lik_region.pdf", height = 10, width = 10)

# calculate percent of CpGs with high quality log_lik 

summary(gm128_meth$log_lik_ratio)
m <- length(which(gm128_meth$log_lik_ratio>2.5))
u <- length(which(gm128_meth$log_lik_ratio<(-2.5)))
a <- length(gm128_meth$log_lik_ratio)
summary(gm128_meth$log_lik_ratio)

gm128_calls <- round((((m+u)/a)*100), digits = 2)
gm128_percentmeth <- round(((m/(m+u))*100), digits=2)
  

summary(v6mcalls$log_lik_ratio)
m <- length(which(v6mcalls$log_lik_ratio>2.5))
u <- length(which(v6mcalls$log_lik_ratio<(-2.5)))
a <- length(v6mcalls$log_lik_ratio)
summary(v6mcalls$log_lik_ratio)

chm13_calls <- round((((m+u)/a)*100), digits=2)
chm13_percentmeth <- round(((m/(m+u))*100), digits=2)

x <- c("gm12878", "chm13")
gm128_stat <- rbind(gm128_calls, gm128_percentmeth)
chm13_stat <- rbind(chm13_calls, chm13_percentmeth)
type <- c("Percent of CpGs called", "Percent methylated")
stat=as.data.frame(cbind(type, gm128_stat,chm13_stat))

names(stat) <- c("", "gm12878", "chm13")
tbl <- tableGrob(stat, rows=NULL)


plot <- ggplot(data = all_meth, aes(x=log_lik_ratio, fill = all_meth$sample))+geom_density(alpha=.5)+xlim(25,-25)+theme_bw()+scale_fill_manual("Assembly",values=c( "red", "green"))+geom_vline(xintercept= c(-2.5, 2.5), linetype="dashed", color="black")

ggarrange(plot, tbl, 
          ncol = 1, nrow = 2,
          heights = c(1, 0.5))

ggsave("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/plots/log_lik.pdf", height = 15, width = 10)
