#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)


meth <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/chrX.methylation_frequency.tsv")
meth$pos <- rep("chrx", nrow(meth))
meth[meth$start >= 57828561 & meth$end  <= 60664792 , ][, "pos"] <- "cenx"
meth[meth$start >=  113868842 & meth$end <= 114116851,][, "pos"] <- "dxz4"


ggplot(data = meth, aes(x = methylated_frequency, fill=meth$pos)) +geom_density(alpha=.5)+facet_wrap(~meth$pos, ncol=1)+theme_bw()+labs(x="methylated frequency (per base)", color="Position")+scale_fill_manual("Region",values=c( "red", "green", "blue"))

binned <- meth %>%
  mutate(bin = ntile(start, (length(meth$start)/50))) %>%
  group_by(bin, pos) %>%
  summarise(avgmeth = mean(methylated_frequency)) %>%
  arrange(desc(avgmeth))

ggplot((data=binned),aes(x=avgmeth, fill = pos))+geom_density(alpha=.5)+theme_bw()+labs(x="Average methylation frequency (bin size 50bps)")+facet_wrap(~pos, ncol=1)+scale_fill_manual("Region",values=c( "red", "green", "blue"))
  

high <- binned %>%
  filter(avgmeth >=.95) %>%
  arrange(desc(avgmeth))
highbins <- high$bin

low <- binned %>%
  filter(avgmeth <=.1) %>%
  filter(pos == "cenx") %>%
  arrange(avgmeth)
lowbins <- low$bin

highmeth <- meth %>%
  mutate(bin = ntile(start, (length(meth$start)/50))) %>%
  filter(bin%in%c(highbins)) %>%
  group_by(bin, pos) 

write.table(highmeth, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/high_meth_freq.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)

lowmeth <- meth %>%
  mutate(bin = ntile(start, (length(meth$start)/50))) %>%
  filter(bin%in%c(lowbins)) %>%
  group_by(bin, pos) %>%
  filter(pos == "cenx")

write.table(lowmeth, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/low_meth_freq.tsv", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)


