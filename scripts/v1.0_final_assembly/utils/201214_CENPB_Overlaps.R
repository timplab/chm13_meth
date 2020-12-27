#!/usr/bin/Rscript

library(knitr)
library(tidyverse)
library(bsseq)
library(Biostrings)
library(ggplot2)
library(png)
library(cowplot)
options(scipen=999)
library(zoo)
options(knitr.duplicate.label = 'allow')
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")
library(mclust)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002"
hg002_gc <- read_csv(paste0(dat, "/nanonome/methylation_calls/pooled/HG002_nanonome_chrX_GpCmethylation_ByCall.csv"))
hg002_meth <- read_csv(paste0(dat, "/nanonome/methylation_calls/pooled/HG002_nanonome_chrX_CpGmethylation_ByCall.csv")) %>%
  GRanges() 


cenpb <- read_delim(paste0(dat, "/reference/chrx_fuzznuc.bed"), delim = " ") %>%
  as.data.frame() %>%
  ungroup() %>%
  mutate("start" = as.numeric(`  Start`), "end" = as.numeric(`    End`)) %>%
  select(-c (`  Start`, `    End`)) %>%
  mutate(center = (end + start)/2, chr = "chrX") %>%
  mutate(id = row_number()) %>%
  mutate("rstart" = start, "rend" = end) %>%
  GRanges()


keepi <- findOverlaps(hg002_meth,cenpb)

freq.matched <- hg002_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(cenpb[subjectHits(keepi)]))

cenpb_meth <- as.data.frame(freq.matched) %>%
  group_by(id, qname, center) %>%
  summarise(meth = sum(mcall), start = rstart, end = rend) %>%
  distinct()

hg002.gr <- hg002_gc %>%
  GRanges()

flankn <- 500

return.dat <- data.frame()
for (i in 1:length(unique(cenpb_meth$qname))){
  
  n=unique(cenpb_meth$qname)[i]
  
  regs.dat <- cenpb_meth %>%
    filter(qname == n) %>%
    mutate(seqnames = "chrX") %>%
    mutate(start = start -flankn, end = end +flankn) %>%
    GRanges()
  
  meth.dat <- hg002_gc %>%
    filter(qname == n) %>%
    GRanges()
  
  keepi <- findOverlaps(meth.dat,regs.dat)
  
  freq.matched <- meth.dat[queryHits(keepi)]
  
  mcols(freq.matched) <- cbind.data.frame(
    mcols(freq.matched),
    mcols(regs.dat[subjectHits(keepi)]))
  
  df <- as.data.frame(freq.matched)
  return.dat <- rbind(df, return.dat)
  
}

saveRDS(return.dat, file = paste0(dat, "/CENPB_GC_readLevel.rds"))
