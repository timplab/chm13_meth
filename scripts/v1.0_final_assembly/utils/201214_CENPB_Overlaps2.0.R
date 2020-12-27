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



flankn <- 500


cenpb_meth <- readRDS(paste0(dat, "/CENPB_MethylatedSites.rds"))
cenpb_GpC <- readRDS(paste0(dat, "/CENPB_GC_readLevelOverlaps.rds"))

return.dat <- data.frame()
for (i in 1:length(unique(cenpb_meth$qname))){
  print(i)
  n=unique(cenpb_meth$qname)[i]
  
  regs.dat <- cenpb_meth %>%
    filter(qname == n) %>%
    mutate(seqnames = "chrX") %>%
    mutate(start = start -flankn, end = end +flankn) %>%
    GRanges()
  
  meth.dat <- cenpb_GpC %>%
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

saveRDS(return.dat, file = paste0(dat, "/CENPB_GC_readLevel_all.rds"))

