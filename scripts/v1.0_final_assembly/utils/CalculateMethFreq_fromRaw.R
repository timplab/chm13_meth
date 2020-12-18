#!/usr/bin/Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(zoo))
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

reads <- tabix_mbed(paste0(dat, "/censat/chrX_cen.mbed"),extcol = "motif",by = "call") 

size_sel <- reads %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= 50000) %>%
  mutate(motif = "CG")

cgcalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)


bis <- cgcalls %>% 
  filter(start > rstart) %>%
  filter(end < rend) %>%
  group_by(chrom, start) %>%
  summarise(meth_freq = mean(mcall), cov = n())

saveRDS(bis, file =paste0(dat, "/methylation_calls/chm13_methylation_50kb_fromRaw.rds"))

smooth <- bis  %>%
  ungroup() %>%
  group_by(chrom) %>%
  summarise(meth_smooth = rollmean(meth_freq, 200, fill = NA), cov_smooth = rollmean(cov, 1000, fill = NA), chrom= chrom, start = start, end = start) %>%
  na.omit()

saveRDS(smooth, file =paste0(dat, "/methylation_calls/chm13_methylation_50kb_200bpRollMean_fromRaw.rds"))



