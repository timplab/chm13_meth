#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(bsseq)
library(Biostrings)
library(ggplot2)
library(png)
library(cowplot)
options(scipen=999)
library(zoo)
library(BSgenome.t2t.v1.0.release)
options(knitr.duplicate.label = 'allow')
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")

# set data output path

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/censat/figures/violin"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

# load and parse censat files

#censat = read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.cenSat_annotation.bed"), col_names = F) #%>%
#  mutate(name = ifelse(grepl("HSat1", X4), "HSAT1", X4)) %>%
#  mutate(name = ifelse(grepl("HSat2", X4), "HSAT2", name)) %>%
#  mutate(name = ifelse(grepl("HSat3", X4), "HSAT3", name)) %>%
#  mutate(name = ifelse(grepl("HSat4", X4), "HSAT4", name)) %>%
#  mutate(name = ifelse(grepl("HSat5", X4), "HSAT5", name)) %>%
#  mutate(name = ifelse(grepl("HSat6", X4), "HSAT6", name)) %>%
#  mutate(name = ifelse(grepl("ct", X4), "CT", name)) %>%
#  mutate(name = ifelse(grepl("bsat", X4), "BSAT", name)) %>%
#  mutate(name = ifelse(grepl("hor", X4), "HOR", name)) %>%
#  mutate(name = ifelse(grepl("mon", X4), "MON", name)) %>%
#  mutate(name = ifelse(grepl("Alu", X4), "TE", name)) %>%
#  mutate(name = ifelse(grepl("SATR", X4), "SATR", name)) %>%
#  mutate(name = ifelse(grepl("ACRO", X4), "ACRO", name)) %>%
#  mutate(name = ifelse(grepl("GSATII", X4), "GSAT", name)) %>%
#  mutate(name = ifelse(grepl("TAR", X4), "TAR", name)) %>%
#  mutate(name = ifelse(grepl("TE", X4), "TE", name)) %>%
#  mutate(name = ifelse(grepl("MER", X4), "TE", name)) %>%
#  mutate(name = ifelse(grepl("MST", X4), "MST", name)) %>%
#  mutate(name = ifelse(grepl("CER", X4), "CER", name)) %>%
#  mutate(name = ifelse(grepl("L1", X4), "TE", name)) %>%
#  mutate(name = ifelse(grepl("SST", X4), "SST", name)) %>%
#  mutate(name = ifelse(grepl("LSAU", X4), "LSAU", name)) %>%
#  mutate(name = ifelse(grepl("GSAT", X4), "GSAT", name)) %>%
#  mutate(name = ifelse(grepl("MSAT", X4), "MSAT", name)) %>%
#  mutate(name = ifelse(grepl("novel", X4), "novel", name)) %>%
#  mutate(name = ifelse(grepl("HERV", X4), "TE", name)) %>%
#  mutate(name = ifelse(grepl("LT", X4), "TE", name)) %>%
#  dplyr::select(c(X1, X2, X3, name)) %>%
#  dplyr::rename("chrom" =1, "start" = 2 ,"end" =3) %>%
#  group_by(name) %>%
#  #filter(n() >= 3) %>%
#  ungroup()
#table(censat$name)
chroms=c("chr1", "chr16")
censat.gr <- read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.cenSat_annotationFormatted.bed"), col_names = c("chr", "start", "end", "name")) %>%
  filter(name=="HSAT2") %>%
  filter(chr %in% chroms) %>%
  mutate(len=end-start) %>%
  filter(len > 1e6) %>%
  GRanges()

# load CG and methylation GRanges data 
chm13_CpG <- readRDS(paste0(dat, "/reference/chm13_sliding_200_CpG.rds"))
#chm13_meth <- readRDS(paste0(dat, "/methylation_calls/chm13_methylation_NanopolishFreq_50kb.rds"))

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  GRanges()


# find overlaps between censat regions and methylation calls

keepi <- findOverlaps(chm13_meth,censat.gr)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

stats <- as.data.frame(freq.matched) %>%
  group_by(name,seqnames) %>%
  summarize(med = median(methylated_frequency))


keepi <- findOverlaps(chm13_CpG,censat.gr)
freq.matched <- chm13_CpG[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

stats.cg <- as.data.frame(freq.matched) %>%
  group_by(name,seqnames) %>%
  summarize(med = median(CpG))

