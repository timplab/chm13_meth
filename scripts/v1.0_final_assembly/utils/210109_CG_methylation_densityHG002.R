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

# load and parse annotation

censat.gr <- read_tsv(paste0(dat, "/annotations/hsat/chm13.hsat3.final.strand.subfams.bed"), col_names = c("chr", "start", "end", "name")) %>%
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


censat_meth <- as.data.frame(freq.matched) 

# violing plot of methylation frequency

stats <- censat_meth %>%
  group_by(name) %>%
  summarize(med = median(methylated_frequency))

violin <- ggplot(data = censat_meth, aes(x = factor(name), y = methylated_frequency, fill = name))+geom_violin()+theme_classic(base_size = 12)+labs(x = "Subtype", y = "Methylation")+geom_boxplot(width=0.1,outlier.shape = NA)+geom_hline(yintercept = .3680, linetype= "dashed") +theme(legend.position = "bottom")

ggsave(
  paste0(figs, "/", "210109_HSat3_methylation_freq.pdf"),
  plot = violin,
  scale = 1,
  width = 10,
  height = 5,
)



keepi <- findOverlaps(chm13_CpG,censat.gr)
freq.matched <- chm13_CpG[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))



cpg <- as.data.frame(freq.matched) 
  
stats.cpg <- cpg %>%
  group_by(name) %>%
  summarize(med = median(CpG))

  # violing plot of methylation frequency

violin <- ggplot(cpg, aes(x = factor(name), y = CpG, fill = name))+geom_violin(adjust=8)+theme_classic(base_size = 20)+labs(x = "Repeat", y = "CpG Density")+geom_boxplot(outlier.shape = NA, width=.1)+coord_cartesian(y=c(0,30))+geom_hline(yintercept = 3, linetype= "dashed") +theme(legend.position = "none")
  
  ggsave(
    paste0(figs, "/210109_Hsat3_CpG_boxplot.pdf"),
    plot = violin,
    scale = 1,
    width = 10,
    height = 5,
  )
 