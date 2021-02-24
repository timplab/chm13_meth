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
library(GenomicRanges)
library(BSgenome.t2t.v1.0.release)
library(ggridges)
options(knitr.duplicate.label = 'allow')
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")

# set data output path

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/censat/figures/violin"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

# load and parse censat files

censat = read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.cenSat_annotation.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("HSat1", X4), "HSAT1", X4)) %>%
  mutate(name = ifelse(grepl("HSat2", X4), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("HSat3", X4), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("HSat4", X4), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("HSat5", X4), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("HSat6", X4), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", X4), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", X4), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("hor", X4), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", X4), "MON", name)) %>%
  mutate(name = ifelse(grepl("Alu", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("SATR", X4), "SATR", name)) %>%
  mutate(name = ifelse(grepl("ACRO", X4), "ACRO", name)) %>%
  mutate(name = ifelse(grepl("GSATII", X4), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("TAR", X4), "TAR", name)) %>%
  mutate(name = ifelse(grepl("TE", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("MER", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("MST", X4), "MST", name)) %>%
  mutate(name = ifelse(grepl("CER", X4), "CER", name)) %>%
  mutate(name = ifelse(grepl("L1", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("SST", X4), "SST", name)) %>%
  mutate(name = ifelse(grepl("LSAU", X4), "LSAU", name)) %>%
  mutate(name = ifelse(grepl("GSAT", X4), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("MSAT", X4), "MSAT", name)) %>%
  mutate(name = ifelse(grepl("novel", X4), "novel", name)) %>%
  mutate(name = ifelse(grepl("HERV", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("LT", X4), "TE", name)) %>%
  dplyr::select(c(X1, X2, X3, name)) %>%
  dplyr::rename("chrom" =1, "start" = 2 ,"end" =3) %>%
  group_by(name) %>%
  filter(n() >= 3) %>%
  ungroup()
table(censat$name)
censat.gr <- GRanges(as.data.frame(censat))

# set repeat colors
repeatColors =c(#"HSATII"="#C19935",
                "HSATI"="#A2A638",
                "HSAT3"="#8CAC3E",
                "(CATTC)n" = "#E87C71",
                "HSAT4"="#53B0E3",
                "HSAT5"="#ED72A5",
                "HSAT6"="#EF768C", 
                "(GAATC)n"="#E28455",
                "ALR/Alpha"="#D78C32",
                "6kbHsap" = "#D8BFD8",
                "BSR/Beta"="#E370AB",
                "CER" = "#CE9334",
                "DNA"="#C19935",
                "DNA?"="#C19935",
                "GSAT"="#4169E1",
                "LINE"="#FFA500",
                "Low_complexity"="#75B042",
                "LSAU"="#54B346",
                "LTR"="#51B756",
                "RC"="#53BB73",
                "Retroposon"="#55BE8D",
                "RNA"="#54C0A5",
                "rRNA"="#52BEBB",
                "SAR"="#51BDCE",
                "SATR2"="#4EB8DF",
                "GSATII"="#6B8E23",
                "SATR1"="#5AA5DA",
                "scRNA"="#6B9AD2",
                "Simple_repeat"="#8992C8",
                "SINE"="#9A8AC1",
                "snRNA"="#A885BC",
                "srpRNA"="#B67EB6",
                "SST1"="#C378B2",
                "HSATII"="#D173AF",
                "tRNA"="#006400",
                "ACRO1"="#9400D3",
                "Unknown"="#BA55D3", 
                "(GAATG)n"="#ff4000",
                "D20S16" = "#ffbf00", 
                "SATR2"= "#0080ff", 
                "TAR1" = "#000080", 
                "SUBTEL2_sat"= "#FFB6C1", 
                "GSATX" = "#D2691E", 
                "MSR1" = "#708090")
# load CG and methylation GRanges data 
chm13_CpG <- readRDS(paste0(dat, "/reference/chm13_sliding_200_CpG.rds"))
chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  GRanges()

# whole genome methylation
all_meth <- as.data.frame(chm13_meth) %>%
  mutate(type = ifelse(seqnames == "chrX", "chrX", "autosome"))



violin <- ggplot(data = all_meth, aes(y = factor(type), x = methylated_frequency, fill = type))+geom_density_ridges(rel_min_height = 0.01,scale = 0.9,quantile_lines = TRUE, quantiles = 2, bandwidth = .015)+theme_classic(base_size = 20)+labs(x = "Methylation")

ggsave(
  paste0(figs, "/WG_methylation_freq.pdf"),
  plot = violin,
  scale = 1,
  width = 8,
  height = 5,
)



#repeatmasker
rm.gr <- read_tsv(paste0(dat, "/annotations/chm13.draft_v1.0.fasta_rm.bed"), col_names = F) %>%
  #mutate(name = ifelse(X7 == "Satellite", X4, X7)) %>%
  rename(X1 = "chr", X2 = "start", X3 = "end") %>%
  GRanges()

keepi <- findOverlaps(chm13_meth,rm.gr)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(rm.gr[subjectHits(keepi)]))

list = c("Simple_repeat", "SINE", "LINE", "LTR", "DNA", "Satellite")
# violing plot of methylation frequency
rep_meth <- as.data.frame(freq.matched) %>%
  filter(X7 %in% list)

  violin <- ggplot(data = rep_meth, aes(y = factor(X7), x = methylated_frequency, fill = X7))+geom_density_ridges(rel_min_height = 0.01,scale = 0.9,quantile_lines = TRUE, quantiles = 2, bandwidth = .015)+theme_classic(base_size = 20)+ scale_fill_manual(values = repeatColors, drop = FALSE)+labs(x = "Methylation", y = "Repeat")

ggsave(
  paste0(figs, "/RM_methylation_freq.pdf"),
  plot = violin,
  scale = 1,
  width = 8,
  height = 5,
)



keepi <- findOverlaps(chm13_CpG,rm.gr)
freq.matched <- chm13_CpG[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(rm.gr[subjectHits(keepi)]))



list = c("Simple_repeat", "SINE", "LINE", "LTR", "DNA", "Satellite")
# violing plot of methylation frequency
rep_cg <- as.data.frame(freq.matched) %>%
  filter(X7 %in% list)

violin <- ggplot(data = rep_cg, aes(y = factor(X7), x = CpG, fill = X7))+geom_density_ridges(rel_min_height = 0.01,scale = 0.9,quantile_lines = TRUE, quantiles = 2, bandwidth = 1)+theme_classic(base_size = 20)+ scale_fill_manual(values = repeatColors, drop = FALSE)+labs(x = "CpG Density", y = "Repeat")+xlim(-3,20)

ggsave(
  paste0(figs, "/RM_CpG_freq.pdf"),
  plot = violin,
  scale = 1,
  width = 8,
  height = 5,
)


rm.gr <- read_tsv(paste0(dat, "/annotations/chm13.draft_v1.0.fasta_rm.bed"), col_names = F) %>%
  mutate(name = ifelse(X7 == "Satellite", X4, X7)) %>%
  rename(X1 = "chr", X2 = "start", X3 = "end") %>%
  GRanges()


keepi <- findOverlaps(chm13_meth,rm.gr)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(rm.gr[subjectHits(keepi)]))

rep_meth <- as.data.frame(freq.matched) %>%
  filter(X7 == "Satellite")
violin <- ggplot(data = rep_meth, aes(y = factor(name), x = methylated_frequency, fill = name))+geom_density_ridges(rel_min_height = 0.01,scale = 0.9,quantile_lines = TRUE, quantiles = 2, bandwidth = .015)+theme_classic(base_size = 20)+ scale_fill_manual(values = repeatColors, drop = FALSE)+labs(x = "Methylation", y = "Repeat")

ggsave(
  paste0(figs, "/all_sat_methylation_freq.pdf"),
  plot = violin,
  scale = 1,
  width = 10,
  height = 10,
)


keepi <- findOverlaps(chm13_CpG,rm.gr)
freq.matched <- chm13_CpG[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(rm.gr[subjectHits(keepi)]))


rep_meth <- as.data.frame(freq.matched) %>%
  filter(X7 == "Satellite")
violin <- ggplot(data = rep_meth, aes(y = factor(name), x = CpG, fill = name))+geom_density_ridges(rel_min_height = 0.01,scale = 0.9,quantile_lines = TRUE, quantiles = 2, bandwidth = .5)+theme_classic(base_size = 20)+ scale_fill_manual(values = repeatColors, drop = FALSE)+labs(x = "CpG frequency", y = "Repeat")

ggsave(
  paste0(figs, "/all_sat_CpG_freq.pdf"),
  plot = violin,
  scale = 1,
  width = 10,
  height = 10,
)



