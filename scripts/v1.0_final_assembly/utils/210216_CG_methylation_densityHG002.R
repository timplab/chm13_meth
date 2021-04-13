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

censat.gr <- read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), col_names = c("chr", "start", "end", "name")) %>%
  GRanges()

#write.table(censat, paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), sep="\t", col.names = F, row.names = F, quote=F)
# set repeat colors
censatColors =c("(CATTC)n" = "#E87C71",
                "(GAATC)n"="#E28455",
                "HOR"="#D78C32",
                "BSAT"="#E370AB",
                "CER" = "#CE9334",
                "HSAT2"="#C19935",
                "HSAT1"="#A2A638",
                "HSAT3"="#8CAC3E",
                "Low_complexity"="#75B042",
                "LSAU"="#54B346",
                "LTR"="#51B756",
                "MST"="#53BB73",
                "GSAT"="#55BE8D",
                "RNA"="#54C0A5",
                "rRNA"="#52BEBB",
                "SAR"="#51BDCE",
                "ACRO1"="#9400D3",
                "HSAT4"="#53B0E3",
                "SATR"="#5AA5DA",
                "CT"="#6B9AD2",
                "Simple_repeat"="#8992C8",
                "SINE"="#9A8AC1",
                "MON"="#A885BC",
                "SST"="#C378B2",
                "HSAT5"="#ED72A5",
                "HSAT6"="#EF768C", 
                "gap-rDNA"="#ff4000",
                "TE" = "#ffbf00", 
                "TAR"= "#0080ff",
                "ACRO"="#9400D3", 
                "DHOR" = "gray")

# load CG and methylation GRanges data 
chm13_CpG <- readRDS(paste0(dat, "/reference/chm13_sliding_200_CpG.rds"))
#chm13_meth <- readRDS(paste0(dat, "/methylation_calls/chm13_methylation_NanopolishFreq_50kb.rds"))

chm13_meth <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/whole_genome/chm13_hg002_reference_pooled/HG002_nanonome_CpGmethylationFrequency_20kb.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  GRanges()


# find overlaps between censat regions and methylation calls
censat.gr <- GRanges(censat)
keepi <- findOverlaps(chm13_meth,censat.gr)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

# split into groups for plotting

SAT = c("GSAT", "CER", "SST", "BSAT","HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6","HOR", "MON", "CT", "DHOR")


censat_meth <- as.data.frame(freq.matched) %>%
 filter(name %in% SAT) %>%
 filter(seqnames != "chrX")

# violing plot of methylation frequency

stats <- censat_meth %>%
  group_by(name) %>%
  summarize(med = median(methylated_frequency))

violin <- ggplot(data = censat_meth, aes(x = factor(name), y = methylated_frequency, fill = name))+geom_violin(width=1,adjust=2.5)+theme_classic(base_size = 12)+ scale_fill_manual(values = censatColors, drop = FALSE)+labs(x = "Methylation", y = "Repeat")+geom_boxplot(width=0.1,outlier.shape = NA)

#ggsave(
#  paste0(figs, "/", "AllSatv2.021621_methylation_freqHG002.pdf"),
#  plot = violin,
#  scale = 1,
#  width = 8,
#  height = 5,
#)

censat = read_tsv(paste0(dat, "/HG002/annotations/t2t_cenAnnotation.hg002_X.v1.bed"), col_names = c("X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8")) %>%
  mutate(name = ifelse(grepl("HSat1", X4), "HSAT1", X4)) %>%
  mutate(name = ifelse(grepl("HSat2", X4), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("HSat3", X4), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("HSat4", X4), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("HSat5", X4), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("HSat6", X4), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", X4), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", X4), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("hor", X4), "HOR", name)) %>%
  mutate(name = ifelse(grepl("dhor", X4), "DHOR", name)) %>%
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
    GRanges()



# find overlaps between censat regions and methylation calls
censat.gr <- GRanges(censat)
keepi <- findOverlaps(chm13_meth,censat.gr)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(censat.gr[subjectHits(keepi)]))

# split into groups for plotting

SAT = c("GSAT", "CER", "SST", "BSAT","HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6","HOR", "MON", "CT", "DHOR")


censat_meth <- as.data.frame(freq.matched) %>%
  filter(name %in% SAT) 

# violing plot of methylation frequency

stats <- censat_meth %>%
  group_by(name) %>%
  summarize(med = median(methylated_frequency))


violin <- ggplot(data = censat_meth, aes(x = factor(name), y = methylated_frequency, fill = name))+geom_violin(width=1,adjust=2.5)+theme_classic(base_size = 12)+ scale_fill_manual(values = censatColors, drop = FALSE)+labs(x = "Methylation", y = "Repeat")+geom_boxplot(width=0.1,outlier.shape = NA)

ggsave(
  paste0(figs, "/", "AllSatv1_chrX_methylation_freqHG002.pdf"),
  plot = violin,
  scale = 1,
  width = 8,
  height = 5,
)