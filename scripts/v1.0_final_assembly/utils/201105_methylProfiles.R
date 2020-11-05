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

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/censat/figures/all_chr_profiles"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

# read methylation GRanges data
chm13_meth <- readRDS(paste0(dat, "/methylation_calls/chm13_methylation_50kb.rds"))

# centromere boundaries
cen <- read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.cenSat_regions.bed"), col_names = c("chr", "start", "end","name"))
cen.gr <- GRanges(cen)

#repeatmasker
rm <- read_tsv(paste0(dat, "/annotations/chm13.draft_v1.0.fasta_rm.bed"), col_names = F)

keepi <- findOverlaps(chm13_meth,cen.gr)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(cen.gr[subjectHits(keepi)]))


# use repeatmasker tracks for annotation, specify rm colors


repeatColors =c("(CATTC)n" = "#E87C71",
                "(GAATC)n"="#E28455",
                "ALR/Alpha"="#D78C32",
                "BSR/Beta"="#E370AB",
                "CER" = "#CE9334",
                "DNA"="#C19935",
                "DNA?"="#C19935",
                "GSAT"="#B3A033",
                "HSAT"="#A2A638",
                "LINE"="#8CAC3E",
                "Low_complexity"="#75B042",
                "LSAU"="#54B346",
                "LTR"="#51B756",
                "RC"="#53BB73",
                "Retroposon"="#55BE8D",
                "RNA"="#54C0A5",
                "rRNA"="#52BEBB",
                "SAR"="#51BDCE",
                "Satellite/acro"="#4EB8DF",
                "Satellite/telo"="#53B0E3",
                "SATR1"="#5AA5DA",
                "scRNA"="#6B9AD2",
                "Simple_repeat"="#8992C8",
                "SINE"="#9A8AC1",
                "snRNA"="#A885BC",
                "srpRNA"="#B67EB6",
                "SST1"="#C378B2",
                "Satellite/subtelo"="#D173AF",
                "tRNA"="#ED72A5",
                "Unknown"="#EF768C", 
                "(GAATG)n"="#ff4000",
                "D20S16" = "#ffbf00", 
                "SATR2"= "#0080ff" )

defaultColor = "#000080"

i=1
for (i in 1:length(cen$chr)){
  chr=cen$chr[i]
  start=cen$start[i]
  end=cen$end[i]

bis <- as.data.frame(freq.matched) %>% 
  filter(seqnames == chr) %>%
  arrange(start) %>%
 # summarise(meth_freq = mean(mcall), cov = n()) %>%
  mutate(smooth = rollmean(meth, 200, fill = NA)) %>%
  mutate(cov_smooth = rollmean(cov, 200, fill = NA))

rm_sub <- rm %>%
  filter(X1 == chr) %>%
  filter(X2 > start) %>%
  filter(X3 < end)

meth <- ggplot(bis, aes(x = (start/1e6), y= smooth))+geom_line(size =1) + labs(x="Genomic coordinates (Mb)", y="Methylation")+theme_classic(base_size = 25)+ylim(0,1)+xlim((start/1e6),(end/1e6))

cov <- ggplot(bis, aes(x = (start/1e6), y= cov_smooth))+geom_line(size =1) + labs(x="Genomic coordinates (Mb)", y="Coverage")+theme_classic(base_size = 25)+xlim(start/1e6,end/1e6)


rep_leg <- ggplot(data=rm_sub, mapping=aes(xmin=(X2/1e6),xmax=(X3/1e6),ymin=0,ymax=.1, fill = X7))+
  geom_rect()+theme(legend.position="top") +labs(y="Axis")+xlim(start/1e6,end/1e6) +labs(y="Axis")+ scale_fill_manual(values = repeatColors, drop = FALSE)+theme(legend.text=element_text(size=rel(1.5)))+theme(legend.title=element_blank()) +theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())

top_row <- plot_grid(rep_leg,meth,cov, ncol = 1, align="v", rel_heights = c(1/4, 1/2, 1/2))
top_row

ggsave(
  paste0(figs, "/", chr, "_HOR_methyl_pattern.pdf"),
  plot = top_row,
  scale = 1,
  width = 12,
  height = 15,
)

}
