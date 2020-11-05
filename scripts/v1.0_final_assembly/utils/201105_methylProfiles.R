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
                "Satellite"="#4EB8DF",
             #   "Satellite"="#53B0E3",
                "SATR1"="#5AA5DA",
                "scRNA"="#6B9AD2",
                "Simple_repeat"="#8992C8",
                "SINE"="#9A8AC1",
                "snRNA"="#A885BC",
                "srpRNA"="#B67EB6",
                "SST1"="#C378B2",
                "Satellite"="#D173AF",
                "tRNA"="#ED72A5",
                "Unknown"="#EF768C", 
                "(GAATG)n"="#ff4000",
                "D20S16" = "#ffbf00", 
                "SATR2"= "#0080ff" )

defaultColor = "#000080"

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
                "ACRO"="#4EB8DF",
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
                "TAR"= "#0080ff" )


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
  mutate(name = ifelse(grepl("GSATII", X4), "GSATII", name)) %>%
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


for (i in 1:length(cen$chr)){
  chr=cen$chr[i]
  start=cen$start[i]
  end=cen$end[i]

bis <- as.data.frame(freq.matched) %>% 
  filter(seqnames == chr) %>%
  arrange(start) %>%
 # summarise(meth_freq = mean(mcall), cov = n()) %>%
  mutate(smooth = rollmean(meth, 500, fill = NA)) %>%
  mutate(cov_smooth = rollmean(cov, 10000, fill = NA))

rm_sub <- rm %>%
  filter(X1 == chr) %>%
  filter(X2 > start) %>%
  filter(X3 < end)

meth <- ggplot(bis, aes(x = (start/1e6), y= smooth))+geom_line(size =1) + labs(x="Genomic coordinates (Mb)", y="Methylation")+theme_classic(base_size = 25)+ylim(0,1)+xlim((start/1e6),(end/1e6))

cov <- ggplot(bis, aes(x = (start/1e6), y= cov_smooth))+geom_line(size =1) + labs(x="Genomic coordinates (Mb)", y="Coverage")+theme_classic(base_size = 25)+xlim(start/1e6,end/1e6)


rep_leg <- ggplot(data=rm_sub, mapping=aes(xmin=(X2/1e6),xmax=(X3/1e6),ymin=0,ymax=.1, fill = X7))+
  geom_rect()+theme(legend.position="top") +labs(y="Axis")+xlim(start/1e6,end/1e6) +labs(y="Axis")+ scale_fill_manual(values = repeatColors, drop = FALSE)+theme(legend.text=element_text(size=rel(1.5)))+theme(legend.title=element_blank()) +theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())


censat_reg <- censat %>%
  filter(chrom == chr) 

rep_leg2 <- ggplot(data=censat_reg, mapping=aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1, fill = name))+
  geom_rect()+theme(legend.position="top") +labs(y="Axis")+xlim(start/1e6,end/1e6) +labs(y="Axis")+ scale_fill_manual(values = censatColors, drop = FALSE)+theme(legend.text=element_text(size=rel(1.5)))+theme(legend.title=element_blank()) +theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())# +  theme_void()


top_row <- plot_grid(rep_leg,rep_leg2,meth,cov, ncol = 1, align="v", rel_heights = c(1/4,1/4, 1/2, 1/2))
top_row


ggsave(
  paste0(figs, "/", chr, "_HOR_methyl_pattern.eps"),
  plot = top_row,
  scale = 1,
  width = 10,
  height = 12,
  dpi = 50
)

}
