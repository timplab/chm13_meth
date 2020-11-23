#!/usr/bin/Rscript

# script to generate censat methylation profiles 

# load libraries
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

# create parser object
parser <- ArgumentParser()
parser$add_argument("-r", "--region", type="character", 
                    help="Genomic region chrx:xx-xx")
parser$add_argument("-w", "--window", type="integer", default=200, 
                    help="Bin size [default %(default)s]")
parser$add_argument("-l", "--readlength", type="integer", default=50000, 
                    help="Minimum read length to include [default %(default)s]")

args <- parser$parse_args()
coords= str_split_fixed(args$region, ":",2)[1,2]
chrom=str_split_fixed(args$region, ":",2)[1,1]
rstart=as.numeric(str_split_fixed(coords, "-",2)[1,1])
rend=as.numeric(str_split_fixed(coords, "-",2)[1,2])
window=as.numeric(args$window)
minlen=as.numeric(args$readlength)
  
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
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
table(censat$name)
censat.gr <- GRanges(as.data.frame(censat))




repeatColors =c("HSAT2"="#C19935",
                "HSAT1"="#A2A638",
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
# read methylation GRanges data
chm13_meth <- readRDS(paste0(dat, "/methylation_calls/chm13_methylation_NanopolishFreq_50kb.rds"))
#chr3:94250000-94300000
rm <- read_tsv(paste0(dat, "/annotations/chm13.draft_v1.0.fasta_rm.bed"), col_names = F) %>%
  mutate(name = ifelse(X7 == "Satellite", X4, X7))


rm_sub <- rm %>%
  filter(X1 == chrom) %>%
  dplyr::filter(X2 > rstart) %>%
  dplyr::filter(X3 < rend)

bis <-  as.data.frame(chm13_meth) %>%
  filter(seqnames == chrom) %>%
  filter(start > rstart) %>%
  filter(end < rend) %>%
  mutate(smooth = rollmean(methylated_frequency, window, fill = NA), cov_smooth = rollmean(called_sites, window, fill = NA))

meth <- ggplot(bis, aes(x = (start/1e6), y= smooth))+geom_line(size =1) + labs(x="Genomic coordinates (Mb)", y="Methylation")+theme_classic(base_size = 25)+ylim(0,1)+xlim((rstart/1e6),(rend/1e6))

cov <- ggplot(bis, aes(x = (start/1e6), y= cov_smooth))+geom_line(size =1) + labs(x="Genomic coordinates (Mb)", y="Coverage")+theme_classic(base_size = 25)+xlim(rstart/1e6,rend/1e6)


rep_leg <- ggplot(data=rm_sub, mapping=aes(xmin=(X2/1e6),xmax=(X3/1e6),ymin=0,ymax=.1, fill = name))+
  geom_rect()+theme(legend.position="top") +labs(y="Axis") +labs(y="Axis")+ scale_fill_manual(values = repeatColors, drop = FALSE)+theme(legend.text=element_text(size=rel(1.5)))+theme(legend.title=element_blank()) +theme(axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))

top_row <- plot_grid(rep_leg,meth,cov, ncol = 1, align="v", rel_heights = c(1/4, 1/2, 1/2))


ggsave(
  paste0(chrom,  "-", rstart, "-", rend, "_methyl_profile.eps"),
  plot = top_row,
  width = 15,
  height = 10
)
