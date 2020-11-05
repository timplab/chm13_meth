#!/usr/bin/Rscript

# script to generate single read plots 

# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(source("/home/isac/Code/ilee/plot/ilee_plot_utils.R"))
suppressPackageStartupMessages(library("ggsci"))
suppressPackageStartupMessages(source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R"))

# load data
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


# create parser object
parser <- ArgumentParser()
parser$add_argument("-r", "--region", type="character", 
                    help="Genomic region chrx:xx-xx")
parser$add_argument("-l", "--readlength", type="integer", default=50000, 
                    help="Minimum read length to include [default %(default)s]")
parser$add_argument("-m", "--maxgap", type="integer", default=20, 
                    help="Maxgap sisze for merging CpGs when plotting [default %(default)s]")

args <- parser$parse_args()
coords= str_split_fixed(args$region, ":",2)[1,2]
chrom=str_split_fixed(args$region, ":",2)[1,1]
rstart=str_split_fixed(coords, "-",2)[1,1]
rend=str_split_fixed(coords, "-",2)[1,2]
window=args$window
minlen=args$readlength
gap=args$maxgap

reads <- tabix_mbed(paste0(dat, "/censat/", chrom, "_cen.mbed"),extcol = "motif",by = "read") 


size_sel <- reads %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= minlen) %>%
  mutate(motif = "CG")

cgcalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)

reg = cgcalls %>%
  dplyr::filter(start >= rstart) %>%
  dplyr::filter(end <= rend)

runs <- getRuns(reg, maxGap = gap)
cpg_runs.ordered <- order_reads(runs)

cpg_runs_reg <- cpg_runs.ordered$x %>%
  mutate(m = ifelse(values == 1, "Methylated","Unmethylated")) %>%
  mutate(mod = "CpG")

pal <- pal_npg("nrc")(10)
meth_pal <- c(pal[8],pal[4]) 
unsmooth <- ggplot(cpg_runs_reg,aes(xmin = start/1e6, xmax = end/1e6, ymin = ymin, ymax = ymax)) + geom_rect(data = cpg_runs.ordered$bounds, fill = "grey80") + 
  geom_rect(aes(fill = m))  +
  #    geom_vline(xintercept = 127638255, linetype == "dashed") +
  scale_fill_manual(name = "State", values = meth_pal) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.spacing = unit(2, "lines"), legend.position = "bottom")


ggsave(
  paste0(chrom,  "-", rstart, "-", rend, "_singleRead.pdf"),
  plot = unsmooth,
  width = 8,
  height = 8
)