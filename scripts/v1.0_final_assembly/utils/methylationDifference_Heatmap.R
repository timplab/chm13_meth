#!/usr/bin/Rscript

# script to generate pairwise methylation heatmap 

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
chm13_meth <- readRDS(paste0(dat, "/methylation_calls/chm13_methylation_50kb.rds"))


# create parser object
parser <- ArgumentParser()
parser$add_argument("-r", "--region", type="character", 
                    help="Genomic region chrx:xx-xx")
parser$add_argument("-w", "--window", type="integer", default=1000, 
                    help="Bin size [default %(default)s]")
parser$add_argument("-l", "--readlength", type="integer", default=50000, 
                    help="Minimum read length to include [default %(default)s]")

args <- parser$parse_args()
coords= str_split_fixed(args$region, ":",2)[1,2]
chrom=str_split_fixed(args$region, ":",2)[1,1]
rstart=str_split_fixed(coords, "-",2)[1,1]
rend=str_split_fixed(coords, "-",2)[1,2]
window=args$window
minlen=as.numeric(args$readlength)

reads <- tabix_mbed(paste0(dat, "/censat/", chrom, "_cen.mbed"),extcol = "motif",by = "read") 

size_sel <- reads %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= minlen) %>%
  mutate(motif = "CG")

cgcalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)


chm13_meth <- cgcalls %>% 
  filter(start > rstart) %>%
  filter(end < rend) %>%
  group_by(chrom, start) %>%
  summarise(meth = mean(mcall)) %>%
  mutate(end = start) %>%
  relocate(end, .before = meth) %>%
  GRanges()


# generate windows, get calls that overlap windows
win.gr <- data.frame(seqnames = chrom, start = rstart, end= rend) %>%
  GRanges()

sliding_blocks <- tile(win.gr, width = window)

ran <- data.frame(sliding_blocks) %>%
  mutate(bin = row_number()) %>%
  mutate(pos_start = start) %>%
  GRanges()


keepi <- findOverlaps(chm13_meth,ran)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(ran[subjectHits(keepi)]))

meth <- as.data.frame(freq.matched) %>%
  mutate(pos = start - pos_start) %>%
  mutate(meth = as.numeric(meth)) %>%
  select(meth, bin, pos) %>%
  group_by(bin) %>%
  summarise(avgmeth = median(meth)) %>%
  na.omit()


# co-occurance matrix formation
N <- meth$avgmeth
i <- seq(meth$avgmeth)
j <- seq(meth$avgmeth)

mat.end <- matrix(nrow = length(N), ncol = length(N))

for (n in i){
  for (y in j) {
    meth1 = N[n]
    meth2 = N[y]
    mat.end[n,y]=(meth1-meth2)} 
}



# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


# for aesthetic purposes - get bin quantiles and corresponding coordinates for labelling the axes with the genomic coordinates instead of the bin number 
coord <- as.data.frame(freq.matched) %>%
  mutate(meth = as.numeric(meth)) %>%
  select(meth, bin, start) %>%
  group_by(bin) %>%
  mutate(avgmeth = median(meth)) %>%
  na.omit() %>%
  select(bin, start) %>%
  distinct() %>%
  group_by(bin) %>%
  summarise(start1 = min(start)) %>%
  rename(bin = "Var1") 

quant <- quantile(coord$Var1,c(0, .25, .5, .75, 1), FALSE, TRUE, 3)

labs <- coord %>%
  filter(Var1 %in% quant) %>%
  mutate(start1 = round(start1/1e6, digits = 2))


upper_tri <- get_upper_tri(mat.end) 

melted_cormat <- melt(upper_tri, na.rm = TRUE) 



# Heatmap
cols <- brewer.pal(5, "Spectral")
rev_cols <- rev(cols)
pal <- colorRampPalette(rev_cols)

p <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradientn(colours =pal(5), name="Methylation difference") + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  theme(axis.text.y = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_fixed()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())+
  scale_x_continuous(breaks  = quant, labels = labs$start1)+scale_y_continuous(breaks  = quant, labels = labs$start1)

# save heatmap

ggsave(
  paste0(chrom,  "-", rstart, "-", rend, "_heatmap.pdf"),
  plot = p,
  width = 8,
  height = 8
)

print("bye bitch. plot done.")