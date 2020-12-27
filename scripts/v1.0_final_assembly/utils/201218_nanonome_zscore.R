library(tidyverse)
library(Biostrings)
options(scipen=999)
library(zoo)
library(GRanges)
library(Repitools)
library(ggpmisc)
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/nanonome/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/HG002/nanonome/methylation_calls"


gccalls <- read_csv(paste0(dat, "/pooled/HG002_nanonome_chrX_GpCmethylation_ByCall.csv"))

calls.gr <-  gccalls %>%
  ungroup() %>%
  group_by(start) %>%
  summarise(num_meth = sum(mcall == 1), num_unmeth = sum(mcall == 0)) %>%
  mutate(chr = "chrX", start = start, end = start) %>%
  GRanges()

library(BSgenome.HG002.chrX)
chrx.gr <- GRanges(seqinfo(BSgenome.HG002.chrX))


blocks <- genomeBlocks(BSgenome.HG002.chrX, chrs = seqnames(BSgenome.HG002.chrX), width = 15000)
binnedSum <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewSums(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

score1 <- coverage(calls.gr, weight="num_meth")
score2 <- coverage(calls.gr, weight="num_unmeth")

binned_meth <- binnedSum(blocks, numvar = score1, "num_meth") %>%
  as.data.frame()
binned_unmeth <-binnedSum(blocks, numvar = score2, "num_unmeth")%>%
  as.data.frame()

meth_bins <- merge(binned_meth, binned_unmeth, by = c("start", "end", "seqnames", "width", "strand")) %>%
  group_by(start, end) %>%
  mutate(sites = num_meth+num_unmeth) %>%
  mutate(freq = num_meth/sites) %>%
  na.omit()

freqmean=mean(meth_bins$freq)
freqsd=sd(meth_bins$freq)

meth_bins <- meth_bins %>%
  mutate(z_score = (freq - freqmean) / freqsd) %>%
  mutate(color = case_when(z_score > 0 ~ "Accessible", 
                           TRUE ~ "Inaccessible")) 
saveRDS(meth_bins, paste0(dat, "/ChrX_accessibilityZscore.rds"))


high <- meth_bins %>%
  filter(z_score > 3 | z_score < -3)

write.table(high, paste0(dat, "/pooled/ChrX_accessibilityZscore.tsv"), quote = F, row.names = F, sep = "\t")
