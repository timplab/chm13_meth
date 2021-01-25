library(BSgenome.t2t.v1.0.release)
library(tidyverse)
library(GenomicRanges)
library(zoo)
options(scipen=999)
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  GRanges()


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

binnedMean <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewMeans(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}




blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = seqnames(BSgenome.t2t.v1.0.release), width = 100000)
score1 <- coverage(chm13_meth, weight="called_sites_methylated")
score2 <- coverage(chm13_meth, weight="called_sites_unmethylated")
score3 <- coverage(chm13_meth, weight="num_motifs_in_group")


binned_meth <- binnedSum(blocks, numvar = score1, "called_sites_methylated") %>%
  as.data.frame()
binned_unmeth <-binnedSum(blocks, numvar = score2, "called_sites_unmethylated")%>%
  as.data.frame()
binned_cov <- binnedMean(blocks, numvar = score3, "num_motifs_in_group") %>%
  as.data.frame()



meth_bins <- merge(binned_meth, binned_unmeth, by = c("start", "end", "seqnames", "width", "strand")) %>%
  merge(binned_cov, by = c("start", "end", "seqnames", "width", "strand")) %>%
  group_by(seqnames,start, end) %>%
  mutate(sites = called_sites_methylated+called_sites_unmethylated) %>%
  mutate(freq =called_sites_methylated/sites) %>%
  ungroup() %>%
  group_by(seqnames) %>%
  arrange(start) %>%
  mutate(smooth = rollmean(freq, 5, fill = NA), site_smooth = rollmean(num_motifs_in_group, 3, fill = NA)) %>%
  ungroup() %>%
  dplyr::select(c(seqnames, start, end,freq, smooth)) %>%
  arrange(seqnames, start) %>%
  mutate(across(where(is.numeric), ~ round(., digits = 2)))

chr1 <- meth_bins %>%
  filter(seqnames=="chr1")

ggplot(chr1, aes(x = (start/1e6), y= smooth))+geom_line(size =1) + labs(x="Genomic coordinates (Mb)", y="Methylation")+theme_classic(base_size = 25)+ylim(0,1)

write.table(meth_bins, file = paste0(dat, "/methylation_calls/210124_CHM13v1.0_MethylationFreq_100kbSmooth5.bed"), quote = F, sep = "\t",dec = ".", row.names = F,
            col.names = F)
