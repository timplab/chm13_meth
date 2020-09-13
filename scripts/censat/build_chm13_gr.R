library(tidyverse)
library(bsseq)
library(biovizBase)
library(Repitools)
library(Biostrings)
library(ggplot2)
library(BSgenome)
library(zoo)
library(GenomicRanges)
library(BSgenome.t2t.v1.1)

dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly"

###functions#####

CalculateCpG <- function (obj, ..., step, width, as.prob, as.array, fast.moving.side,with.labels) 
{
  require(BSgenome)
  seqs <- getSeq(obj, ..., as.character = FALSE)
  if (missing(step)) {
    res <- oligonucleotideFrequency(seqs, 2, as.prob = as.prob)[,7]
  }
  else {
    res <- lapply(seqs, function(x) oligonucleotideFrequency(seqs, width = 2, step = step, as.prob = as.prob, as.array = as.array,fast.moving.side=fast.moving.side,with.labels=with.labels)[,7])
    if (length(res) == 1) 
      res <- unlist(res)
  }
  res
}
########


# build t2t BSgenome package and load it as BSgenome.t2t.v1.1

# split genome into 200bp bins and calculate CpG density per bin, save as GRanges object for easy reloading, run once comment out

chm13.gr <- GRanges(seqinfo(BSgenome.t2t.v1.1))
blocks <- genomeBlocks(BSgenome.t2t.v1.1, chrs = seqnames(BSgenome.t2t.v1.1), width = 200)
chm13_CG <- CalculateCpG(BSgenome.t2t.v1.1, blocks, as.prob = F)

chm13 <- GRanges(blocks, CpG = chm13_CG)
saveRDS(chm13, file =paste0(dat, "/ref/chm13_step_200_CpG_num.rds"))
#chm13 <- readRDS(paste0(dat, "/ref/chm13_step_200_CpG_freq.rds"))

freq = read_tsv(paste0(dat, "/whole_genome/CpGmethylation50kb.freq"), col_names = F)
bis <- read.bismark(paste0(dat, "/whole_genome/CpGmethylation50kb.freq"))
cov <- getCoverage(bis)
getmeth <- getMeth(bis, type = "raw")
getmeth <- dplyr::as_tibble(cbind(freq$X1, freq$X2, getmeth,cov))%>%
  dplyr::rename("chrom" = 1, "start" = 2, "meth" = 3, "cov" = 4) %>%
  mutate("meth" = as.numeric(meth)) %>%
  mutate("start" = as.numeric(start)) %>%
  mutate("chrom" = as.factor(chrom))

freq_sub <- getmeth %>%
  mutate("end" = start)

freq.gr <- GRanges(freq_sub)

keepi <- findOverlaps(freq.gr,chm13)
freq.matched <- freq.gr[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
    mcols(freq.matched),
    mcols(chm13[subjectHits(keepi)]))

freq.matched
saveRDS(freq.matched, file =paste0(dat, "/ref/chm13_200bpCpG_cov_200_CpG_num.rds"))
