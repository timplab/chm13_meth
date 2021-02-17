#!/usr/bin/Rscript

# load tidyverse
library(tidyverse)

# path to data directory and path to output figure directory
figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/TE/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

# load methylation data, save as grobj
chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  GRanges()

# load repeatmasker data
r1 <- read_tsv(paste0(dat, "/annotations/chm13.draft_v1.0.fasta_rm.bed"), col_names = c("chr", "start", "end", "rep", "width", "direction", "rep_type", "rep_class")) %>%
  dplyr::rename("chr" = 1, "start"=2, "end"=3, "rep"=4, "width"=5, "direction"=6, "rep_type"=7, "rep_class"=8)

# filter repeatmasker data
subtypes=c("L1P", "L1M")
rm <- r1 %>%
  filter(rep_type == "LINE") %>% # keep only L1 LINE elements
  filter(rep_class == "L1") %>%
  mutate(name = ifelse(grepl("L1M", rep), "L1M", rep)) %>% # make a name column that says what type of L1 (M or P)
  mutate(name = ifelse(grepl("L1P", rep), "L1P", name)) %>% 
  mutate(ID = row_number()) %>% # give a unique ID for each element
  filter(name %in% subtypes) %>% # remove LINEs that don't fit in these subtypes
  GRanges()


# use GRanges to find overlaps and save as dataframe
keepi <- findOverlaps(chm13_meth,rm)
freq.matched <- chm13_meth[queryHits(keepi)]

mcols(freq.matched) <- cbind.data.frame(
  mcols(freq.matched),
  mcols(rm[subjectHits(keepi)]))

rm.meth <- as.data.frame(freq.matched)

# group by ID, get start and end of the repeat, average methylation and number of CG sites
meth.avg <- rm.meth %>%
  group_by(ID) %>%
  mutate(start=min(start), end=max(end)) %>%
  group_by(ID,seqnames, start, end, name) %>%
  summarise(methylation=mean(methylated_frequency), num_cpg=n())

# density plot of repeat methylation, greater than and less than .5 might be good for these
ggplot(meth.avg, aes(x=methylation, fill=name))+geom_density()

# does number of CG site influence methylation - doesn't look like it, but L1P has more CGs than L1M
ggplot(meth.avg, aes(x=methylation, y=num_cpg, color=name, alpha=.2))+geom_point()
  