#!/usr/bin/env Rscript
library(tidyverse)
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures/read_length_analysis"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

chroms=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


for (i in 1:length(chroms)){
  
chr=chroms[i]

reads <- tabix_mbed(paste0(dat, "/censat/", chr, "_cen.mbed"),extcol = "motif",by = "read") 

sizes=c(0,10000,20000,30000,40000,50000)

cgcalls.list <- data.frame()
for (i in 1:length(sizes)){
  size_sel <- reads %>%
    mutate(rlen = end-start) %>%
    filter(rlen >= sizes[i]) %>%
    mutate(motif = "CG")
  
  cgcalls <- mbedByCall(size_sel) %>%
    drop_na(mcall) %>%
    group_by(start) %>%
    summarise(num_meth = sum(mcall == 1), num_unmeth = sum(mcall == 0)) %>%
    mutate(rlen = sizes[i])
  cgcalls.list <- rbind(cgcalls, cgcalls.list)
  
}


cgcalls.list.sum <- cgcalls.list %>%
  mutate(cov = num_meth+num_unmeth, freq = num_meth/cov) %>%
  group_by(rlen) %>% 
  mutate(z_score = (cov - mean(cov)) / sd(cov))

ggplot(cgcalls.list.sum, aes(y=z_score, x=as.factor(rlen), fill=as.factor(rlen)))+geom_violin()+geom_boxplot(width=.1)+theme_classic()+geom_hline(yintercept=c(3,-3), linetype = "dashed")+theme(text = element_text(size=20))+labs(ylab = "CG call Z-score", xlab = "Read Length")

# find percentage that is more than 3 standard deviations from the mean

cgcalls.list.sum <- cgcalls.list %>%
  mutate(cov = num_meth+num_unmeth, freq = num_meth/cov) %>%
  group_by(rlen) %>% 
  mutate(z_score = (cov - mean(cov)) / sd(cov))

qual <- cgcalls.list.sum %>%
  summarise( bad= sum(z_score > 3 | z_score < -3 ), good = sum(z_score <= 3 | z_score <= -3)) %>%
  mutate(total = good+bad, perc = (bad/total)*100)

ggsave(
  paste0(figs, "/", chr, "_CoverageZscore.pdf"),
  plot = last_plot(),
  scale = 1,
  width = 6,
  height = 5,
)

write.table(qual, file = paste0(figs, "/", chr, "_CoverageScore.tsv" ), quote=F, sep = "\t", row.names = F, col.names = T)
write.table(cgcalls.list, file = paste0(figs, "/", chr, "_CGCoverage.tsv" ), quote=F, sep = "\t", row.names = F, col.names = T)
}

