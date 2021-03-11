
library(tidyverse)
library(zoo)
library(stats)
source("/home/isac/Code/ilee/plot/ilee_plot_utils.R")
library("ggsci")
source("/home/isac/Code/nanopore-methylation-utilities/methylation_R_utils.R")


figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/censat/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"


chr="chr15"

reads <- tabix_mbed(paste0(dat, "/censat/", chr, "_cen.mbed"),extcol = "motif",by = "read") 
#129004453       142241828 



size_sel <- reads %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= 50000) %>%
  mutate(motif = "CG")

cgcalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)

region_start=10000000
region_end= region_start+50000
reg = cgcalls %>%
  dplyr::filter(start >= region_start) %>%
  dplyr::filter(end <= region_end)

bis <- reg %>% 
  group_by(chrom, start) %>%
  summarise(meth_freq = mean(mcall), cov = n()) 


n_windows=1000
bis$cut = cut(bis$start, breaks=n_windows)

meth.stat <- bis %>%
  group_by(cut) %>%
  summarise(med = median(meth_freq), top = quantile(meth_freq, 0.75), bot = quantile(meth_freq, 0.25), n_genes = length(meth_freq)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  mutate(med_smooth=rollmean(med, 10, NA))


meth <- meth.stat %>%
  select(c(min,med_smooth)) %>%
  column_to_rownames("min") %>%
  na.omit() %>%
  as.matrix()

colnames(meth)<-NULL
pdf(paste0(figs, "/", chr, "_autocorrelationPlot.pdf")) 
acf(meth,lag=100)
dev.off()



chr="chr1"

reads <- tabix_mbed(paste0(dat, "/censat/", chr, "_cen.mbed"),extcol = "motif",by = "read") 
#129004453       142241828 



size_sel <- reads %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= 50000) %>%
  mutate(motif = "CG")

cgcalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)

region_start=130000000
region_end= region_start+50000

reg = cgcalls %>%
  dplyr::filter(start >= region_start) %>%
  dplyr::filter(end <= region_end)

bis <- reg %>% 
  group_by(chrom, start) %>%
  summarise(meth_freq = mean(mcall), cov = n()) 


n_windows=1000
bis$cut = cut(bis$start, breaks=n_windows)

meth.stat <- bis %>%
  group_by(cut) %>%
  summarise(med = median(meth_freq), top = quantile(meth_freq, 0.75), bot = quantile(meth_freq, 0.25), n_genes = length(meth_freq)) %>%
  mutate(x_tmp = str_sub(cut, 2, -2)) %>% 
  separate(x_tmp, c("min", "max"), sep = ",") %>% 
  mutate_at(c("min", "max"), as.double) %>%
  mutate(med_smooth=rollmean(med, 10, NA))

meth <- meth.stat %>%
  select(c(min,med_smooth)) %>%
  column_to_rownames("min") %>%
  na.omit() %>%
  as.matrix()
colnames(meth)<-NULL

pdf(paste0(figs, "/", chr, "_autocorrelationPlot.pdf")) 
acf(meth,lag=100)
dev.off()

