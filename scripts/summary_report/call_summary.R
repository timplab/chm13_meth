    
#!/usr/bin/env Rscript 

library(knitr)
library(markdown)
library(rmarkdown)
library(tidyverse)

setwd('/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/summary_files')
list.files[1:2] <- list.files(pattern=".tsv")
list.files[3] <- list.files(pattern="bismark")
list.files[4] <- list.files(pattern="bedgraph")
roi <- "chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow:57828561-60664792"
smoothed <- 500
wide <- 50000

summary <- list()
print(list.files)
title="Methylation Report"
files = list(summary = paste("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/summary_files/", list.files, sep = ''))

name=substr(list.files[2],1,nchar(list.files[2])-4)

rmarkdown::render(input = '/home/gmoney/Code/chm13/chm13_meth/scripts/summary_report/summary.rmd',
output_format = "pdf_document",
output_file = paste(name, ".pdf", sep = ''),
output_dir = '/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/summary_files')

