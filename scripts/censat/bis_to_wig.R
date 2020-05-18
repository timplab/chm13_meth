library(knitr)
library(tidyverse)
library(bsseq)

# script for taking bismark output and making into a methyl wig file for loading into IGV

dirs <- list.dirs(path = "/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/meth_calls", full.names = F, recursive = F)

dirs <- grep("tig", dirs, value = T)

for (tig in dirs){
  print(tig)
  BS <- read.bismark(paste0("/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/meth_calls/",tig,"/bismark.out"))
  rawmeth <- getMeth(BS, type = "raw")
  bismark$X2 <- format(bismark$X2, scientific=F)
  bismark$X2 <- str_trim(bismark$X2, side = c("both", "left", "right"))
  rawmeth <- as_tibble(cbind(bismark$X2, rawmeth))
  
  header=paste0("variableStep chrom=",tig," span=1")
  
  cat(header,"\n", file = paste0("/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/meth_calls/wig/",tig, "_raw_meth_calls.wig"))
  
  write.table(rawmeth, file = paste0("/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/meth_calls/wig/",tig, "_raw_meth_calls.wig"), col.names = F, row.names = F, quote = F, sep = "\t", append=TRUE)
  
}


