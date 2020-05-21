library(knitr)
library(tidyverse)
library(bsseq)

# script for taking bismark output and making into a methyl wig file for loading into IGV

dirs <- list.dirs(path = "/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2/", full.names = F, recursive = F)

dirs <- grep("tig", dirs, value = T)

for (tig in dirs){
  print(tig)
  bismark <- read_tsv(paste0("/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2/",tig,"/bismark.out"), col_names = F)
  BS <- read.bismark(paste0("/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2/",tig,"/bismark.out"))
  rawmeth <- getMeth(BS, type = "raw")
  bismark$X2 <- format(bismark$X2, scientific=F)
  bismark$X2 <- str_trim(bismark$X2, side = c("both", "left", "right"))
  
  wig <- as_tibble(cbind(bismark$X2, rawmeth)) %>%
    dplyr::rename("meth" = paste0("/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2/",tig,"/bismark.out")) %>%
    mutate(meth = as.numeric(meth)) %>%
    mutate_at(vars(meth),  funs(round(., 1)))
  
  wig$meth <- format(wig$meth, digits = 2)
  header=paste("track type=wiggle_0", ' name="CpG Methylation"', 'description="variableStep format"', "visibility=full", "autoScale=off", "viewLimits=0.0:1.0", "yLineMark=0.5")
  
  line1 = paste0("variableStep chrom=",tig," span=1")
  
  cat(header,"\n",line1,"\n",file = paste0("/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2/wig/",tig, "_raw_meth_calls.wig",  sep = ""))
  
  write.table(wig, file = paste0("/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2/wig/",tig, "_raw_meth_calls.wig"), col.names = F, row.names = F, quote = F, sep = "\t", append=TRUE)
  
}


