library(karyoploteR)
library(BSgenome.t2t.v1.0.release)
cen <- read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.cenSat_regions.bed"), col_names = c("chr", "start", "end","name")) %>%
  mutate(stain="acen") 

pdf("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures/chm13_ideogram.pdf",width=10,height=10 ) 

kp <- plotKaryotype(genome="BSgenome.t2t.v1.0.release")
kpAddBaseNumbers(kp, tick.dist = 50000000, tick.len = 10, tick.col="red", cex=1,
                 minor.tick.dist = 5000000, minor.tick.len = 5, minor.tick.col = "gray")
kpRect(kp, chr=seqnames(BSgenome.t2t.v1.0.release), x0=cen$start, x1=cen$end, y0=0, y1=1, col="red", data.panel="ideogram", border=NA)

dev.off()
