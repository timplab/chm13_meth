source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/ilee_plot_utils.R")
library("ggsci")
source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/methylation_R_utils.R")
library(tidyverse)
library(GenomicRanges)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
list=c("chr13", "chr14", "chr15", "chr21", "chr22")



#22:5749418-5966044
chrom="chr22"
zstart=3634322
zend=4793756

issues <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/chm13_v1.0_lifted.issuesandchm13gaps.bed", col_names = F) %>%
  dplyr::rename("seqnames"=X1, "start"=X2, "end"=X3) %>%
  filter(seqnames == chrom) 
unmapped <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/chm13_v1.0_unmapped_byGRCh38.bed", col_names = F) %>%
  dplyr::rename("seqnames"=X1, "start"=X2, "end"=X3) %>%
  filter(seqnames == chrom)

transcriptColors =c("coding" = "navy",
                "lncRNA"="darkgreen", 
                "pseudo" = "#cc3399", 
                "smallRNA" = "goldenrod", 
                "unsure" = "blueviolet")


cat <- paste0(dat, "/chm13_final_beds/CHM13_CAT_all_exons.bed")
cat.gr <- import(cat)
exp <- read_tsv(paste0(dat, "/chm13_final_beds/CHM13.category-colors.bed.gz"), col_names = c("chr","start", "end", "name", "num", "strand", "num2", "num3", "color")) %>%
  mutate(type=case_when(color == "0,100,0" ~ "lncRNA", 
                        color == "0,0,128" ~ "coding", 
                        color == "255,0,255" ~ "pseudo", 
                        color == "218,165,32" ~ "smallRNA", 
                        color =="138,43,226" ~ "unsure")) %>%
  dplyr::select(c(name, type))


cat.sub <- cat.gr[seqnames(cat.gr) == chrom & start(cat.gr) > zstart & end(cat.gr) < zend]

transcript.df <- as.data.frame(cat.sub) %>%
  group_by(name) %>%
  summarise(start=min(start), end=max(end)) %>%
  mutate(y = seq_len(length(unique((as.data.frame(cat.sub))$name))))
transcript.df <- merge(transcript.df, exp, by = "name")

exon.df <- merge(transcript.df, as.data.frame(cat.sub),by = "name")
exon.df <- merge(exon.df, exp, by="name")




reads <- tabix_mbed(paste0(dat, "/censat/", chrom, "_cen.mbed"),extcol = "motif",by = "read") 
#129004453       142241828 



size_sel <- reads %>%
  mutate(rlen = end-start) %>%
  filter(rlen >= 100000) %>%
  mutate(motif = "CG")

cgcalls <- mbedByCall(size_sel) %>%
  drop_na(mcall)

reg = cgcalls %>%
  dplyr::filter(start >= zstart) %>%
  dplyr::filter(end <= zend)


cpg_runs <-getRuns(reg, maxGap = 500)


cpg_runs.ordered <- order_reads(cpg_runs)

cpg_runs <- cpg_runs.ordered$x %>%
  mutate(m = ifelse(values == 1, "Methylated","Unmethylated")) %>%
  mutate(mod = "CpG")

pal <- pal_npg("nrc")(10)
meth_pal <- c(pal[8],pal[4]) 
#ids <- order(sapply(seq_along(gcruns.list),function(i){length(unique(gcruns.list$qname))}))


bis <- reg %>% 
  group_by(chrom, start) %>%
  summarise(meth_freq = mean(mcall), cov = n()) 

meth <- ggplot(bis, aes(x = (start/1e6), y= meth_freq))+geom_smooth(size =1, se = F, method = "loess", span = .1)+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank(),panel.border=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+scale_y_continuous(breaks=seq(0, 1.50, .50),limits = c(0, 1))+labs(y="Methylation Frequency")+coord_cartesian(xlim=c(zstart/1e6,zend/1e6))


g <- ggplot(cpg_runs,aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax)) + geom_rect(data = cpg_runs.ordered$bounds, fill = "grey80") + 
  geom_rect(data=cpg_runs,aes(fill = m))+scale_fill_manual(name = "State", values = meth_pal) + theme(axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),legend.position= "none",axis.ticks.y = element_blank(), panel.spacing = unit(2, "lines"))+theme(legend.title=element_blank(),panel.border=element_blank())+coord_cartesian(xlim=c(zstart,zend))

plot <- ggplot()+
  geom_rect(data=transcript.df, mapping=aes(xmin=start/1e6,xmax=end/1e6,ymin=y,ymax=y+.1,fill=type))+theme(legend.position="none")+theme(legend.text=element_text(size=rel(1)))+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+geom_rect(exon.df, aes(xmin=start.y/1e6, xmax=end.y/1e6, ymin=y-.3, ymax=y+.3,fill=type.y))+theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+scale_fill_manual(values = transcriptColors, drop = FALSE)+theme(legend.title=element_blank(),panel.border=element_blank())+coord_cartesian(xlim=c(zstart/1e6,zend/1e6))
plot

map.plot <- ggplot()+geom_rect(data=issues, aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=X4, color="red", alpha=.5, linetype="dashed"))+scale_color_manual(values = censatColors, drop = FALSE)+geom_rect(data=unmapped, aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill="black", color="black", alpha=.2, linetype="dashed"))+theme(legend.title=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+scale_fill_manual(values = repeatColors, drop = FALSE)+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.border=element_blank(), legend.position = "none")+coord_cartesian(xlim=c(zstart/1e6,zend/1e6))


row <- plot_grid(meth,g,plot, map.plot,align="v", ncol=1,rel_heights = c(1/8,1/4,1/8,1/20))
row

ggsave(
  paste0(figs, "/methyl_profiles/", chrom, "-", zstart, "-", zend, "_GeneModelPlot.pdf"),
  plot = row,
  scale = 1,
  width = 8,
  height = 6
)