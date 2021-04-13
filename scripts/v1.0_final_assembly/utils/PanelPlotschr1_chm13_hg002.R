
source("~/T2T-Epigenetics/utils/ilee_plot_utils.R")
source("~/T2T-Epigenetics/utils/methylation_R_utils.R")
library(tidyverse)
library(cowplot)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

# read methylation GRanges data
hg002_meth <- read_tsv(paste0(dat, "/HG002/nanonome/methylation_calls/whole_genome/chm13_hg002_reference_pooled/HG002_nanonome_CpGmethylationFrequency_20kb.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  filter(chromosome != "chrY") %>%
  GRanges()

chm13_meth <- read_tsv(paste0(dat, "/methylation_calls/methylation_frequency_50kb_split.tsv")) %>%
  mutate(called_sites_unmethylated = called_sites - called_sites_methylated) %>%
  GRanges()

cen <- read_tsv(paste0(dat, "/annotations/t2t_cenRegions.v2.021621.bed"), col_names = c("chr", "start", "end","name")) %>%
  group_by(chr) %>%
  summarise(start=min(start), end=max(end))

i=20
chrom=cen$chr[i]
rstart=cen$start[i]
rend=cen$end[i]


high=quantile(hg002_meth$called_sites, .99)
low=quantile(hg002_meth$called_sites, .02)

flags.hg002 <- as.data.frame(hg002_meth) %>% 
  filter(called_sites < 10 | called_sites > 100) %>%
  mutate(line="HG002") %>%
  filter(seqnames == chrom) %>%
  mutate(num=1) %>%
  GRanges()

high=quantile(chm13_meth$called_sites, .99)
low=quantile(chm13_meth$called_sites, .02)

flags.chm13 <- as.data.frame(chm13_meth) %>% 
  filter(called_sites < 10 | called_sites > 100) %>%
  mutate(line="CHM13") %>%
  mutate(num=1) %>%
  GRanges()


blocks <- genomeBlocks(BSgenome.t2t.v1.0.release, chrs = seqlevels(BSgenome.t2t.v1.0.release), width = 10000)
score1 <- coverage(hg002_meth, weight="called_sites_methylated")
score2 <- coverage(hg002_meth, weight="called_sites_unmethylated")
score3 <- coverage(hg002_meth, weight="num_motifs_in_group")
score4 <- coverage(flags.hg002, weight="num")
score5 <- coverage(flags.chm13, weight="num")
score6 <- coverage(flags.bsseq, weight="num")


binned_meth <- binnedSum(blocks, numvar = score1, "called_sites_methylated") %>%
  as.data.frame()

binned_unmeth <-binnedSum(blocks, numvar = score2, "called_sites_unmethylated")%>%
  as.data.frame()

binned_cov <- binnedSum(blocks, numvar = score3, "num_motifs_in_group") %>%
  as.data.frame()

binned_flag <- binnedSum(blocks, numvar = score4, "num")%>%
  as.data.frame()
binned_flag.chm <- binnedSum(blocks, numvar = score5, "num")%>%
  as.data.frame() %>%
  dplyr::rename("num_chm13"=num)
binned_flag.bsseq <- binnedSum(blocks, numvar = score6, "num")%>%
  as.data.frame() %>%
  dplyr::rename("num_bsseq"=num)

meth_bins <- merge(binned_meth, binned_unmeth, by = c("start", "end", "seqnames", "width", "strand")) %>%
  merge(binned_cov, by = c("start", "end", "seqnames", "width", "strand")) %>%
  merge(binned_flag, by = c("start", "end", "seqnames", "width", "strand")) %>%
  merge(binned_flag.chm, by = c("start", "end", "seqnames", "width", "strand")) %>%
  filter(seqnames == chrom) %>%
  # filter(start > rstart) %>%
  # filter(end < rend) %>%
  group_by(start, end, seqnames) %>%
  mutate(sites = called_sites_methylated+called_sites_unmethylated) %>%
  mutate(freq = called_sites_methylated/sites) %>%
  ungroup() %>%
  group_by(seqnames) %>%
  arrange(start,seqnames) %>%
  mutate(hg002_smooth = rollmean(freq, 3, fill = NA), HG002_site_smooth = rollmean(num_motifs_in_group, 3, fill = NA)) %>%
  ungroup() %>%
  arrange(seqnames, start) %>%
  mutate(perc_bad.hg002=num/num_motifs_in_group) %>%
  mutate(perc_bad.chm13=num_chm13/num_motifs_in_group) %>%
  select(hg002_smooth,HG002_site_smooth,perc_bad.hg002,perc_bad.chm13) 



all.dat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/210126_Allchr_10kbBinned_ALL.bed") %>%
  filter(seqnames == chrom) 

binned.all <- cbind(all.dat, meth_bins)%>%
  filter(start > rstart) %>%
  filter(end < rend)

SAT = c("GSAT", "DHOR", "BSAT","HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6","HOR", "MON", "CT", "gap-rDNA")
censat = read_tsv(paste0(dat, "/annotations/t2t_cenAnnotation.v2.021621FORMATTED.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("gap", X4), "gap-rDNA", X4)) %>%
  dplyr::filter(X4 %in% SAT) %>%
  dplyr::rename("chr" =1, "start" = 2 ,"end" =3) 

binned.hg002 <- binned.all %>%
  filter(perc_bad.hg002 < .9)

binned.chm13 <- binned.all %>%
  filter(perc_bad.chm13 < .9)

meth <- ggplot()+geom_line(data=binned.all, aes(x = start/1e6, y= smooth),size =1,color="dodgerblue")+geom_line(data=binned.all, aes(x = (start/1e6), y= hg002_smooth),size =1,color="purple")  + labs(y="Methylation")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+ylim(0,1)

cov <- ggplot()+geom_line(data=binned.all, aes(x = start/1e6, y= ont_cov),size =1,color="dodgerblue")+geom_line(data=binned.all, aes(x = (start/1e6), y= HG002_site_smooth),size =1,color="purple")  + labs(y="Methylation")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

censat_sub <- censat %>%
  dplyr::filter(chr == chrom) 

censat.plot <- ggplot(data=censat_sub, mapping=aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=name))+
  geom_rect()+theme(legend.position="none") +labs(y="CenSat")+theme(legend.text=element_text(size=rel(1)))+theme(legend.title=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))+scale_fill_manual(values = censatColors, drop = FALSE)+theme(legend.title=element_blank(),panel.border=element_blank())

hg002.bad <- binned.all %>%
  filter(seqnames == chrom) %>%
  filter(perc_bad.hg002 > .9)
  
chm13.bad <- binned.all %>%
  filter(seqnames == chrom) %>%
  filter(perc_bad.chm13 > .9)


map.plot1 <- ggplot()+geom_point(data=hg002.bad,aes(x=(start/1e6),y=.1, alpha=.5),color="purple")+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.border=element_blank(), legend.position = "none",axis.title.x=element_blank())+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))
map.plot2 <- ggplot()+geom_point(data=chm13.bad,aes(x=(start/1e6),y=.1, alpha=.5),color="dodgerblue")+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.border=element_blank(), legend.position = "none",axis.title.x=element_blank())+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))

top_row <- plot_grid(meth,map.plot2, map.plot1,censat.plot,ncol = 1, align="v",rel_heights = c(1/3,1/25,1/25,1/20))
top_row
ggsave(
  paste0(figs, "/methyl_profiles/", chrom,"_PanelPlotv2_chm13_hg002_mapstat.pdf"),
  plot = top_row,
  scale = 1,
  width = 10,
  height = 3
)


