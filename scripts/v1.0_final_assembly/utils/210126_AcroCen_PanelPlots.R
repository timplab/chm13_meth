source("~/projects/chm13_meth/scripts/v1.0_final_assembly/utils/ilee_plot_utils.R")
library(tidyverse)
library(cowplot)

figs="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/figures"
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"
list=c("chr13", "chr14", "chr15", "chr21", "chr22")


cen <- read_tsv(paste0(dat, "/chm13_final_beds/cenRegions.bed"), col_names = c("chr", "start", "end","name")) %>%
  filter(chr %in% list)

for (i in 1:length(cen$chr)){
chrom=cen$chr[i]
rstart=cen$start[i]
rend=cen$end[i]

all.dat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/210126_AcroCen_10kbBinned_ALL.bed") %>%
  filter(seqnames == chrom) %>%
  filter(start >= rstart) %>%
  filter(end <= rend)

repeatColors =c("DNA"="#C19936",
                "DNA?"="#C19935",
                "LINE"="#FFA500",
                "Low_complexity"="#75B043",
                "LTR"="#51B756",
                "RC"="#53BB74",
                "Retroposon"="#55BE9D",
                "RNA"="#ff4000",
                "rRNA"="#52BEBB",
                "scRNA"="#6B9AD3",
                "Simple_repeat"="#8992C9",
                "SINE"="#9A8AC1",
                "snRNA"="#A886BC",
                "srpRNA"="#B67EB6",
                "Unspecified"="#C378B2",
                "tRNA"="#006400",
                "Unknown"="#BA55D3",
                "Satellite"="#53B0E4",
                "red"="black", 
                "back"="black", 
                "#647FA4" = "black")

defaultColor = "#000080"

censatColors =c("TE" = "#E87C71",
                "MER"="#E28455",
                "HOR"="#D78C32",
                "BSAT"="#E370AB",
                "CER" = "#CE9334",
                "HSAT2"="#C19935",
                "HSAT1"="#A2A638",
                "HSAT3"="#8CAC3E",
                "L1"="#75B042",
                "LSAU"="#54B346",
                "LTR"="#51B756",
                "MST"="#53BB73",
                "GSAT"="#55BE8D",
                "GSATII"="#54C0A5",
                "rRNA"="#52BEBB",
                "SAR"="#51BDCE",
                "novel"="#9400D3",
                "HSAT4"="#53B0E4",
                "SATR"="#5AA5DA",
                "CT"="#6B9AD2",
                "HERV"="#8992C8",
                "MSAT"="#9A8AC2",
                "MON"="#A885BC",
                "SST"="#C378B2",
                "HSAT5"="#ED72A5",
                "HSAT6"="#EF768C", 
                "gap-rDNA"="#ff4000",
                "L1" = "#ffbf00", 
                "TAR"= "#0080ff",
                "ACRO"="#9400D4",
                "Alu"="#9A8AC3")

sats=c("HSAT1", "HSAT2", "HSAT3", "HSAT4", "HSAT5", "HSAT6", "CT", "BSAT", "HOR", "MON", "Alu", "SATR", "ACRO", "GSATII", "TAR","TE","MER","MST","CER","L1","SST","LSAU","GSAT","MSAT","novel","HERV","LTR", "gap-rDNA")
censat = read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.1110.cenSat_annotationTAB.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("HSat1", X4), "HSAT1", X4)) %>%
  mutate(name = ifelse(grepl("HSat2", X4), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("HSat3", X4), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("HSat4", X4), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("HSat5", X4), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("HSat6", X4), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", X4), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", X4), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("hor", X4), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", X4), "MON", name)) %>%
  mutate(name = ifelse(grepl("Alu", X4), "Alu", name)) %>%
  mutate(name = ifelse(grepl("SATR", X4), "SATR", name)) %>%
  mutate(name = ifelse(grepl("ACRO", X4), "ACRO", name)) %>%
  mutate(name = ifelse(grepl("GSATII", X4), "GSATII", name)) %>%
  mutate(name = ifelse(grepl("TAR", X4), "TAR", name)) %>%
  mutate(name = ifelse(grepl("TE", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("MER", X4), "MER", name)) %>%
  mutate(name = ifelse(grepl("MST", X4), "MST", name)) %>%
  mutate(name = ifelse(grepl("CER", X4), "CER", name)) %>%
  mutate(name = ifelse(grepl("L1", X4), "L1", name)) %>%
  mutate(name = ifelse(grepl("SST", X4), "SST", name)) %>%
  mutate(name = ifelse(grepl("LSAU", X4), "LSAU", name)) %>%
  mutate(name = ifelse(grepl("GSAT", X4), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("MSAT", X4), "MSAT", name)) %>%
  mutate(name = ifelse(grepl("novel", X4), "novel", name)) %>%
  mutate(name = ifelse(grepl("HERV", X4), "HERV", name)) %>%
  mutate(name = ifelse(grepl("LTR", X4), "LTR", name)) %>%
  dplyr::select(c(X1, X2, X3, name)) %>%
  dplyr::filter(name %in% sats) %>%
  dplyr::rename("chr" =1, "start" = 2 ,"end" =3) 

issues <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/chm13_v1.0_lifted.issuesandchm13gaps.bed", col_names = F) %>%
  dplyr::rename("seqnames"=X1, "start"=X2, "end"=X3) %>%
  filter(seqnames == chrom) 
unmapped <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/chm13_v1.0_unmapped_byGRCh38.bed", col_names = F) %>%
  dplyr::rename("seqnames"=X1, "start"=X2, "end"=X3) %>%
  filter(seqnames == chrom)


meth <- ggplot(all.dat, aes(x = (start/1e6), y= smooth))+geom_line(size =1) + labs(y="Methylation")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+ylim(0,1)

kmer.plot <- ggplot(all.dat, aes(x = (start/1e6), y= kmer_num))+geom_bar(size =1, stat="identity",color='#1B9E77') + labs(y="51mers")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

genes.plot <- ggplot()+geom_bar(all.dat, aes(x = (start/1e6), y= gene_num),size =1, stat="identity",color='#D95F02') + labs(y="Genes")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+geom_bar(all.dat, aes(x = (start/1e6), y= iso_num, alpha=.5),size =1, stat="identity",color='#7570B3', alpha=.5)+geom_bar(all.dat, aes(x = (start/1e6), y=cat_num),size =1, stat="identity",color='#E7298A')

pro.plot <- ggplot()+geom_bar(all.dat, aes(x = (start/1e6), y= proA_cov),size =1, stat="identity",color='#E6AB02') + labs(y="PRO-Seq",x= "Genomic Coordinates")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+geom_bar(all.dat, aes(x = (start/1e6), y= proB_cov),size =1, stat="identity",color='#66A61E')


sd.plot <- ggplot(all.dat, aes(x = (start/1e6), y= sd_num))+geom_bar(size =1, stat="identity",color='darkgreen') + labs(y="SDs")+theme_classic(base_size = 10)+xlim((rstart/1e6),(rend/1e6))+ theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


censat_sub <- censat %>%
  dplyr::filter(chr == chrom) %>%
  filter(start >= rstart) %>%
  filter(end <= rend)

censat.plot <- ggplot(data=censat_sub, mapping=aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=name))+
  geom_rect()+theme(legend.position="none")+theme(legend.text=element_text(size=rel(1)))+theme(legend.title=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))+scale_fill_manual(values = censatColors, drop = FALSE)+theme(legend.title=element_blank(),panel.border=element_blank())+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.border=element_blank())

map.plot <- ggplot()+geom_rect(data=issues, aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1), fill ="black",color="black", alpha=.5)+geom_rect(data=unmapped, aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1), fill = "red",color="red", alpha=.5)+theme(legend.title=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+theme(legend.title=element_blank(),panel.border=element_blank())+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))+scale_fill_manual(values = repeatColors, drop = FALSE)

gc.plot <- ggplot()+geom_rect(all.dat, aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=gc.content*100))+scale_fill_gradient2(low="#132B43", high="dodgerblue", limits=c(20,80))+theme(legend.title=element_blank())+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))+theme(legend.title=element_blank(),panel.border=element_blank())+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.border=element_blank(), legend.position = "none")

top_row <- plot_grid(gc.plot, censat.plot,meth,genes.plot,sd.plot,pro.plot,ncol = 1, align="v",rel_heights = c(1/25,1/25,1/6,1/6,1/6,1/6))

top_row
ggsave(
  paste0(figs, "/methyl_profiles/", chrom, "AcroCen_PanelPlot.pdf"),
  plot = top_row,
  scale = 1,
  width = 10,
  height = 5
)}

# legends
all.dat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/210126_AcroCen_10kbBinned_ALL.bed") 

gc.plot <- ggplot()+geom_rect(all.dat, aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=gc.content*100))+scale_fill_gradient2(low="#132B43", high="dodgerblue", limits=c(20,80))+theme(legend.title=element_blank())+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))+theme(legend.title=element_blank(),panel.border=element_blank())+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),panel.border=element_blank(), legend.position = "bottom")


censat.all <- ggplot(data=censat, mapping=aes(xmin=(start/1e6),xmax=(end/1e6),ymin=0,ymax=.1,fill=name))+
  geom_rect()+theme(legend.position="bottom") +labs(y="CenSat")+theme(legend.text=element_text(size=rel(1)))+theme(legend.title=element_blank()) +theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+coord_cartesian(xlim=c(rstart/1e6,rend/1e6))+scale_fill_manual(values = censatColors, drop = FALSE)



legend <- cowplot::get_legend(gc.plot)
legend2 <- cowplot::get_legend(censat.all)
leg <- plot_grid(NULL, legend, legend2, ncol=1)

ggsave(
  paste0(figs, "/methyl_profiles/AcroCen_PanelPlotLegend.pdf"),
  plot = leg,
  scale = 1,
  width = 10,
  height = 10
)

all.dat <- read_tsv("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/210126_AllCen_10kbBinned_ALL.bed") 

gc.bed <- all.dat %>%
  select(c(seqnames, start,end, gc.content)) %>%
  mutate(gc.content=round(gc.content,2))

write.table(gc.bed, "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/chm13_final_beds/chm13_v1.0_GCcontent_10kb.bed", quote=F, sep="\t", col.names = F, row.names = F)
