library(bsseq)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(Sushi)

# load data 
BSchm13 <- read.bismark("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/bismark.out")
bed <- read.table("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/chm13_v0.6.nanopolish2.arrow2_10x2.chrX.bg")
meth_cov <- read.table("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/chm13_meth.chrX.bg")

# set parameters depending on region
#par1
start=1563 
stop=2601563 
chrom="chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow"
wide=0
smoothed=1000

#params for dxz4
#start=113868842
#stop=114116851
#chrom="chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow"
#wide=50000
#smoothed=500

#params for cenx
#start=57828561
#stop=60664792
#chrom="chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow"
#wide=500000
#smoothed=1000

# smooth data
chm13_smooth <- BSmooth(BSchm13,
                        ns=smoothed)

# plot bsmooth
pdf(file = "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/figs/par1_smooth.pdf",   
    width =15,
    height = 8)

cenx_df <- data.frame(start = start, end = stop, chr = chrom)
#roi <- data.frame(start = 113868842, end = 113968842, chr = chrom)
plotRegion(chm13_smooth,cenx_df, extend = wide, addRegions = cenx_df)

dev.off()

# plot coverage 
pdf(file = "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/figs/par1_cov.pdf",   
    width =15,
    height = 8)

ONT <- plotBedgraph(bed,"chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow",(start-wide),(stop+wide), transparency =.5, color= SushiColors(2)(2)[1])


meth <- plotBedgraph(meth_cov,"chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow",(start-wide),(stop+wide), transparency =.5, color= SushiColors(2)(2)[2], overlay = TRUE)


labelgenome("cenx",start,stop,n=3,scale="Kb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
abline(v=c(start,stop), col="black", lwd=3, lty=2)
legend("topright",xpd=T, bty="n", inset=c(-5,0),legend=c("All ONT reads","ONT with high quality meth calls"), fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2,cex=1.0)

dev.off()
