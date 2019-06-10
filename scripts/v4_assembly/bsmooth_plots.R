library(bsseq)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(Sushi)

# load data 
BSchm13 <- read.bismark("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/bismark.out")

bed <- read.table("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v4_assembly/v4_all.chrX.bg")
meth_cov <- read.table("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v4_assembly/v4_meth.chrX.bg")

# set parameters depending on region
#par1
start=1563
stop=141388505
chrom="chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow"
wide=50000
smoothed=50000

#params for dxz4
#start=113868842
#stop=114116851
#chrom="chrX_fixedBionanoSV_centromereV3"
#wide=50000
#smoothed=500

#params for cenx
#start=57828561
#stop=60664792
#chrom="chrX_fixedBionanoSV_centromereV3"
#wide=500000
#smoothed=1000

# smooth data
chm13_smooth <- BSmooth(BSchm13,
                        ns=smoothed)

# plot bsmooth
pdf(file = "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v4_assembly/par1_smooth.pdf",   
    width =15,
    height = 8)

cenx_df <- data.frame(start = start, end = stop, chr = chrom)
#qroi <- data.frame(start = 113868842, end = 113968842, chr = chrom)
plotRegion(chm13_smooth,cenx_df, extend = wide, addRegions = roi)

dev.off()

# plot coverage 
pdf(file = "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v4_assembly/cenx_cov.pdf",   
    width =15,
    height = 8)

ONT <- plotBedgraph(bed,"chrX_fixedBionanoSV_centromereV3",(start-wide),(stop+wide), transparency =.5, color= SushiColors(2)(2)[1])


meth <- plotBedgraph(meth_cov,"chrX_fixedBionanoSV_centromereV3",(start-wide),(stop+wide), transparency =.5, color= SushiColors(2)(2)[2], overlay = TRUE)


labelgenome("cenx",start,stop,n=3,scale="Kb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
abline(v=c(start,stop), col="black", lwd=3, lty=2)
legend("topright",xpd=T, bty="n", inset=c(-5,0),legend=c("All ONT reads","ONT with high quality meth calls"), fill=opaque(SushiColors(2)(2)),border=SushiColors(2)(2),text.font=2,cex=1.0)

dev.off()
