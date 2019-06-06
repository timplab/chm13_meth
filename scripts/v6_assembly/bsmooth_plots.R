library(bsseq)
library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(Sushi)

# load data 
BSchm13 <- read.bismark("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/bismark.out")
bed <- read.table("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call/chm13_v0.6.nanopolish2.arrow2_10x2.chrX.bg")

# set parameters depending on region
start=1563
stop=141388505
chrom="chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow"
wide=50000
smoothed=50000

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
pdf(file = "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/figs/par1_smooth1000.pdf",   
    width =15,
    height = 8)

cenx_df <- data.frame(start = start, end = stop, chr = chrom)
plotRegion(chm13_smooth,cenx_df, extend = wide, addRegions = cenx_df)

dev.off()

# plot coverage 
pdf(file = "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/figs/par1_cov.pdf",   
    width =15,
    height = 8)

plotBedgraph(bed,"chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow",(start-wide),(stop+wide),colorbycol= SushiColors("darkblue"))
labelgenome("cenx",start,stop,n=3,scale="Kb")
mtext("Read Depth",side=2,line=1.75,cex=1,font=2)
axis(side=2,las=2,tcl=.2)
abline(v=c(start,stop), col="red", lwd=3, lty=2)
dev.off()
