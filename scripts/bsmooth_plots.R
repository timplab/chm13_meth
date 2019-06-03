library(bsseq)
library(Biostrings)
library(GenomicRanges)
library(Gviz)
library(ggbio)

BSchm13 <- read.bismark("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/bismark.out")
chm13_smooth <- BSmooth(BSchm13,
                        ns=1000)
cenx_df <- data.frame(start = 57812643, end = 60534379, chr = "chrX_fixedBionanoSV_centromereV3")
cenx_meth <- plotRegion(chm13_smooth,cenx_df, extend = 500000, addRegions = cenx_df)

chm13_smooth <- BSmooth(BSchm13,
                        ns=200)
cenx_open <- data.frame(start = 59060000, end = 59250000, chr = "chrX_fixedBionanoSV_centromereV3")
cenx_open <- plotRegion(chm13_smooth, cenx_open, extend = 50000, addRegions = cenx_open)

#113740520-113987898
chm13_smooth <- BSmooth(BSchm13,
                        ns=500)
dxz4_df <- data.frame(start = 113740520, end = 113987898, chr = "chrX_fixedBionanoSV_centromereV3")
dxz4_meth <- plotRegion(chm13_smooth, dxz4_df, extend = 50000, addRegions = dxz4_df)

chm13_smooth <- BSmooth(BSchm13,
                        ns=50000)
all_df <- data.frame(start = 0, end = 153867332, chr = "chrX_fixedBionanoSV_centromereV3")
regions_df <- rbind.data.frame(dxz4_df, cenx_df)
all_meth <- plotRegion(chm13_smooth, all_df, addRegions = regions_df)


# can't do custom ideotracks without cytoband info :(  -- using hg19 
itrack <- IdeogramTrack(genome = "hg19", chromosome = "chrX")
gtrack <- GenomeAxisTrack()
dxz4 <- plotTracks(itrack, from = 113740520, to = 113987898)
cenx <- plotTracks(itrack, from = 57812643, to = 60534379)

