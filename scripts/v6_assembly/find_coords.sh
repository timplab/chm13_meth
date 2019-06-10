#!/bin/bash 

# dxz4
#./find_in_chm13.sh chrX:115840395-115969110 

#par1
#./find_in_chm13.sh chrX:10001-2781479 

#qsubtel 
#tel=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/subtelomere/chrx_subtel.fasta
#ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/ref/chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow.fasta
#minimap2 $ref $tel > minimap_subtel.tsv


#IL1
#./find_in_chm13.sh chrX:155997581-156013017

# high meth region (switched search)
./find_in_chm13.sh chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow:153842058-153910699
