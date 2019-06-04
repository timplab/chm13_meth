#!/bin/bash

# script takes coordinates from chrX grch38 and outputs that region in chrX jain gm12878 assembly

ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/chrX_fixedBionanoSV_centromereV3.fasta
gr38=/kyber/Data/Nanopore/Analysis/gmoney/ref/GRCh38.primary_assembly.genome.fa

if [ $# -gt 0 ]; then
	echo "command contains $# arguments"
else 
	echo "command contains no arguments, please include coordinates ie chr#:xxxx-xxxx"
fi
coords=$1

coord1=$(echo $coords | awk '{split($0,a, "-");print a[2]}') 
chr=$(echo $coords | awk '{split($0,a, ":");print a[1]}')
coord2=$(echo $coords | awk '{split($0,a, ":");print a[2]}' |  awk '{split($0,a, "-");print a[1]}')


printf '%s\t%s\t%s\n' $chr $coord2 $coord1 > region_tmp.bed

bedtools getfasta -fi $gr38 -bed region_tmp.bed -fo region_tmp.fasta
minimap2 $ref region_tmp.fasta | awk '($12>0)' > minimap_tmp.tsv
#chroms=$(cut -f6 minimap_tmp.tsv | uniq)

newcoord1=$(cut -f8 minimap_tmp.tsv | sort -n | head -n 1)
newcoord2=$(cut -f9 minimap_tmp.tsv | sort -n | tail -n 1)
newchr=$(cut -f6 minimap_tmp.tsv | head -n 1)

echo "coordinates in chm13 assebly:" 
echo $newchr 
echo $newcoord1 
echo $newcoord2
echo "length of qry region:"
echo "$(($coord1-$coord2))"
echo "length of chm13 region:"
echo "$(($newcoord2-$newcoord1))"

#rm *tmp*
