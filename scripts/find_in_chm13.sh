#!/bin/bash

# script takes coordinates from grch38 and outputs that region in jain gm12878 assembly

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

bedtools getfasta -fi $gr38 -bed region_tmp.bed -fo region.fasta
minimap2 $ref region.fasta | awk '($12>0)' > minimap.tsv
newcoord1=$(cut -f8 minimap.tsv | sort -n | head -n 1)
newcoord2=$(cut -f9 minimap.tsv | sort -n | tail -n 1)
newchr=$(cut -f6 minimap_coords.tsv | head -n 1)

echo "coordinates in chm13 assebly:" 
echo $newchr 
echo $newcoord1 
echo $newcoord2
echo "length of qry region:"
echo "$(($coord1-$coord2))"
echo "length of chm13 region:"
echo "$(($newcoord2-$newcoord1))"

#rm *tmp*
