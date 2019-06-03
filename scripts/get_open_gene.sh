#!/bin/bash

# make bed file of region of interest
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/chrX_fixedBionanoSV_centromereV3.fasta
gr38=/home/gmoney/c9orf72/ref/chrx.fasta

printf 'chrX_fixedBionanoSV_centromereV3\t59060000\t59250000%s\n' > region_tmp.bed
bedtools getfasta -fi $ref -bed region_tmp.bed -fo open_region_tmp.fasta
minimap2 $gr38 open_region_tmp.fasta > minimap_region_ref.tsv
