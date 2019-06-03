#!/bin/bash

# make bed file of region of interest
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/chrX_fixedBionanoSV_centromereV3.fasta
gr38=/home/gmoney/c9orf72/ref/GRCh38.primary_assembly.genome.fa

printf 'chrX\t115840395\t115969110\n' > region_tmp.bed
bedtools getfasta -fi $gr38 -bed region_tmp.bed -fo region_tmp.fasta
minimap2 $ref region_tmp.fasta > minimap_dzx4_ref.tsv
#rm *tmp*
