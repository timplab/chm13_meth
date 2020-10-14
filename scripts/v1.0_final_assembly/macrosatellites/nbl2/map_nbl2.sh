#!/bin/bash 
seq=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/macrosatellites/nbl2/nbl2.fasta
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/macrosatellites/nbl2
minimap2 -cx asm5 $seq $ref > ${out}/nbl2_chm13_minimap.paf
