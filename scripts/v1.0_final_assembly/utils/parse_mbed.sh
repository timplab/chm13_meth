#!/bin/bash 
in=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/annotations
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/censat
file=${in}/t2t-chm13.v1.0.cenSat_regions.bed
tab=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/methylation_calls/CpGmethylation.bed.gz
while IFS= read -r line; do
	coord=$(echo "$line" | awk '{print $1,":",$2,"-",$3}' | awk '{ sub(/^[ \t]+/, ""); print}' | sed -e "s/ //g")
	chr=$(echo "$line" | awk '{print $1}')
	echo $chr
	tabix $tab $coord > ${out}/${chr}_cen.mbed

done < $file
