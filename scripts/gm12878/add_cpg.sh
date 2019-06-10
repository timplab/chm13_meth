#!/bin/bash

ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/ref/chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow.fasta
bam=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/data/bam/GM12878_pooled.bam
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/mcall
mbed=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/data/mcall/GM12878.cpg.meth.bed.gz

if [ "$1" == "bam" ]; then
	python3 /home/isac/Code/nanopore-methylation-utilities/convert_bam_for_methylation.py -t 45 --remove_poor --verbose -b $bam \
                -c $mbed -w chrX:1-156040895 |\
                samtools sort -o ${outdir}/gm12878_meth_chrX_sorted.bam
        samtools index ${outdir}/gm12878_chrX_meth_sorted.bam
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/methylation.bed.gz > ${outdir}/bismark.out

fi

if [ "$1" == "freq" ]; then
	/home/gmoney/repos/nanopolish/scripts/calculate_methylation_frequency.py -i /kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/data/nanopolish_calls/GM12878.chrX.cpg.meth.tsv >  /kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/data/nanopolish_calls/chrX.methylation_frequency.tsv
fi
