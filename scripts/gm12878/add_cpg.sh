#!/bin/bash

ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/ref/chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow.fasta
bam=/kyber/Data/Nanopore/projects/nanonome/analysis/data/nanonome/pooled/bam/190516_GM12878_nanoNOMe.pooled.bam
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call
mbed=/kyber/Data/Nanopore/projects/nanonome/analysis/data/nanonome/pooled/mbed/GM12878_nanoNOMe.pooled.cpg.meth.bed.gz

if [ "$1" == "bam" ]; then
    python3 /home/gmoney/repos/nanopore-methylation-utilities/convert_bam_for_methylation.py -t 45 --verbose -b $bam \
                -c $mbed -w chrX:1-156040896 |\
                samtools sort -o ${outdir}/gm12878_meth_chrX_sorted.bam
        samtools index ${outdir}/gm12878_chrX_meth_sorted.bam
fi

if [ "$1" == "bismark" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/methylation.bed.gz > ${outdir}/bismark.out

fi
