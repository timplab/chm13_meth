#!/bin/bash
root=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/chr8
ref=${root}/data/chm13_chr8_nanopolish2.fasta
bam=${root}/data/chm13_chr8_nanopolish2.bam
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/test

if [ "$1" == "meth" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -i ${root}/chr8.methylation.tsv |\
		sort -k1,1 -k2,2n | bgzip > ${outdir}/methylation.bed.gz
	tabix -p bed ${outdir}/methylation.bed.gz

    /home/gmoney/miniconda3/bin/python3 /home/isac/Code/nanopore-methylation-utilities/convert_bam_for_methylation.py --remove_poor --verbose -b $bam \
                -c ${root}/mcalls/methylation.bed.gz -w "chr8_v0.1_unpolished:1-10000" |\ samtools sort -o ${outdir}/test_meth.bam
samtools index  ${outdir}/test_meth.bam
fi

if [ "$1" == "bismark" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/methylation.bed.gz > ${outdir}/test.bismark
	bedtools genomecov -ibam ${outdir}/test_meth.bam -bg > ${outdir}/test_meth.bedgraph
fi
