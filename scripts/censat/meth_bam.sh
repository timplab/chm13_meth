#!/bin/bash
ref=/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_1/meth_calls/reference/chm13_20k_hicanu_hifi.fasta
bam=/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2
outdir=/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2


if [ "$1" == "meth" ]; then
	dir="$(ls -d ${bam}/tig*)"
	echo $dir
	echo $dir | xargs -n 1 | parallel -P 10 "python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -i {}/methylation_calls.tsv -g $ref > {.}/methylation.tmp"
	echo $dir | xargs -n 1 | parallel -P 10	"sort {}/methylation.tmp -k1,1 -k2,2n | bgzip > {.}/methylation.bed.gz"
	echo $dir |  xargs -n 1 | parallel -P 10 "tabix -p bed {}/methylation.bed.gz"
fi


if [ "$1" == "bam" ]; then
		dir=
		base="$(basename $dir/*bam .bam)"
		echo $base
		echo ${dir}/${base}.bam
		/home/gmoney/miniconda3/bin/python3 /home/isac/Code/nanopore-methylation-utilities/convert_bam_for_methylation.py -t 35 --remove_poor --verbose -b ${dir}/${base}.bam \
			-c ${dir}/methylation.bed.gz -f $ref |\
			samtools sort -o ${dir}/${base}_meth.bam
		samtools index ${dir}/${base}_meth.bam
fi


if [ "$1" == "bismark" ]; then

	dir="$(ls -d ${bam}/tig*)"
	echo $dir
	echo $dir | xargs -n 1 | parallel -P 10 "python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i {}/methylation.bed.gz > {.}/bismark.out"
fi

if [ "$1" == "bedgraph" ]; then
	dir="$(ls -d ${bam}/tig*/*.bam)"
	echo $dir
	echo $dir | xargs -n 1 | parallel -P 10 "bedtools genomecov -ibam {} -bg > {.}_cov.bg"
fi
