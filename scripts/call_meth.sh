#!/bin/bash
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/ref/t2t-chm13.20200727.fasta
tsv=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/methylation_calls.tsv
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly
repo=/home/isac/Software/nanopolish-cpggpc
dir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly



if [ "$1" == "bed" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -q cpg -c 1.5 -i ${outdir}/methylation_calls.tsv -g $ref > ${outdir}/methylationCpG.tmp
	sort ${outdir}/methylationCpG.tmp -k1,1 -k2,2n | bgzip > ${outdir}/CpGmethylation.bed.gz
	tabix -p bed ${outdir}/CpGmethylation.bed.gz
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/CpGmethylation.bed.gz > ${outdir}/bismark.out
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

if [ "$1" == "bam" ]; then
	dir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/bam
	base=output
	samtools view -h -b -F 272 ${dir}/${base}.bam > ${dir}/${base}_filtered.bam
	samtools index ${dir}/${base}_filtered.bam
	/usr/bin/python /home/gmoney/repos/nanopore-methylation-utilities/convert_bam_for_methylation.py -t 65 --verbose -b ${dir}/${base}_filtered.bam \
		-c /kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/whole_genome/CpGmethylation.bed.gz -f $ref |\
		samtools sort -o ${dir}/${base}_filtered_meth.bam
	samtools index ${dir}/${base}_filtered_meth.bam
fi
