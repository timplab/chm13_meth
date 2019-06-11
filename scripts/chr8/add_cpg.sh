#!/bin/bash
root=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/chr8/data
ref=${root}/chm13_chr8_nanopolish2.fasta
bam=${root}/chm13_chr8_nanopolish2.bam
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/chr8/mcalls

if [ "$1" == "meth" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -i ${root}/chr8.methylation.tsv |\
		sort -k1,1 -k2,2n | bgzip > ${outdir}/methylation.bed.gz
	tabix -p bed ${outdir}/methylation.bed.gz

    /home/gmoney/miniconda3/bin/python3 /home/isac/Code/nanopore-methylation-utilities/convert_bam_for_methylation.py --remove_poor --verbose -b $bam \
                -c ${outdir}/methylation.bed.gz -f $ref |\
samtools sort -o ${outdir}/chm13_chr8_meth.bam
samtools index  ${outdir}/chm13_chr8_meth.bam
fi

if [ "$1" == "bismark" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/methylation.bed.gz > ${outdir}/chm13_chr8.bismark
	bedtools genomecov -ibam ${outdir}/chm13_chr8_meth.bam -bg > ${outdir}/chm13_ch8_meth.bedgraph
fi
