#!/bin/bash


ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/ref/chrX_fixedBionanoSV_centromereV4_racon_patch139_arrow_arrow.fasta
bam=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/bam/subset.bam
outdir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v6_assembly/meth_call

if [ "$1" == "meth" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -i ${outdir}/chrX.methylation.tsv |\
		sort -k1,1 -k2,2n | bgzip > methylation.bed.gz
	tabix -p bed methylation.bed.gz

fi

if [ "$1" == "bam" ]; then
    python3.7 /home/gmoney/repos/nanopore-methylation-utilities/convert_bam_for_methylation.py --verbose -b $bam \
                -c ${outdir}/methylation.bed.gz -f $ref > ${outdir}/test.bam #|\
               # samtools sort -o ${outdir}/test.bam
   # samtools index  ${outdir}/test.bam
fi

if [ "$1" == "bismark" ]; then
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/methylation.bed.gz > ${outdir}/bismark.out

fi

if [ "$1" == "bedgraph" ]; then
	bedtools genomecov -ibam $bam -bg > ${outdir}/chm13_v0.6.nanopolish2.arrow2_10x2.chrX.bg
fi
