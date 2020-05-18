#!/bin/bash
ref=/uru/Data/Nanopore/Analysis/gmoney/chm13/draft_v0.7/chrX_v0.7.nanopolish2.arrow2_10x.fasta
bam=/uru/Data/Nanopore/Analysis/gmoney/chm13/draft_v0.7/chrX.bam
outdir=/uru/Data/Nanopore/Analysis/gmoney/chm13/draft_v0.7/meth_call
indir=/uru/Data/Nanopore/Analysis/gmoney/chm13/draft_v0.7

if [ "$1" == "meth" ]; then
        python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph.py -i ${indir}/methylation_calls.tsv -g $ref |\
                sort -k1,1 -k2,2n | bgzip > ${outdir}/methylation.bed.gz
        tabix -p bed ${outdir}/methylation.bed.gz

fi

if [ "$1" == "bam" ]; then
    /home/gmoney/miniconda3/bin/python3 /home/isac/Code/nanopore-methylation-utilities/convert_bam_for_methylation.py --remove_poor --verbose -b $bam \
                -c ${outdir}/methylation.bed.gz -f $ref |\
samtools sort -o ${outdir}/chm13_chrX_meth.bam
samtools index  ${outdir}/chm13_chrX_meth.bam
fi

if [ "$1" == "bismark" ]; then
        python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${outdir}/methylation.bed.gz > ${outdir}/bismark.out

fi

if [ "$1" == "bedgraph" ]; then
#       bedtools genomecov -ibam $bam -bg > ${outdir}/chm13_v0.6.nanopolish2.arrow2_10x2.chrX.bg
        bedtools genomecov -ibam ${outdir}/chm13_chrX_meth.bam -bg > ${outdir}/chm13_meth.chrX.bg
fi
