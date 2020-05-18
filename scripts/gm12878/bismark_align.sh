#!/bin/bash
root=/kyber/Data/Nanopore/Analysis/gmoney/alu/gm12878_jain/bsseq
#fqdir=/dilithium/Data/NGS/projects/gm12878/bsseq/fqtrim
samp=GM12878_BSseq_ENCLB794YYH
dir=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/bsseq
ref=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/gm12878/bsseq/reference
bamdir=$root/bam
tmpdir=$dir/tmp
[ -e ${dir}/bam ]||mkdir ${dir}/bam
[ -e $tmpdir ]||mkdir $tmpdir
log=$dir/align.log

/home/isac/Code/miniconda3/bin/bismark_genome_preparation --verbose --path_to_bowtie /home/gmoney/miniconda3/bin/ $ref &> ${dir}/genome_prep.log

/home/isac/Code/miniconda3/bin/bismark --bam --gzip --non_directional --path_to_bowtie /home/gmoney/miniconda3/bin/ \
	--un --multicore 4 -p 3 --genome $ref \
	-1 $bamdir/$samp*val_1.fq.gz \
	-2 $bamdir/$samp*val_2.fq.gz \
	--temp_dir $tmpdir \
	--output_dir $bamdir &> $log
