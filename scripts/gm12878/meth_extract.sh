#!/bin/bash
bam=/kyber/Data/Nanopore/Analysis/gmoney/alu/gm12878_jain/bsseq/bam/GM12878_BSseq_ENCLB794YYH_R1_val_1_bismark_bt2_pe.bam
dir=/kyber/Data/Nanopore/Analysis/gmoney/alu/gm12878_jain/bsseq/meth_extract

#samtools sort $bam -o ${dir}/GM12878_BSseq_ENCLB794YYH_R1_val_1_bismark_bt2_pe_sorted.bam
#samtools index ${dir}/GM12878_BSseq_ENCLB794YYH_R1_val_1_bismark_bt2_pe_sorted.bam
#samtools view -b ${dir}/GM12878_BSseq_ENCLB794YYH_R1_val_1_bismark_bt2_pe_sorted.bam "tig00000309_pilon_pilon:2358419-2439449" > ${dir}/GM12878_BSseq_brca1_bismark.bam
#samtools sort -n ${dir}/GM12878_BSseq_brca1_bismark.bam -o ${dir}/GM12878_BSseq_brca1_bismark_unsort.bam
#/home/isac/Code/miniconda3/bin/deduplicate_bismark --bam \
#	--output_dir $dir \
#	$bam 

/home/isac/Code/miniconda3/bin/bismark_methylation_extractor --gzip \
	--genome /kyber/Data/Nanopore/Analysis/gmoney/alu/gm12878_jain/bsseq/reference \
	--cytosine_report \
	--multicore 30 \
	--gazillion \
	-o ${dir} \
	${dir}/GM12878_BSseq_ENCLB794YYH_R1_val_1_bismark_bt2_pe.deduplicated.bam &> ${dir}/meth_extract.out 	
