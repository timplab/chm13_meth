#!/bin/bash 

[ -e /atium/Data/Nanopore/projects/210124_mason_bisulfite/tmp/TruSeq_HG002_LAB01_REP02_ILM_merge ]||mkdir -p /atium/Data/Nanopore/projects/210124_mason_bisulfite/tmp/TruSeq_HG002_LAB01_REP02_ILM_merge && /home/roham/software/bismark/Bismark-0.22.2/bismark --bam --non_directional --bowtie2 \
	-p 4 --genome /atium/Data/Nanopore/projects/210124_mason_bisulfite/reference 
        -1 /atium/Data/Nanopore/projects/210124_mason_bisulfite/TruSeq_HG002_LAB01_REP02_ILM_merge.R1.fastq.gz \
		-2 /atium/Data/Nanopore/projects/210124_mason_bisulfite/TruSeq_HG002_LAB01_REP02_ILM_merge.R2.fastq.gz \
                --temp_dir /atium/Data/Nanopore/projects/210124_mason_bisulfite/tmp/TruSeq_HG002_LAB01_REP02_ILM_merge \
                --output_dir /atium/Data/Nanopore/projects/210124_mason_bisulfite/bsseq/bismark &> /atium/Data/Nanopore/projects/210124_mason_bisulfite/bismark/TruSeq_HG002_LAB01_REP02_ILM_merge.align.log && touch /atium/Data/Nanopore/projects/210124_mason_bisulfite/bsseq/bismark/TruSeq_HG002_LAB01_REP02_ILM_merge_bismark_bt2_pe.bam
