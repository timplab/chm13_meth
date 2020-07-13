#!/bin/bash 
fastq=/uru/Data/Nanopore/Analysis/gmoney/chm13/hg002/cenX_reads_hg002_lt100kb.fastq
fast5=/uru/Data/Nanopore/Analysis/gmoney/chm13/hg002/filteredFiles
bam=/uru/Data/Nanopore/Analysis/gmoney/chm13/hg002/cenX_reads_hg002_lt100kb.asm8_primary.tig00021361.sorted.bam
out=/uru/Data/Nanopore/Analysis/gmoney/chm13/hg002/nanopolish
repo=/home/isac/Software/nanopolish
ref=/uru/Data/Nanopore/Analysis/gmoney/chm13/hg002/asm8_primary.tig00021361.fa
base=asm8_primary.tig00021361

if [ "$1" == "nanopolish" ]; then

        ${repo}/nanopolish index -d ${fast5} ${fastq}
	${repo}/nanopolish call-methylation -t 72 --reads ${fastq} --bam ${bam} --genome ${ref} > ${out}/${base}_methylation_calls.tsv
        ${repo}/scripts/calculate_methylation_frequency.py ${out}/${base}_methylation_calls.tsv > ${out}/${base}_methylation_frequency.tsv


        python3 /home/gmoney/repos/nanopore-methylation-utilities/mtsv2bedGraph_upperlower.py -c 1.5  -i ${out}/${base}_methylation_calls.tsv -g $ref > ${out}/${base}_methylation.tmp
        cat ${out}/${base}_methylation.tmp | sort -k1,1 -k2,2n | bgzip > ${out}/${base}_methylation.bed.gz
        tabix -p bed ${out}/${base}_methylation.bed.gz
	rm ${out}/${base}_methylation.tmp
	
	python3 /home/gmoney/repos/nanopore-methylation-utilities/parseMethylbed.py frequency -i ${out}/${base}_methylation.bed.gz > ${out}/${base}_bismark.out
        python3 /home/gmoney/repos/nanopore-methylation-utilities/convert_bam_for_methylation.py -t 48 --remove_poor --verbose -b ${bam} \
		-c ${out}/${base}_methylation.bed.gz -f ${ref} |\
                samtools sort -o ${out}/${base}_meth.bam
        samtools index ${out}/${base}_meth.bam
	${repo}/scripts/calculate_methylation_frequency.py ${out}/${base}_methylation_calls.tsv > ${out}/${base}_methylation_frequency.tsv


fi
