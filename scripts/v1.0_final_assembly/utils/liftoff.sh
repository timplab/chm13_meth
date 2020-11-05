#!/bin/bash 
ref=/atium/Data/old_mithril_ref/human38/GCF_000001405.39_GRCh38.p13_genomic.fna
asm=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta
gtf=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/liftoff/GCF_000001405.39_GRCh38.p13_genomic.gff
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/annotations
if [ "$1" == "liftoff" ]; then
	liftoff -t $asm -r $ref -g $gtf -o ${out}/hg38_genes_t2t_liftoff.gtf -p 48 -m ~/repos/minimap2/minimap2
fi
