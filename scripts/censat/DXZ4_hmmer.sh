#!/bin/bash 
path=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/DXZ4
seq=${path}/HQ659113.1_DXZ4_formatted.fasta
# first cat together DXZ4 repeat from chm13 and the DXZ4 clone sequence from GenBank 
# align them as a msa with muscle then convert into a stockholm format for hmmer
#fasta_formatter -i HQ659113.1_clone4_DXZ4.fasta a -w 80 -o HQ659113.1_clone4_DXZ4_formatted.fasta format fastas first
#clustalo --in=$seq --out=${path}/msa_DXZ4_clustalo.sto -t DNA --force --outfmt=st --wrap=80 -v --threads 48
muscle -in $seq -out ${path}/msa_DXZ4_muscle_aln.fasta
hmmbuild msa_DXZ4_muscle_hmm.out msa_DXZ4_muscle_aln.fasta
hmmpress DXZ4_alone.hmm.out
nhmmscan --cpu 48 --tblout DXZ4_alone_hmmsearch.tsv DXZ4_alone.hmm.out DXZ4_formatted.fasta

# load output into R to parse with hmmer parser package
# now put together into multi fasta and do a msa
bedtools getfasta -fi ../ref/t2t-chm13.20200727.fasta -bed DXZ4_nhmmer_parsed.tsv -split -name > DXZ4_monomers.fasta
