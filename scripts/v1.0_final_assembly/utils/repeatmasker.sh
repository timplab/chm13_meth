#!/bin/bash 
input_fa=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/reference/chm13.draft_v1.0.fasta
out=/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly/annotations
/usr/local/RepeatMasker/RepeatMasker -parallel 68 -dir ${out} -s -species human $input_fa
