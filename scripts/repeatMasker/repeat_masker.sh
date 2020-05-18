#!/bin/bash

# FLAG descrip : 

  #   -s :   use the 'slow' search option --  ( more sensitive )
  #   -species :  specifies species (or clade) of input seq 
  #   -a (alignment) :  writes alignment in .align output file 
  #   -small  :   returns the WHOLE masked sequence - in lower case  (versus .. ithink. .. replacing repeats with 'xxxxxx' )  
  #   -u  :  Creates an additional annotation file not processed by ProcessRepeats
  #   -gff  :  Creates an additional Gene Feature Finding format output
  #   -xm   :  Creates an additional output file in cross_match format (for parsing)
  #   -alu  :   Only masks Alus (and 7SLRNA, SVA and LTR5)(only for primate DNA)
 
 # [ note : default output is in the query directory  ] 


# flip one assembly (reverse compliment)
#bedtools getfasta -fi /mithril/Data/NGS/Reference/human38/GRCH38.fa  -bed ${maindir}/hg38/hg38_alu3.bed > ${maindir}/hg38/alu3.fasta
#seqtk seq -r $maindir/hap1/*.fasta > $maindir/hap1_flip/hap1_flip.fasta

bam=/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2
ref=/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_1/meth_calls/reference
for dir in ${bam}/tig*
do
	echo $dir
        base="$(basename $dir)"
        echo $base
	samtools faidx $ref/chm13_20k_hicanu_hifi.fasta $base > $dir/$base.fasta
	
	input_fa=$dir/$base.fasta
	echo 'this is input fasta:'
	echo $input_fa
	
	bioawk -c fastx '{print $name"\t"length($seq)}' $input_fa > $dir/contig_sizes.tsv
	echo 'startinn upp with that mutha f@@ckin repeat a-masskkkkin'
	
	
	/usr/local/RepeatMasker/RepeatMasker -parallel 50 -s -species human $input_fa
	
	RM_output=$dir/*fasta.out
	
	echo 'this is repeat masker output file that we use'
	echo $RM_output
	
	new_tsv=$dir/${base}_forPlot.tsv
	
	echo 'reformatting repeat masker output into a file that is actually useable' 
	
	python ./helper_scripts/format_RM_output.py -i $RM_output -o  $new_tsv -a $base

done 

