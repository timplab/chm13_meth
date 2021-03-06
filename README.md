# chm13_meth
Methylation analysis of centromere vs whole chromosome in chm13 assembly

# dependencies
knitr
tidyverse
bsseq
Biostrings
GenomicRanges
Rsamtools
Sushi
nanopore-methylation-utilities

# generating input files
Call methylation for chromsome
```
nanopolish call-methylation -t 8 -r albacore_output.fastq -b albacore_output.sorted.bam -g reference.fasta -w "chr8" > chr8_methylation_calls.tsv

scripts/calculate_methylation_frequency.py -i methylation_calls.tsv > chr8_methylation_frequency.tsv

```
use Isac's repo nanopore-methylation-utilities to generate bam, bed, bismark and bedtools to generate bedgraph

```
python3 mtsv2bedGraph.py -i chr8.methylation.tsv |\
                sort -k1,1 -k2,2n | bgzip > ${outdir}/methylation.bed.gz
tabix -p bed methylation.bed.gz

convert_bam_for_methylation.py --remove_poor --verbose -b bam \
                -c mcalls/methylation.bed.gz -f ref.fasta |\ 
                samtools sort -o test_meth.bam
samtools index test_meth.bam

 python3 parseMethylbed.py frequency -i methylation.bed.gz > test.bismark
 bedtools genomecov -ibam test_meth.bam -bg > test_meth.bedgraph
```
keep all these files in one directory

# generating report

```
Rscript chm13_meth/summary_report/call_summary.R -d /path/to/files -c chr8:1000000-2000000
```
compare centromere methylation to entire chromosome with -c [cetromere region]
```
Usage: call_summary.R [options]


Options:
        -d DIRECTORY, --directory=DIRECTORY
                path to directory with bedgraph, methylation tsv, methylation frequency, bismark files

        -s NUMBER, --smoothed=NUMBER
                number of nucleotides for bsseq smoothing

        -w NUMBER, --wide=NUMBER
                size of flanking regions

        -p NUMBER, --png=NUMBER
                path to igv png images

        -c CHARACTER, --coordinates=CHARACTER
                chrx:xx-xx

        -o CHARACTER, --ouput=CHARACTER
                output directory

        -r CHARACTER, --roi=CHARACTER
                region in the png

        -h, --help
                Show this help message and exit


```
![example_report](https://github.com/gmoneyomics/chm13_meth/tree/master/plots/chr8.methylation.pdf)
