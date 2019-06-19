# chm13_meth
Methylation analysis of chrx in chm13 assembly

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

use Isac's repo nanopore-methylation-utilities to generate bam, bed, bismark and bedtools to generate bedgraph

python3 mtsv2bedGraph.py -i ${root}/chr8.methylation.tsv |\
                sort -k1,1 -k2,2n | bgzip > ${outdir}/methylation.bed.gz
tabix -p bed methylation.bed.gz

convert_bam_for_methylation.py --remove_poor --verbose -b bam \
                -c ${root}/mcalls/methylation.bed.gz -w "chr8:1-10000" |\ samtools sort -o test_meth.bam
samtools index test_meth.bam

 python3 parseMethylbed.py frequency -i ${outdir}/methylation.bed.gz > test.bismark
 bedtools genomecov -ibam ${outdir}/test_meth.bam -bg > test_meth.bedgraph

keep all these files in one directory

# generating report

Rscript call_summary.R -d /path/to/files -c chr8:10-20

use whole chromosome as input -c is region of centromere
Usage: call_summary.R [options]


Options:
        -d DIRECTORY, --directory=DIRECTORY
                path to directory with bedgraph, methylation tsv, methylation frequency, bismark files

        -s NUMBER, --smoothed=NUMBER
                number of nucleotides for bsseq smoothing

        -w NUMBER, --wide=NUMBER
                size of flanking regions

        -p NUMBER, --png=NUMBER
                size of flanking regions

        -c CHARACTER, --coordinates=CHARACTER
                chrx:xx:xx

        -o CHARACTER, --ouput=CHARACTER
                output directory

        -r CHARACTER, --roi=CHARACTER
                region in the png

        -h, --help
                Show this help message and exit
