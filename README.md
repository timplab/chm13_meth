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

# generating report
Rscript --vanilla call_summary.R [path/to/dir/np_meth_utilities_outputs] [chrx:xxxx-xxxxx]
