    
#!/usr/bin/env Rscript 

# supply /path/to/sumdir chrx:xxx-xxx 
# optional smoothed wide

#parse arguments 

library("optparse")

option_list = list(
   make_option(c("-d", "--directory"), type="character", default=NULL, 
               help="path to directory with bedgraph, methylation tsv, methylation frequency, bismark files"),
   make_option(c("-s", "--smoothed"), type="integer", default=500, 
               help="number of nucleotides for bsseq smoothing", metavar="number"), 
   make_option(c("-w", "--wide"), type="integer", default=500000, 
               help="size of flanking regions", metavar="number"), 
   make_option(c("-p", "--png"), type="integer", default=NULL, 
               help="size of flanking regions", metavar="number"),
   make_option(c("-c", "--coordinates"), type="character", default=NULL, 
               help="chrx:xx:xx", metavar="character"),
   make_option(c("-o", "--ouput"), type="character", default=getwd(), 
               help="output directory", metavar="character"),
   make_option(c("-r", "--roi"), type="character", default=NULL, 
               help="region in the png", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$directory)

if (is.null(opt$directory)){
   print_help(opt_parser)
   stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

if (is.null(opt$coordinates)){
   print_help(opt_parser)
   stop("At least one argument must be supplied (input coordinates).\n", call.=FALSE)
}

if (is.null(opt$png)){
   pic = FALSE
}else{
         pic = TRUE
   }
   


#args = commandArgs(trailingOnly=TRUE)
#
#if (length(args)<2) {
#   stop("provide path to summary files and coordinates for region", call.=FALSE)
#} else if (length(args)==2) {
#   # default output file
#   args[3] = 500
#   args[4] = 500000
#}

#print(args)
library(knitr)
library(markdown)
library(rmarkdown)
library(tidyverse)

setwd(opt$directory)

file_list <- list()

file_list[1:2] <- list.files(path=opt$directory, pattern=".tsv")
file_list[3] <- list.files(path=opt$directory, pattern="bismark")
file_list[4] <- list.files(path=opt$directory, pattern="bedgraph")
roi <- opt$coordinates
smoothed <- as.numeric(opt$smoothed)
wide <- as.numeric(opt$wide)

summary <- list()
print(file_list)
title="Methylation Report"
files = list(summary = paste(opt$directory, "/", file_list, sep = ''))
name=substr(file_list[2],1,nchar(file_list[2])-4)

rmarkdown::render(input = '/home/gmoney/Code/chm13/chm13_meth/scripts/summary_report/summary.rmd',
output_format = "pdf_document",
output_file = paste(name, ".pdf", sep = ''),
output_dir = opt$output)

