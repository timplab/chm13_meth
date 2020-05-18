
#!/usr/bin/env Rscript

library(knitr)
library(markdown)
library(rmarkdown)
suppressMessages(library(tidyverse))
library(optparse)

option_list = list(
  make_option(c("-d", "--directory"), type = "character", default = NULL, 
              help = "top level directory where contig directories located"),
  make_option(c("-c", "--coords"), type="character", default=NULL,
              help="tsv file with coordinates"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="path to output directory"),
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$directory)){
  print_help(opt_parser)
  stop("provide directory.\n", call.=FALSE)
}



if (is.null(opt$coords)){
  print_help(opt_parser)
  stop("provide file with coordinates of interest.\n", call.=FALSE)
}


if (is.null(opt$out)){
  print_help(opt_parser)

    stop("Specify an output directory.\n", call.=FALSE)
  }
  if (is.null(opt$markdown)){
    print_help(opt_parser)
    stop("provide path to rmarkdown.\n", call.=FALSE)
  }
  

opt = list()
opt$directory="/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2"
opt$coords="/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2/censat.tsv"
opt$out="/uru/Data/Nanopore/Analysis/gmoney/chm13/censat/phase_2/figures"

coords= read_tsv(opt$coords)

  for (i in 1: length(coords$hicanu_ctg)){
    print(coords$hicanu_ctg[i])
    title = coords$hicanu_ctg[i]
    tig = coords$hicanu_ctg[i]
    start = coords$ctg_s[i]
    end = coords$ctg_e[i]
    title = coords$name[i]
    
    params = list(tig = tig, 
                  bismark = paste(opt$directory,"/", tig, "/bismark.out",sep = ''), 
                  cov = paste(opt$directory,"/", tig, "/", tig, "_meth_cov.bg",sep = ''), 
                  repeats = paste(opt$directory,"/", tig, "/", tig, "_forPlot.tsv",sep = ''),
                  newtitle = title, 
                  cen1 = start,
                  cen2 = end)
    
    rmarkdown::render("/home/gmoney/Code/chm13/chm13_meth/scripts/censat/censat_phase2.rmd",
                      output_format = "pdf_document",
                      output_file = paste(title, ".pdf", sep = ''),
                      output_dir = opt$out)
  }
  