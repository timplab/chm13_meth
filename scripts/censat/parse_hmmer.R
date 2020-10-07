library(rhmmer)
domtblout <- system.file('extdata', 'example.domtblout.txt', package='rhmmer')
tbl <- read_tblout("/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/DXZ4/DXZ4_alone_hmmsearch.tsv") %>%
  dplyr::rename("targ_name" = 1, "accession" = 2, "qry_name" = 3, "query_accession" =4, "hmmfrom" = 5, "hmmto" = 6, "alifrom" =7, "alito" = 8, "envfrom" = 9, "envto" = 10, "len" = 11, "strand" = 12) %>%
  select(alifrom, alito) %>%
  mutate(chr = "chrX") %>%
  mutate(start = alifrom+113996828) %>%
  mutate(end = alito +113996828) %>%
  mutate(len = end-start) %>%
  mutate(strand = case_when(
      len > 0 ~ "+",
      len < 0 ~ "-")) %>%
  select(chr, start, end, strand) 

ctcf <- tbl %>%
  mutate(ctcf_start = end-500) %>%
  mutate(ctcf_end = end-1000) %>%
  select(chr,ctcf_start, ctcf_end, strand)
  

write.table(ctcf, file = "/kyber/Data/Nanopore/Analysis/gmoney/CHM13/final_assembly/DXZ4/DXZ4_nhmmer_ctcf.tsv", append = FALSE, quote = F, sep = "\t", row.names = F,
            col.names = F)