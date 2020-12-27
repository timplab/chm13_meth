
library(tidyverse)
dat="/kyber/Data/Nanopore/Analysis/gmoney/CHM13/v1.0_final_assembly"

censat = read_tsv(paste0(dat, "/annotations/t2t-chm13.v1.0.cenSat_annotationTab.bed"), col_names = F) %>%
  mutate(name = ifelse(grepl("HSat1", X4), "HSAT1", X4)) %>%
  mutate(name = ifelse(grepl("HSat2", X4), "HSAT2", name)) %>%
  mutate(name = ifelse(grepl("HSat3", X4), "HSAT3", name)) %>%
  mutate(name = ifelse(grepl("HSat4", X4), "HSAT4", name)) %>%
  mutate(name = ifelse(grepl("HSat5", X4), "HSAT5", name)) %>%
  mutate(name = ifelse(grepl("HSat6", X4), "HSAT6", name)) %>%
  mutate(name = ifelse(grepl("ct", X4), "CT", name)) %>%
  mutate(name = ifelse(grepl("bsat", X4), "BSAT", name)) %>%
  mutate(name = ifelse(grepl("hor", X4), "HOR", name)) %>%
  mutate(name = ifelse(grepl("mon", X4), "MON", name)) %>%
  mutate(name = ifelse(grepl("Alu", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("SATR", X4), "SATR", name)) %>%
  mutate(name = ifelse(grepl("ACRO", X4), "ACRO", name)) %>%
  mutate(name = ifelse(grepl("GSATII", X4), "GSATII", name)) %>%
  mutate(name = ifelse(grepl("TAR", X4), "TAR", name)) %>%
  mutate(name = ifelse(grepl("TE", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("MER", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("MST", X4), "MST", name)) %>%
  mutate(name = ifelse(grepl("CER", X4), "CER", name)) %>%
  mutate(name = ifelse(grepl("L1", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("SST", X4), "SST", name)) %>%
  mutate(name = ifelse(grepl("LSAU", X4), "LSAU", name)) %>%
  mutate(name = ifelse(grepl("GSAT", X4), "GSAT", name)) %>%
  mutate(name = ifelse(grepl("MSAT", X4), "MSAT", name)) %>%
  mutate(name = ifelse(grepl("novel", X4), "novel", name)) %>%
  mutate(name = ifelse(grepl("HERV", X4), "TE", name)) %>%
  mutate(name = ifelse(grepl("LT", X4), "TE", name)) %>%
  dplyr::select(c(X1, X2, X3, name)) %>%
  dplyr::rename("chrom" =1, "start" = 2 ,"end" =3)

write.table(censat, file = paste0(dat, "/annotations/t2t-chm13.v1.0.cenSat_annotationFormatted.bed"), append = FALSE, quote = F, sep = "\t", row.names = F,
            col.names = F)