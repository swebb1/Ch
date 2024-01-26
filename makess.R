library(tidyverse)
library(readxl)
t<-read_xlsx("ChIP-seq_mm_samples.xlsx",col_names = T) %>% 
  write_tsv(file = "samples.tsv",col_names = T)
