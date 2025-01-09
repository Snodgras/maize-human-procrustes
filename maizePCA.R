#!/usr/bin/env Rscript

args<-commandArgs(T)

library(tidyverse)

input<-read_tsv(args[1])

#use if already mean-centered and subsetted samples
#pca<-input %>% select(contains("SEEDGWAS")) %>% as.matrix() %>% prcomp(center = F)

#use if samples haven't been subsetted
sample_key<-read_tsv(args[4])
pca<-input %>% select(contains(sample_key$vcf_sample_id)) %>% as.matrix() %>% prcomp(center = F)


pca$rotation[,1:20] %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  write_tsv(args[2])

pca$sdev^2 %>% write_lines(args[3])
