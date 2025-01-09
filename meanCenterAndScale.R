#!/usr/bin/env Rscript

args<-commandArgs(T)
tag<-args[2] #this will be the path to the directory the output file will be written

library(tidyverse)

#vcf without the header data, but the column names still included
#grep -v "^##" full.vcf > full.headerless.vcf #run in terminal
input_vcf<-read_tsv(args[1],na="[.]")

output_vcf<-input_vcf %>% mutate(genotype_mean = rowMeans(.[,10:length(input_vcf)]))
output_vcf$genotype_sd<-apply(output_vcf[,10:length(input_vcf)], 1, sd)

output_vcf<-output_vcf %>% mutate(across(matches("SEEDGWAS"),~(genotype_mean - .x)/genotype_sd))

write_tsv(select(output_vcf,-c(genotype_mean, genotype_sd)), file = paste0(tag,args[3],"centeredMaizeGBS.vcf"))
