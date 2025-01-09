#Corrlations of maize with human PCs
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)

####Read in files ####
maize_vcf<-read_tsv(args[1], na = c("","NA",".")) #21km_maize_centeredMaizeGBS.vcf

human_paired_pcs<-read_tsv(args[2]) #maybe change the column names? 

####Join files####

#test.maize<-tibble(SNP_ID=c("snp1","snp2","snp3","snp4","snp5"),
#                   sample1=c(1,0,1,1,1),
#                   sample2=c(0,1,1,0,0),
#                   sample3=c(0,0,1,".",0),
#                   sample4=c(1,0,0,0,1),
#                   sample5=c(1,1,0,0,0))
#test.human<-tibble(maize_sample=c("sample1","sample2","sample3","sample4","sample5"),
#                   human_PC1=c(84.80282,113.19731,118.36005,80.00427,75.96019),
#                   human_PC2=c(-216.9386,-246.6669,-252.8669,-178.5138,-182.0036))

#change "." into NAs in vcf
#better to specify NA strings when reading in the files
#test.maize$sample3<-test.maize$sample3 %>% as.numeric()

####create correlations dataframe####

#this doesn't assume the dfs start in the same order wrt maize samples

#create loop
test.final<-tibble(SNP_ID=NA, PC1.cor=NA, PC2.cor=NA)
for(i in 1:nrow(maize_vcf)){ #must change the "cols=starts_with" argument to match vcf
  test.joined<-maize_vcf[i,] %>% pivot_longer(cols=starts_with("SEEDGWAS"), names_to = "maize_sample", values_to = "genotype") %>% 
    inner_join(.,human_paired_pcs, by="maize_sample")
  test.final<-add_row(test.final,SNP_ID=test.joined$SNP_ID[1],
                      PC1.cor=cor(test.joined$genotype,test.joined$human_PC1,use = "complete.obs"),
                      PC2.cor=cor(test.joined$genotype,test.joined$human_PC2,use = "complete.obs"))
}

#write the output to text file
test.final %>% na.omit() %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/PCcorrelations.tsv")
