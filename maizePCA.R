#!/usr/bin/env Rscript

args<-commandArgs(T)

library(tidyverse)

input<-read_tsv(args[1])

#remove all non-genotype columns
input<-input %>% select(contains("SEEDGWAS"))

#transpose so samples are rows and snps are columns
input<-t(input)

#impute NA with means across genotypes of a given sample (rowwise)
means<-apply(input, 1, mean,na.rm=TRUE)
for(i in 1:nrow(input)){
  input[i,is.na(input[i,])]=means[i]
}

#use if already mean-centered and subsetted samples
#pca<-as.matrix(input) %>% prcomp(center = F)

#use if samples haven't been subsetted
#sample_key<-read_tsv(args[4])
#pca<-as.matrix(input) %>% prcomp(center = F)

#use if NOT mean centered and subsetted
pca<- prcomp(input,center = T)

pca$x %>%
  as.data.frame() %>% 
  rownames_to_column() %>% 
  write_tsv(args[2])
  
  #correlation with missing data 
#pov of snps
#for(i in 1:10){
#  test<-cor.test(apply(test_snps,1,function(x) sum(is.na(x))), pca$x[,i])
#  print(paste("PC",i, "correlation estimate:",test$estimate," , correlation p-value:", test$p.value))
#}

#pca$sdev^2 %>% write_lines(args[3])
summary(pca)$importance[2,] %>% enframe() %>% write_tsv(args[3])
