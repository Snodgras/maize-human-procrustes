#!/usr/bin/env Rscript

library(tidyverse)

####Test using random data of right type, but no basis in biology####
#pilot_vcf<-tibble(SNP_ID = c(paste0("SNP",1:50)),
#                  Mex_AF = runif(n=50,min=0,max=1),
#                  IND_1 = round(runif(n=50, min=0,max=2)),
#                  IND_2 = round(runif(n=50, min=0,max=2)),
#                  IND_3 = round(runif(n=50, min=0,max=2)),
#                  IND_4 = round(runif(n=50, min=0,max=2)),
#                  IND_5 = round(runif(n=50, min=0,max=2)),
#                  IND_6 = round(runif(n=50, min=0,max=2)),
#                  IND_7 = round(runif(n=50, min=0,max=2)),
#                  IND_8 = round(runif(n=50, min=0,max=2)),
#                  IND_9 = round(runif(n=50, min=0,max=2)),
#                  IND_10 = round(runif(n=50, min=0,max=2)))
#pilot_admix<-tibble(Individual = c("IND_1","IND_2","IND_3","IND_4","IND_5","IND_6","IND_7","IND_8","IND_9","IND_10"),
#                    alpha = runif(n=10,min=0,max=1))
#
#(pilot_vcf[,3] - 2*pilot_vcf[,2]*pull(pilot_admix[which(pilot_admix[,1] == "IND_1"),2]))/(1-pull(pilot_admix[which(pilot_admix[,1] == "IND_1"),2]))
##test first entry:
#(1 - 2*0.172*0.00514)/(1-0.00514)
#
#pilot_vcf_transformed<-select(pilot_vcf, SNP_ID)
#for(i in pilot_admix$Individual){
#  alpha<-pilot_admix %>% filter(Individual == i) %>% select(alpha) %>% pull()
#  pilot_vcf_transformed<-add_column(pilot_vcf_transformed, !!i := pull((pilot_vcf[,i] - 2*pilot_vcf[,2]*alpha)/(1-alpha)))
#}
#remove(ls(pattern = "pilot"))

#### Test on a subset of the real data ####
#pilot_vcf<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/test.vcf")
#pilot_admix<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/admixProp.GBSsamples.txt", col_names = c("ind","admix"))

### split out mex AF from other INFO fields
#pilot_vcf<-mutate(pilot_vcf, Mex_AF = str_split(string = INFO, pattern = ";", simplify = T)[,1] %>%
#                            str_remove_all(.,"AF=")
#         )
### Convert genotypes to numeric
#pilot_vcf<-lapply(pilot_vcf, function(x) gsub("0/0","0",x))
#pilot_vcf<-lapply(pilot_vcf, function(x) gsub("0/1|1/0","1",x))
#pilot_vcf<-lapply(pilot_vcf, function(x) gsub("1/1","2",x))
#pilot_vcf<-lapply(pilot_vcf, function(x) gsub("./.",NA,x))
#
#pilot_vcf<-as_tibble(pilot_vcf) %>% mutate_at(vars(matches("SEEDGWAS")), as.numeric) 
#pilot_vcf<-mutate_at(pilot_vcf, vars("Mex_AF"), as.numeric)
### attempt the transformation
#
#### create the id variables
#pilot_transformed<-select(pilot_vcf, c("#CHROM","POS","ID","REF","ALT","Mex_AF"))
#### go through each individual (genotype column) and apply transformation
#for(i in pilot_admix$ind){
#  alpha<-pilot_admix %>% filter(ind == i) %>% select(admix) %>% pull()
#  pilot_transformed<-add_column(pilot_transformed, !!i := pull((pilot_vcf[,i] - 2*pilot_vcf[,"Mex_AF"]*alpha)/(1-alpha)))
#}
#
#### double checking that the transformation makes sense
#pca<- select(pilot_transformed, contains("SEEDGWAS")) %>% na.omit() %>% as.matrix() %>% prcomp(center = T)
#pilot_pcs<-pca$rotation[,1:2] 
#ggplot(pilot_pcs, aes(x=pilot_pcs[,1],y=pilot_pcs[,2]))+geom_point()+xlab("PC1")+ylab("PC2")
#summary(pilot_transformed[,1:20])
#pilot_transformed %>% summarize(across(contains("SEEDGWAS"), ~ mean(.x, na.rm =T))) %>% pivot_longer(cols = everything(),names_to = "ind",values_to = "mean")
#p<-pilot_transformed %>% summarize(across(contains("SEEDGWAS"), ~ min(.x, na.rm =T))) %>% 
#  pivot_longer(cols = everything(),names_to = "ind",values_to = "min") %>% 
#  ggplot()+geom_density(aes(x=min),color="dodgerblue", fill="dodgerblue", alpha = 0.75)
#p<-p+ geom_density( data = pilot_transformed %>% summarize(across(contains("SEEDGWAS"), ~ max(.x, na.rm =T))) %>% 
#  pivot_longer(cols = everything(),names_to = "ind",values_to = "max"),
#  aes(x=max), color = "red", fill="red", alpha =0.75)
#p+theme_bw()+xlab("Min and Max Genotype")+ylab("Density of individuals (4845)")

#####what about an individual snp
#snp_orig<-pilot_vcf[2,] %>% select(starts_with("SEEDGWAS"))
#snp_trans<-pilot_transformed %>% filter(ID == pull(pilot_vcf[2,3])) %>% select(starts_with("SEEDGWAS"))
#
#snp_orig<-snp_orig %>% pivot_longer(cols = everything(), values_to = "Orig_Genotype", names_to = "Ind")
#snp_trans<-snp_trans %>% pivot_longer(cols = everything(), values_to = "Transformed_Genotype", names_to = "Ind")

#inner_join(x=snp_orig, y=snp_trans, by="Ind") %>%
#  pivot_longer(cols = ends_with("Genotype"), names_to = "State", values_to = "Genotype") %>% 
#  ggplot(aes(x=State, y=Genotype, group = Ind))+
#  geom_point(aes(color = State))+
#  geom_line(color = "black", alpha = 0.1)+
#  theme_bw()+ggtitle("Single SNP genotype transformation")+guides(color = "none")

#remove(list = c("pilot_admix","pilot_pcs","pilot_transformed","pilot_vcf", "pca","p","snp_orig","snp_trans"))

####Final draft####
vcf<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf")
admix<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/admixProp.GBSsamples.txt", col_names = c("ind","admix"))

### split out mex AF from other INFO fields
vcf<-mutate(vcf, Mex_AF = str_split(string = INFO, pattern = ";", simplify = T)[,1] %>%
                            str_remove_all(.,"AF=")
         )
### Convert genotypes to numeric
vcf<-lapply(vcf, function(x) gsub("0/0","0",x))
vcf<-lapply(vcf, function(x) gsub("0/1|1/0","1",x))
vcf<-lapply(vcf, function(x) gsub("1/1","2",x))
vcf<-lapply(vcf, function(x) gsub("./.",NA,x))

vcf<-as_tibble(vcf) %>% mutate_at(vars(matches("SEEDGWAS")), as.numeric) 
vcf<-mutate_at(vcf, vars("Mex_AF"), as.numeric)
## attempt the transformation

### create the id variables
transformed<-select(vcf, c("#CHROM","POS","ID","REF","ALT","Mex_AF"))
### go through each individual (genotype column) and apply transformation
for(i in admix$ind){
  alpha<-admix %>% filter(ind == i) %>% select(admix) %>% pull()
  transformed<-add_column(transformed, !!i := pull((vcf[,i] - 2*vcf[,"Mex_AF"]*alpha)/(1-alpha)))
}

#compare to manual calculation
#gprime = function(g,pmex,alpha){(g-alpha*2*pmex)/(1-alpha)}


write_tsv(transformed, file = "/group/jrigrp11/snodgras_maizePopulations/admix_removal/transformed_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv")