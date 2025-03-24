#!/usr/bin/env Rscript

library(tidyverse)
library(vroom)
library("janitor")

impute_mean = function(X) {
  X_mean = colMeans(X,na.rm=T)
  X[is.na(X)] = X_mean[rep(1:ncol(X),colSums(is.na(X)))]
  X
}

admix<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/admixProp.GBSsamples.txt", col_names = c("ind","admix"))

####Test on subset####
#Coded up by Jeff Ross-Ibarra

#get vcf, make numeric, add allele freq column, transpose
untransformed<-vroom("/group/jrigrp11/snodgras_maizePopulations/admix_removal/test/jri_20KSNPs.numeric.tsv")
undata<-mutate(untransformed, Mex_AF = str_split(string = INFO, pattern = ";", simplify = T)[,1] %>%
              str_remove_all(.,"AF=")) %>% filter(Mex_AF != ".") %>% 
  mutate(Mex_AF=as.numeric(Mex_AF)) %>%
  select(Mex_AF, ID, starts_with("SEEDGWAS"))

tundata.names<-t(undata[,2:ncol(undata)])  %>% row_to_names(row_number = 1)

names <- rownames(tundata.names)
rownames(tundata.names) <- NULL
tundata.id.names <- data.frame(cbind(names,tundata.names))
colnames(tundata.id.names)[1]="name"
tundata.admix.id.names<-merge(tundata.id.names,admix,by.x="name",by.y="ind" )


tundata.id.names <- mutate(tundata.admix.id.names,name=str_split(name, "\\.",simplify = T)[,1])

dim(tundata.id.names) #check

# get metadata
meta<-vroom("/group/jrigrp11/snodgras_maizePopulations/clean_seedspassport.csv")

#merge untransformed, meta
mtundata<-merge(tundata.id.names,meta,by.x="name",by.y="Sample_ID_of_DNA_from_single_plants_used_in_GWAS")
forpca.mtundata=select(mtundata,2:ncol(tundata.id.names)) %>% select(-admix)
forpca.mtundata <- forpca.mtundata %>% mutate_if(is.character, as.numeric)


#PCA
#impute NA with means
tforpca.mtundata=impute_mean(as.matrix(forpca.mtundata))

#PCA
pca.mtundata <- prcomp(tforpca.mtundata,center=T,scale=F)
summary.pca<-summary(pca.mtundata)
summary.pca$importance[2,1:10] # % variance of first 10 PCs

#get eigenvecs
eigens<-pca.mtundata$x
new_eigens<-data.frame(cbind(eigens,mtundata$locations_elevation))

#plot
ggplot(new_eigens,aes(x=PC1,y=PC2,color=V2910))+geom_point(alpha=0.1)+
  xlab(paste("PC1 ",summary.pca$importance[2,1]*100, "%"))+
  ylab(paste("PC2 ",summary.pca$importance[2,2]*100, "%"))+
  theme_minimal()

cor(new_eigens$PC1,new_eigens$V2910)
cor(new_eigens$PC1,mtundata$locations_elevation)
cor(new_eigens$PC1,mtundata$admix)

summary(lm(new_eigens$PC1~new_eigens$V2910))
summary(lm(new_eigens$PC1~mtundata$admix))

#Runcie transform
genos<-as.matrix(forpca.mtundata)
genos_imputed<-impute_mean(genos)
n = nrow(genos_imputed)
a = cbind(1,mtundata$admix)
P = diag(1,n) - a %*% solve(crossprod(a),t(a))
X_resid = P %*% genos_imputed

Xpca<-prcomp(X_resid,center=T,scale=F)

summary.Xpca<-summary(Xpca)
summary.Xpca$importance[2,1:10] # % variance of first 10 PCs

#get eigenvecs
Xeigens<-Xpca$x

#plot
ggplot(Xeigens,aes(x=PC1,y=PC2,color=mtundata$admix))+geom_point(alpha=0.2)+
  xlab(paste("PC1 ",summary.Xpca$importance[2,1]*100, "%"))+
  ylab(paste("PC2 ",summary.Xpca$importance[2,2]*100, "%"))+
  theme_minimal()

cor(Xeigens[,1],Xeigens[,2])
cor(Xeigens[,1:10],mtundata$location_elevation)

summary(lm(new_eigens$PC1~new_eigens$V2910))
summary(lm(new_eigens$PC1~mtundata$admix))

####Run on full####
#merge of Sam's code with Jeff's

vcf<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf")

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

#transpose
undata<-select(vcf, Mex_AF, ID, starts_with("SEEDGWAS")) #non-transformed/original data

tundata.names<-t(undata[,2:ncol(undata)])  %>% row_to_names(row_number = 1)

names <- rownames(tundata.names)
rownames(tundata.names) <- NULL
tundata.id.names <- data.frame(cbind(names,tundata.names))
colnames(tundata.id.names)[1]="name"
tundata.admix.id.names<-merge(tundata.id.names,admix,by.x="name",by.y="ind" )


tundata.id.names <- mutate(tundata.admix.id.names,name=str_split(name, "\\.",simplify = T)[,1])

dim(tundata.id.names) #check

# add in the metadata to filter out samples that don't have geographical information

meta<-vroom("/group/jrigrp11/snodgras_maizePopulations/clean_seedspassport.csv")

#merge untransformed, meta
mtundata<-merge(tundata.id.names,meta,by.x="name",by.y="Sample_ID_of_DNA_from_single_plants_used_in_GWAS")
forpca.mtundata=select(mtundata,2:ncol(tundata.id.names)) %>% select(-admix)
forpca.mtundata <- forpca.mtundata %>% mutate_if(is.character, as.numeric)

#PCA
#impute NA with means
tforpca.mtundata=impute_mean(as.matrix(forpca.mtundata))

#PCA
pca.mtundata <- prcomp(tforpca.mtundata,center=T,scale=F)
summary.pca<-summary(pca.mtundata)

#write the pcs for the nontransformed data
pca.mtundata$x %>% rownames_to_column() %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/non_transformed_MAF0.01.pca")
summary.pca$importance[2,] %>% enframe() %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/non_transformed_MAF.01.eigens")

#transform
#Dan Runcie consulted on math by Jeff
#goal is to regress genotypes on admixture and then use the residuals as input for the PCA

genos<-as.matrix(forpca.mtundata)
genos_imputed<-impute_mean(genos)
n = nrow(genos_imputed)
a = cbind(1,mtundata$admix)
P = diag(1,n) - a %*% solve(crossprod(a),t(a))
X_resid = P %*% genos_imputed

Xpca<-prcomp(X_resid,center=T,scale=F)

summary.Xpca<-summary(Xpca)

#check correlation with admixture
cor(Xpca$x[,1:10], mtundata$admix)
cor(Xpca$x[,1:10], mtundata$location_elevation)

#write transformed PCs 
 
Xpca$x %>% rownames_to_column() %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/transformed_MAF0.01.pca")
summary.Xpca$importance[2,] %>% enframe() %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/transformed_MAF0.01.eigens")