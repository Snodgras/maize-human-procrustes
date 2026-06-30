#!/usr/bin/env Rscript

library(tidyverse)
library(vroom)
library(janitor,lib.loc="/home/snodgras/R/x86_64-conda-linux-gnu/r-4.3.3/")

impute_mean = function(X) {
  X_mean = colMeans(X,na.rm=T)
  X[is.na(X)] = X_mean[rep(1:ncol(X),colSums(is.na(X)))]
  X
}

admix<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/admixProp.GBSsamples.txt", col_names = c("ind","admix"))

####Run on full####
#merge of Sam's code with Jeff's

vcf<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv")

### split out mex AF from other INFO fields
vcf<-mutate(vcf, Mex_AF = str_split(string = INFO, pattern = ";", simplify = T)[,1] %>%
                            str_remove_all(.,"AF=")
         )
### Convert genotypes to numeric

vcf<-as_tibble(vcf) %>% mutate_at(vars(matches("SEEDGWAS")), as.numeric) 
if(!is.numeric(vcf$Mex_AF)){vcf<-mutate_at(vcf, vars("Mex_AF"), as.numeric)}

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
print("Merging genotypes with meta so that only doing PCA on samples with lat/long data")
mtundata<-merge(tundata.id.names,meta,by.x="name",by.y="Sample_ID_of_DNA_from_single_plants_used_in_GWAS")
forpca.mtundata=select(mtundata,2:ncol(tundata.id.names)) %>% select(-admix)
forpca.mtundata <- forpca.mtundata %>% mutate_if(is.character, as.numeric)

#PCA
#impute NA with means
print("imputing missing data as means")
tforpca.mtundata=impute_mean(as.matrix(forpca.mtundata))

#PCA
#print("Running non-transformed PCA")
#pca.mtundata <- prcomp(tforpca.mtundata,center=TRUE,scale=FALSE)
#summary.pca<-summary(pca.mtundata)

#write the pcs for the nontransformed data
#pcs.untransformed<-pca.mtundata$x
#pcs.untransformed %>% as.data.frame() %>% add_column(ind_id = mtundata[,1]) %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_non_transformed_LDprune_MAF0.01.pca")
#summary.pca$importance[2,] %>% enframe() %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_non_transformed_LDprune_MAF.01.eigens")

#transform
#Dan Runcie consulted on math by Jeff
#goal is to regress genotypes on admixture and then use the residuals as input for the PCA

print("Running transformation to remove admixture")
genos<-as.matrix(forpca.mtundata)
genos_imputed<-impute_mean(genos)
n = nrow(genos_imputed)
a = cbind(1,mtundata$admix)
P = diag(1,n) - a %*% solve(crossprod(a),t(a))
X_resid = P %*% genos_imputed

print("Writing transformed genotypes for LM analysis...")
write_tsv(as.data.frame(X_resid), "/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.genotypes.tsv")

print("Running transformed PCA")
Xpca<-prcomp(X_resid,center=TRUE,scale=FALSE)

summary.Xpca<-summary(Xpca)

#check correlation with admixture
print("Checking correlations with admixture and elevation across 1-10 PCs")
cor(Xpca$x[,1:10], mtundata$admix)
cor(Xpca$x[,1:10], mtundata$location_elevation)

#write transformed PCs 
 
Xpca$x %>% as.data.frame() %>% add_column(ind_id = mtundata[,1]) %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.pca")
summary.Xpca$importance[2,] %>% enframe() %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.eigens")
