library(vroom)
library("janitor")

impute_mean = function(X) {
  X_mean = colMeans(X,na.rm=T)
  X[is.na(X)] = X_mean[rep(1:ncol(X),colSums(is.na(X)))]
  X
}


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







#do transformation
transformed=data.frame(undata$ID)

admix<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/admixProp.GBSsamples.txt", col_names = c("ind","admix"))


for(i in admix$ind){
  alpha<-filter(admix, ind == i) %>% 
    select(admix) %>% 
    pull()
  transformed<-add_column(transformed, !!i := pull((undata[,i] - 2*undata[,"Mex_AF"]*alpha)/(1-alpha)))
}

new.transformed=transformed
colnames(new.transformed)=str_split(colnames(transformed), "\\.",simplify = T)[,1]
forpca.tranfsormed=select(new.transformed,mtundata$name)

tforpca.transformed=t(forpca.tranfsormed)

#PCA
#impute NA with means
#forpca.mtundata <- forpca.mtundata %>% mutate_if(is.character, as.numeric)

means<-apply(tforpca.transformed,1, mean,na.rm=TRUE)
forpca.mtdata=impute_mean(as.matrix(tforpca.transformed))

#PCA
pca.mtdata <- prcomp(forpca.mtdata,center=T,scale=F)
summary.pcamt<-summary(pca.mtdata)
summary.pcamt$importance[2,1:10] # % variance of first 10 PCs

#get eigenvecs
teigens<-pca.mtdata$x
new_teigens<-data.frame(cbind(teigens,mtundata$locations_elevation))

#plot
ggplot(eigens,aes(x=PC1,y=PC2))+geom_point(alpha=0.1)+
  xlab(paste("PC1 ",summary.pca$importance[2,1]*100, "%"))+
  ylab(paste("PC2 ",summary.pca$importance[2,2]*100, "%"))+
  theme_minimal()




# test to check transformation
testdude<-sample(colnames(transformed),1)
testsnp<-sample(transformed$undata.ID,1)
bob<-filter(undata,ID %in% testsnp) %>% 
  select(testdude,Mex_AF)
testmix<-filter(admix,ind %in% testdude)




filter(transformed,undata.ID %in% testsnp) %>% 
  select(testdude)







pca<-read.table("~/jritrysam.eigenvec",header=F)
eigenval<-read.table("~/jritrysam.eigenval",header=F)
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
plot(pca$PC2~pca$PC1)
cor(pca$PC1,pca$PC2)


transformed<-vroom("/group/jrigrp11/snodgras_maizePopulations/admix_removal/test/20000RandomeSNPs_transformed_withHeader.tsv")
transformed<-vroom("/group/jrigrp11/snodgras_maizePopulations/admix_removal/test/test.transformed.tsv")
transformed.imputed<-transformed[,7:ncol(transformed)]
transformed.imputed<-t(transformed.imputed)

#impute NA with means
means<-apply(transformed.imputed,1, mean,na.rm=TRUE)
for(i in 1:nrow(transformed.imputed)){
  transformed.imputed[i,is.na(transformed.imputed[i,])]=means[i]
}

#PCA
transformed.pca <- prcomp(transformed.imputed,center=T,scale=T)
summary.pca<-summary(transformed.pca)
summary.pca$importance[2,1:10] # % variance of first 10 PCs

#get eigenvecs
eigens<-transformed.pca$x
#no correlation of first two PCs
cor(eigens[,1],eigens[,2])

#plot
ggplot(eigens,aes(x=PC1,y=PC2))+geom_point(alpha=0.1,color="blue")+
  xlab(paste("PC1 ",summary.pca$importance[2,1]*100, "%"))+
  ylab(paste("PC2 ",summary.pca$importance[2,2]*100, "%"))+
  theme_minimal()

#missigness check
missing<-sapply(7:ncol(transformed),function(i) sum(is.na(transformed[7:ncol(transformed),i])))
cor(missing,eigens[,1:10])
  
  
