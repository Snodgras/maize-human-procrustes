library(tidyverse)
library(vegan) #vegan 2.6-8

print(paste("Reading in files...",Sys.time()))
human_21km_pairs_PCs<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/human_21_pairs_PCs.tsv")

maize_wHumanPCcentroids<-human_21km_pairs_PCs %>% 
  select(starts_with("PC_"),sample_id) %>%
  group_by(sample_id) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)))

maize_genotypes<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/test.tsv") #read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv")
maize_genotypes<-select(maize_genotypes, ID, contains("SEEDGWAS"))

#transpose the genotypes so the SNPs are columns and individuals are rows and conserve names
print(paste("Transposing genotypes...", Sys.time()))
t_maize_genotypes<-t(select(maize_genotypes, contains("GWAS")))  %>% as_tibble(rownames = NA, .name_repair = TRUE) %>% rownames_to_column()
colnames(t_maize_genotypes)<-c("sample_id",maize_genotypes$ID)


#Modeled off of https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd

#Run PCA to get the structure?
pca<-rda(t_maize_genotypes[, -1], scale = TRUE)

#screeplot
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")

#need an object with the genotypes as a list column, and then each human PC and lat/long/elevation(?)

variables<-data.frame(maize_sample_ids = t_maize_genotypes[,1], maize_genotypes = t_maize_genotypes[,-1])
variables<-mutate(variables, sample_id = str_split(sample_id, ":",simplify=TRUE)[,1])
variables<-left_join(variables, maize_wHumanPCcentroids, by = c("sample_id"))
variables<-left_join(variables, select(human_21km_pairs_PCs, latitude.human, longitude.human, locations_elevation.maize, Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize), by = c("sample_id"="Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize"))

variables_samples<-select(variables, sample_id)
variables_genotypes<-select(variables, starts_with("maize_genotypes"))
variables_environments<-select(variables, -starts_with("maize_genotypes"), -sample_id)

variables<-data.frame(variables_samples, variables_genotypes, variables_environments)
#try rda
rda0<-rda(maize_genotypes ~ 1, variables)
