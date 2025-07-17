#!/usr/bin/env Rscript

library(tidyverse)

#args<-commandArgs(T)

#human_21km_pairs_PCs object from plottingPCbyGeography.R which is only Mexican samples used to calculate PCs

#nrow(human_21km_pairs_PCs) == length(unique(human_21km_pairs_PCs$sample_id))

#maize_per_human<-human_21km_pairs_PCs %>% select(sample_id.human) %>% group_by(sample_id.human) %>%count()
#human_per_maize<-human_21km_pairs_PCs %>% select(sample_id) %>% group_by(sample_id) %>% count()

#left_join(human_21km_pairs_PCs, maize_per_human, by="sample_id.human") %>% 
#  select(longitude.human, latitude.human, Region.human, n) %>% 
#  ggplot(aes(x=longitude.human, y=latitude.human, color = Region.human))+
#  geom_jitter(width = 0.1, height = 0.1, aes(size=n), alpha = 0.5)+
#  theme_bw() + ggtitle("Maize per Human")

#left_join(human_21km_pairs_PCs, human_per_maize, by="sample_id") %>% 
#  select(locations_longitude.maize, locations_latitude.maize, n, Region.human) %>% 
#  ggplot(aes(x=locations_longitude.maize, y=locations_latitude.maize, color = Region.human ))+
#  geom_jitter(width = 0.1, height = 0.1, aes(size=n), alpha = 0.5)+
#  theme_bw() + ggtitle("Humans per Maize")
print(paste("Reading in files...",Sys.time()))
human_21km_pairs_PCs<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/human_21_pairs_PCs.tsv")

maize_wHumanPCcentroids<-human_21km_pairs_PCs %>% 
  select(starts_with("PC_"),sample_id) %>%
  group_by(sample_id) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)))

#non-transformed genotypes
#maize_genotypes<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv")

#maize_genotypes<-select(maize_genotypes, ID, contains("SEEDGWAS"))

#to run a general linear regression model, then I need:
#genotype as a column
#each human PC as a column
#only running it on a given SNP locus

#transpose the genotypes so the SNPs are columns and individuals are rows and conserve names
#print(paste("Transposing genotypes...", Sys.time()))
#t_maize_genotypes<-t(select(maize_genotypes, contains("GWAS")))  %>% as_tibble(rownames = NA) %>% rownames_to_column()
#colnames(t_maize_genotypes)<-c("sample_id",maize_genotypes$ID)

#transformed genotypes (because they're already transposed)
t_maize_genotypes<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.genotypes.tsv")
t_maize_genotypes$sample_id<-read_lines("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.genotypes.sampleIDs.txt")

#change sample ID to match with PC sample IDs and then join with human centroid PCs
print(paste("Joining transposed genotypes to PCs...", Sys.time()))

#for non-transformed genotypes
#t_maize_genotypes<-mutate(t_maize_genotypes, sample_id = str_split(sample_id, ":", simplify = T)[,1]) %>% 
#  left_join(., maize_wHumanPCcentroids, by="sample_id")

#for transformed genotypes
t_maize_genotypes<-mutate(t_maize_genotypes, sample_id = str_split(sample_id, "\\.", simplify = T)[,1]) %>% 
  left_join(., maize_wHumanPCcentroids, by="sample_id")

#loop through each SNP and regress the human PCs on it and record the resulting r squareds
print(paste("Starting LM loop", Sys.time()))

#lm_results<-tibble(snp_ID = NA, r.squared = NA, adj.r.squared = NA)
# This for loop statement for non-transfromed genotypes: for(snp in 2:(ncol(t_maize_genotypes)-10)){
#This for loop and if statement for transformed genotypes
#for(snp in 1:ncol(t_maize_genotypes)){
#	if(str_detect(string = colnames(t_maize_genotypes)[snp], pattern = "S[0-9]")){ 
#		placeholder.summary<-summary(lm(data = t_maize_genotypes, 
#                                  pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10))
#  	lm_results<-add_row(lm_results, 
#                      snp_ID = colnames(t_maize_genotypes)[snp],
#                      r.squared = placeholder.summary$r.squared,
#                      adj.r.squared = placeholder.summary$adj.r.squared)
#  	if(snp %% 1000 == 0){print(paste("On SNP",snp,Sys.time()))}
#	}
#}

lm_results<-tibble(snp_ID = NA, r.squared = NA, adj.r.squared = NA, model = NA)
for(snp in 1:ncol(t_maize_genotypes)){
  #so long as the column is a SNP column
  if(str_detect(string = colnames(t_maize_genotypes)[snp], pattern = "S[0-9]")){
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9 + PC_10))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_10"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8 + PC_9))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_9"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7 + PC_8))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_8"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6 + PC_7))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_7"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5 + PC_6))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_6"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3 + PC_4 + PC_5))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_5"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3 + PC_4))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_4"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2 + PC_3))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_3"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1 + PC_2))
    lm_results<-add_row(lm_results, 
                        snp_ID = colnames(t_maize_genotypes)[snp],
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model = c("PC1_2"))
    
    placeholder.summary<-summary(lm(data = t_maize_genotypes, 
                                    pull(t_maize_genotypes[,snp]) ~ PC_1))
    
    lm_results<-add_row(lm_results, 
                          snp_ID = colnames(t_maize_genotypes)[snp],
                          r.squared = placeholder.summary$r.squared,
                          adj.r.squared = placeholder.summary$adj.r.squared,
                          model = c("PC1"))
    if(snp %% 1000 == 0){print(paste("On SNP",snp,Sys.time()))}
  }
}
lm_results<-na.omit(lm_results)

print(paste("Finished LM loop and now writing...",Sys.time()))

#ggplot(lm_results, aes(x=r.squared))+geom_histogram()
#ggplot(lm_results, aes(x=adj.r.squared))+geom_histogram()

#non-transformed_genotypes
#write_tsv(lm_results, "/group/jrigrp11/snodgras_maizePopulations/21km_maize/lm_results.tsv")

#transformed genotypes
write_tsv(lm_results, "/group/jrigrp11/snodgras_maizePopulations/21km_maize/cumulative_transformed_lm_results.tsv")
