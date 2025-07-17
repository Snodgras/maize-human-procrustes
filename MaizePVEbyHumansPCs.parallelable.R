#!/usr/bin/env Rscript

library(tidyverse)

args<-commandArgs(T)

#ARG 1 = transformed or not

#ARG 2 = genotype file path

#ARG 3 = FORMULA?

#ARG 4 = output file name/path

#### READ IN CENTROID FILE ####
print(paste("Reading in files...",Sys.time()))
human_21km_pairs_PCs<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/human_21_pairs_PCs.tsv")

maize_wHumanPCcentroids<-human_21km_pairs_PCs %>% 
  select(starts_with("PC_"),Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize) %>%
  group_by(Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)))

colnames(maize_wHumanPCcentroids)[1]<-"sample_id"

#### READ IN GENOTYPE FILES AND PREP FOR LM ####

if(args[1] == "raw"){
  print("Reading in raw genotypes")
  maize_genotypes<-read_tsv(args[2])
  maize_genotypes<-select(maize_genotypes, ID, contains("SEEDGWAS"))
  #transpose the genotypes so the SNPs are columns and individuals are rows and conserve names
  print(paste("Transposing genotypes...", Sys.time()))
  t_maize_genotypes<-t(select(maize_genotypes, contains("GWAS")))  %>% as_tibble(rownames = NA) %>% rownames_to_column()
  colnames(t_maize_genotypes)<-c("sample_id",maize_genotypes$ID)
  print(paste("Joining transposed genotypes to PCs...", Sys.time()))
  t_maize_genotypes<-mutate(t_maize_genotypes, sample_id = str_split(sample_id, "\\.", simplify = T)[,1]) %>% left_join(., maize_wHumanPCcentroids, by="sample_id")
}else{
  if(args[1] == "transformed"){
    print("Reading in transformed genotypes")
    #transformed genotypes (because they're already transposed)
    print("WARNING: sample IDs are hard coded in, so make sure that's OK")
    t_maize_genotypes<-read_tsv(args[2])
    t_maize_genotypes$sample_id<-read_lines("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.genotypes.sampleIDs.txt")
    print(paste("Joining transposed genotypes to PCs...", Sys.time()))
    t_maize_genotypes<-mutate(t_maize_genotypes, sample_id = str_split(sample_id, "\\.", simplify = T)[,1]) %>% left_join(., maize_wHumanPCcentroids, by="sample_id")
}else{
  print("ERROR: Don't know what type of genotypes these are! Check arguement 1")
}
  }

#### LM LOOP ####
print(paste("Starting LM loop", Sys.time()))

lm_model<-paste0("pull(t_maize_genotypes[,snp]) ~ ", args[3])

print(paste0("Using this formula: ", lm_model))

lm_results<-tibble(snp_ID = NA, r.squared = NA, adj.r.squared = NA, model_formula =NA)

for(snp in 1:ncol(t_maize_genotypes)){
	if(str_detect(string = colnames(t_maize_genotypes)[snp], pattern = "S[0-9]")){ 
		placeholder.summary<-summary(lm(data = t_maize_genotypes, lm_model))
  	lm_results<-add_row(lm_results, 
                      snp_ID = colnames(t_maize_genotypes)[snp],
                      r.squared = placeholder.summary$r.squared,
                      adj.r.squared = placeholder.summary$adj.r.squared,
                      model_formula = args[3])
  	if(snp %% 1000 == 0){print(paste("On SNP",snp,Sys.time()))}
	}
}
lm_results<-na.omit(lm_results)

print(paste("Finished LM loop and now writing...",Sys.time()))

#### WRITE OUTPUT ####
write_tsv(lm_results, args[4])


