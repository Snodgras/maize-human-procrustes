#!/usr/bin/env Rscript

library(tidyverse)

args<-commandArgs(T)

#ARG 1 = transformed or not

#ARG 2 = genotype file path

#ARG 3 = FORMULA?

#ARG 4 = output file name/path

#### READ IN CENTROID FILE ####
print(paste("Reading in files...",Sys.time()))
human_21km_pairs_PCs<-read_tsv("/group/jrigrp11/snodgras/snodgras_maizePopulations/21km_maize/human_21_pairs_PCs.tsv")

#create centriods of human PCs for maize individuals
maize_wHumanPCcentroids<-human_21km_pairs_PCs %>% 
  select(starts_with("PC_"),Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize, 
         locations_latitude.maize,
         locations_longitude.maize, 
         locations_elevation.maize) %>%
  group_by(Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize) %>% 
  summarize(across(contains(c("PC","locations")), ~ mean(.x,na.rm = T)))

colnames(maize_wHumanPCcentroids)[1]<-"sample_id"

#add in admixture estimates
global_maize_admix<-read_tsv("/group/jrigrp11/snodgras/snodgras_maizePopulations/admix_removal/admixProp.GBSsamples.txt", col_names = c("ID","Admix"))
global_maize_admix<-mutate(global_maize_admix, sample_id = str_split(ID, "\\.",simplify = T)[,1])

maize_wHumanPCcentroids<-inner_join(maize_wHumanPCcentroids, select(global_maize_admix, -ID), by="sample_id")

#gather the Mexican maize PCs that will be used as response variables
maize_meta<-read_csv("/group/jrigrp11/snodgras/snodgras_maizePopulations/clean_seedspassport.csv")
Mexmaize_pca<-read_tsv("/group/jrigrp11/snodgras/snodgras_maizePopulations/admix_removal/Mex_subset/Mex_non_transformed_LDprune_MAF0.01.pca")
Mexmaize<-mutate(Mexmaize_pca, sample_id = str_split(ind_id, ":",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

maize_wHumanPCcentroids<-inner_join(maize_wHumanPCcentroids, select(Mexmaize,c(paste0("PC",1:20), "sample_id")), by = "sample_id")

# Run the linear model using each maize PC as response variable instead of genotype
#human pcs: PC_number
#maize pcs: PCnumber

lm_model<-paste0("pull(maize_wHumanPCcentroids[,pc]) ~ ", args[3])

print(paste0("Using this formula: ", lm_model))

#scratch test
##this works
lm(data = maize_wHumanPCcentroids, PC1 ~ locations_longitude.maize+locations_latitude.maize+locations_elevation.maize+Admix)
## this doesn't work
lm(data = maize_wHumanPCcentroids, PC1 ~ locations_longitude.maize+locations_latitude.maize+locations_elevation.maize+Admix+PC_1)

#have to change the variable names so it doesn't confuse them
colnames(maize_wHumanPCcentroids)<-c("sample_id", paste0("PC_",1:10,".human"), "latitude.maize","longitude.maize","elevation.maize","Admix",paste0("PC",1:20,".maize"))

lm(data = maize_wHumanPCcentroids, PC1.maize ~ longitude.maize+latitude.maize+elevation.maize+Admix+PC_1.human)
lm_trial<-lm(data = maize_wHumanPCcentroids, 
   PC1.maize ~ longitude.maize+latitude.maize+elevation.maize+Admix+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human+PC_7.human+PC_8.human+PC_9.human+PC_10.human)

summary(lm_trial)

cor(maize_wHumanPCcentroids[,2:ncol(maize_wHumanPCcentroids)]) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  pivot_longer(cols = paste0("PC",1:20,".maize"), names_to = "maize_PC", values_to = "R") %>%
  select(rowname, maize_PC, R) %>%
  filter(!rowname %in% maize_PC) %>%
  mutate(maize_PC = factor(maize_PC, levels = c(paste0("PC",1:20,".maize"))),
         rowname = factor(rowname, levels = c("latitude.maize","longitude.maize","elevation.maize","Admix",paste0("PC_",1:10,".human"))))%>%
  ggplot(aes(x=rowname, y=maize_PC))+
  geom_tile(aes(fill = R^2))+
  geom_text(aes(label = signif(R^2, 2)))+
  colorspace::scale_fill_continuous_diverging()

lm_results<-tibble(maizePC = NA, r.squared = NA, adj.r.squared = NA, model_formula =NA)

for(n in 1:20){
  placeholder.summary<-summary(
    lm(data = maize_wHumanPCcentroids, 
       pull(maize_wHumanPCcentroids[,paste0("PC",n,".maize")]) ~ longitude.maize+latitude.maize+elevation.maize+Admix+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human+PC_7.human+PC_8.human+PC_9.human+PC_10.human)
  )
  
  lm_results<-add_row(lm_results, 
                      maizePC = n,
                      r.squared = placeholder.summary$r.squared,
                      adj.r.squared = placeholder.summary$adj.r.squared,
                      model_formula = "LLEA_10")
}

lm_results<-na.omit(lm_results)

#But what if we want to know which parameters are driving each maize model?
for(i in c(" ",
           "+PC_1.human",
           "+PC_1.human+PC_2.human",
           "+PC_1.human+PC_2.human+PC_3.human",
           "+PC_1.human+PC_2.human+PC_3.human+PC_4.human",
           "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human",
           "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human",
           "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human+PC_7.human",
           "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human+PC_7.human+PC_8.human",
           "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human+PC_7.human+PC_8.human+PC_9.human")){ #for each human PC
  for(n in 1:20){ #for each maize PC
    lm_model<-paste0("pull(maize_wHumanPCcentroids[,'PC",n,".maize']) ~ longitude.maize+latitude.maize+elevation.maize+Admix",i)
    placeholder.summary<-summary(lm(data = maize_wHumanPCcentroids, lm_model))
    if(i == " "){model_formula_text = "LLEA"}
    if(i == "+PC_1.human"){model_formula_text = "LLEA_1"}
    if(i == "+PC_1.human+PC_2.human"){model_formula_text = "LLEA_2"}
    if(i == "+PC_1.human+PC_2.human+PC_3.human"){model_formula_text = "LLEA_3"}
    if(i == "+PC_1.human+PC_2.human+PC_3.human+PC_4.human"){model_formula_text = "LLEA_4"}
    if(i == "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human"){model_formula_text = "LLEA_5"}
    if(i == "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human"){model_formula_text = "LLEA_6"}
    if(i == "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human+PC_7.human"){model_formula_text = "LLEA_7"}
    if(i == "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human+PC_7.human+PC_8.human"){model_formula_text = "LLEA_8"}
    if(i == "+PC_1.human+PC_2.human+PC_3.human+PC_4.human+PC_5.human+PC_6.human+PC_7.human+PC_8.human+PC_9.human"){model_formula_text = "LLEA_9"}
    lm_results<-add_row(lm_results, 
                        maizePC = n,
                        r.squared = placeholder.summary$r.squared,
                        adj.r.squared = placeholder.summary$adj.r.squared,
                        model_formula = model_formula_text)
    }
}

lm_results$model_formula<-factor(lm_results$model_formula,levels = c("LLEA",paste0("LLEA_",1:10)))
lm_results$maizePC<-factor(lm_results$maizePC,levels = c(1:20))

ggplot(lm_results, aes(x=model_formula, y=adj.r.squared))+
  geom_line(aes(group = maizePC))+
  geom_point(aes(color = maizePC))+
  geom_text(aes(label = maizePC, color = maizePC), nudge_y = 0.02)+
  scale_color_viridis_d(direction = -1, option = "H")+
  theme_bw()

ggplot(lm_results, aes(x=model_formula, y=maizePC))+
  geom_tile(aes(fill = adj.r.squared))+
  geom_label(aes(label = signif(adj.r.squared, 2)))+
  theme_bw()+
  scale_fill_viridis_c()+
  xlab("Model Formula")+
  ggtitle("Mex. Maize PCs as Response Variable")
ggsave("/group/jrigrp11/snodgras/snodgras_maizePopulations/Plots/2026-06-08-MexMaizePCsResponseVar-heatmap.pdf",
       device = "pdf", dpi = 300, width = 2400, height = 2400, units = "px")

filter(lm_results, maizePC %in% c(12, 8, 9, 6, 2)) %>% #key PCs to look at
  ggplot(aes(x=model_formula, y=adj.r.squared))+
  geom_line(aes(group = maizePC, color = maizePC))+
  geom_point(aes(color = maizePC))+
  #geom_text(aes(label = maizePC, color = maizePC))+
  scale_color_viridis_d(direction = -1, option = "H")+
  theme_bw()

# quantitative cutoff for looking at changes? 
filter(lm_results, model_formula %in% c("LLEA", "LLEA_10")) %>%
  pivot_wider(id_cols = maizePC, names_from = model_formula, values_from = adj.r.squared) %>%
  mutate(adj.r.squared.change = LLEA_10 - LLEA) %>% 
  summary()
#top quarter show a 0.18 increase
filter(lm_results, model_formula %in% c("LLEA", "LLEA_10")) %>%
  pivot_wider(id_cols = maizePC, names_from = model_formula, values_from = adj.r.squared) %>%
  mutate(adj.r.squared.change = LLEA_10 - LLEA) %>%
  filter(adj.r.squared.change > 0.18) #5, 8, 9, 12, 14, 15

filter(lm_results, maizePC %in% c(5, 8, 9, 12, 14, 15)) %>% #key PCs to look at
  ggplot(aes(x=model_formula, y=adj.r.squared))+
  geom_line(aes(group = maizePC, color = maizePC))+
  geom_point(aes(color = maizePC))+
  scale_color_viridis_d(direction = -1, option = "H")+
  theme_bw()+
  xlab("Model Formula")+ylab("Adj. R2")+
  ggtitle("Mex. Maize PCs as response variables")
ggsave("/group/jrigrp11/snodgras/snodgras_maizePopulations/Plots/2026-06-08-MexMaizePCsResponseVar-linegraph.pdf",
       device = "pdf", width = 2400, height = 1200, units = "px")

#for(snp in 1:ncol(t_maize_genotypes)){
#	if(str_detect(string = colnames(t_maize_genotypes)[snp], pattern = "S[0-9]")){ 
#		placeholder.summary<-summary(lm(data = t_maize_genotypes, lm_model))
#  	lm_results<-add_row(lm_results, 
#                      snp_ID = colnames(t_maize_genotypes)[snp],
#                      r.squared = placeholder.summary$r.squared,
#                      adj.r.squared = placeholder.summary$adj.r.squared,
#                      model_formula = args[3])
#  	if(snp %% 1000 == 0){print(paste("On SNP",snp,Sys.time()))}
#	}
#}
#lm_results<-na.omit(lm_results)

#print(paste("Finished LM loop and now writing...",Sys.time()))

#### WRITE OUTPUT ####
#write_tsv(lm_results, args[4])



