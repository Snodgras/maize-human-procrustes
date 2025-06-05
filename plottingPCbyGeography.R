#libraries
library(tidyverse)
library(vegan)

#### 0. Import data ####

#human pc data
human_pcs<-read_csv("/group/jrigrp11/snodgras_maizePopulations/human_PCs/all-indigenous.pca.csv")

#human location data
human_meta<-read_csv("/group/jrigrp11/snodgras_maizePopulations/human_PCs/indigenous-american-genomes-popinfo.csv")

#joined dataframe
human<-inner_join(human_pcs, human_meta, by = c("sample_name" = "sample_id" ))

#### 1. Identify outliers and remove ####

### HUMAN DATA ###

human.outlier.stats<-human %>% group_by(Region) %>% summarize(mean.PC1 = mean(PC_1, na.rm = T),
                                                              mean.PC2 = mean(PC_2, na.rm = T),
                                                              sd.PC1 = sd(PC_1, na.rm = T),
                                                              sd.PC2 = sd(PC_2, na.rm = T))
#filter out any samples where they're beyond 3 sd from mean for their regional group
outlier<-inner_join(human, human.outlier.stats, by = "Region") %>% 
  mutate(PC1_lower.cutoff = mean.PC1 - 3*sd.PC1, 
         PC1_upper.cutoff = mean.PC1 + 3*sd.PC1,
         PC2_lower.cutoff = mean.PC2 - 3*sd.PC2, 
         PC2_upper.cutoff = mean.PC2 + 3*sd.PC2) %>%
  filter(PC_1 < PC1_lower.cutoff | PC_1 > PC1_upper.cutoff | PC_2 < PC2_lower.cutoff | PC_2 > PC2_upper.cutoff)

outlier %>% ggplot(aes(x=longitude, y=latitude, color = Region))+geom_point()
outlier %>% ggplot(aes(x=PC_1, y=PC_2, color = Region))+geom_point()

human_noOutliers<-filter(human, !sample_name %in% outlier$sample_name)
# remove recently admixed human samples
#these sample IDs come from Santiago 9/4/24
human_noOutliers<-filter(human_noOutliers, !sample_name %in% c("GA006377","GA006434","GA004777","GA006473",
                                                               "GA004778","GA006357","GA004764","GA004803","GA006340","GA004774",
                                                               'GA006557',"GA004768","GA006380","GA006447","GA004798","GA006553"))



#plot location
ggplot(human_noOutliers, aes(x=longitude,y=latitude))+
  geom_point()+
  theme_bw()+
  xlab("Longitude")+ylab("Latitude")+
  ggtitle("Indigenous human location")

### MAIZE DATA ###
#PCA from non-transformed data
maize_pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/non_transformed_LDprune_MAF0.01.pca")
maize_eigenval<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/non_transformed_LDprune_MAF.01.eigens")

#pve<- data.frame(PC = 1:length(maize_eigenval), pve = (maize_eigenval/sum(maize_eigenval))*100)

#load in the maize location data
maize_meta<-read_csv("/group/jrigrp11/snodgras_maizePopulations/clean_seedspassport.csv")
ggplot(maize_meta, aes(x=locations_longitude,y=locations_latitude))+
  geom_point()+
  theme_bw()+
  xlab("Longitude")+ylab("Latitude")+
  ggtitle("Maize Location")

colnames(maize_pca)[1]<-"ind"

#join meta and pcs
maize<-mutate(maize_pca, sample_id = str_split(ind_id, ":",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

#maize.outlier.stats<-maize %>% group_by(countries_country_name) %>% summarize(mean.PC1 = mean(PC1, na.rm = T),
#                                                                              mean.PC2 = mean(PC2, na.rm = T),
#                                                                              sd.PC1 = sd(PC1, na.rm = T),
#                                                                              sd.PC2 = sd(PC2, na.rm = T))
#filter out any samples where they're beyond 3 sd from mean for their regional group
#maize.outlier<-inner_join(maize, maize.outlier.stats, by = "countries_country_name") %>% 
#  mutate(PC1_lower.cutoff = mean.PC1 - 3*sd.PC1, 
#         PC1_upper.cutoff = mean.PC1 + 3*sd.PC1,
#         PC2_lower.cutoff = mean.PC2 - 3*sd.PC2, 
#         PC2_upper.cutoff = mean.PC2 + 3*sd.PC2) %>%
#  filter(PC1 < PC1_lower.cutoff | PC1 > PC1_upper.cutoff | PC2 < PC2_lower.cutoff | PC2 > PC2_upper.cutoff)

#maize.outlier %>% ggplot(aes(x=locations_longitude, y=locations_latitude, color = countries_country_name))+geom_point()
#maize.outlier %>% ggplot(aes(x=PC1, y=PC2, color = countries_country_name))+geom_point()

#maize_noOutliers<-filter(maize, !ind %in% maize.outlier$ind)


### maize Mexican only ###

Mexmaize_pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/Mex_subset/Mex_non_transformed_LDprune_MAF0.01.pca")
Mexmaize_eigenval<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/Mex_subset/Mex_non_transformed_LDprune_MAF.01.eigens") 

#Mexpve<- data.frame(PC = 1:length(Mexmaize_eigenval), pve = (Mexmaize_eigenval/sum(Mexmaize_eigenval))*100)

#join meta and pcs
Mexmaize<-mutate(Mexmaize_pca, sample_id = str_split(ind_id, ":",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

#### Pair humans with maize PC centroids based on 20km buffer ####
ind_human_maize_20km_buffer<-read_csv("/group/jrigrp11/snodgras_maizePopulations/Individual-human-20km-maize-samples-joined.csv")
colnames(ind_human_maize_20km_buffer)<-c("sample_id.human","Ethnicity.human","Country.human","Region.human","latitude.human","longitude.human","Indigenous_Proportion.human","id.maize","general_identifier.maize",
                                         "name.maize","bank_number.maize","taxonomy_id.maize","collnumb.maize","colldate.maize","location_id.maize","created_on.maize","locations_id.maize","locations_region.maize",
                                         "locations_site_name.maize","locations_elevation.maize","locations_latitude.maize","locations_longitude.maize","countries_id.maize","countries_country_code2.maize",
                                         "countries_country_code3.maize","countries_country_name.maize","taxonomies_id.maize","taxonomies_genus.maize","taxonomies_species.maize","taxonomies_crop_name.maize",
                                         "taxonomies_ploidy.maize","Sample_ID_of_DNA_from_composite_samples.maize","Sample_ID_of_DNA_from_most_recent_CML_regenerations.maize",
                                         "Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize","Tester_GID.maize","Tester_pedigree.maize","Testcross_GID.maize","PrimaryRace.maize","PrimaryPurity.maize",
                                         "SecondaryRace.maize","Pedigree.maize","GrainType1.maize","GrainType2.maize","GrainType3.maize","GrainColor1.maize","GrainColor2.maize","GrainColor3.maize","PopulationType.maize")
ind_human_maize_20km_buffer %>% group_by(sample_id) %>% count() %>% 
  ggplot(aes(x=n))+
  geom_histogram(binwidth = 1)

ggplot(ind_human_maize_20km_buffer)+
  geom_point(aes(x=locations_longitude.maize, y=locations_latitude.maize), shape = 3)+
  geom_point(aes(x=longitude.human,y=latitude.human, color = Region.human), alpha = 0.6)+
  theme_bw()+xlab("Longitude")+ylab("Latitude")+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c( "Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_maizeAndHumans_geography.png",device = "png",dpi=300, width = 6, height = 4.5)

#write out the ids so that pca can be rerun on just those samples
selected_passport_ids<-read_csv("/group/jrigrp11/snodgras_maizePopulations/selected_genotypeIDs.csv", col_names = c("vcf_sample_id","gwas_sample_id"))
filter(selected_passport_ids, gwas_sample_id %in% ind_human_maize_20km_buffer$Sample_ID_of_DNA_from_single_plants_used_in_GWAS) %>% 
  write_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_maize_samples_ids.tsv")

#after running pca on maize samples

ind_human_maize_20km_buffer.pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_non_transformed_LDprune_MAF0.01.pca")
ind_human_maize_20km_buffer.eigenval<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_non_transformed_LDprune_MAF.01.eigens")
colnames(ind_human_maize_20km_buffer.eigenval)<-c("PC","PVE")

ggplot(ind_human_maize_20km_buffer.pca, aes(x=PC1, y=PC2))+
  geom_point()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(ind_human_maize_20km_buffer.eigenval$PVE[1]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(ind_human_maize_20km_buffer.eigenval$PVE[2]*100,3),"%)"))+
  ggtitle("21km radius maize PCA")

ind_human_maize_20km_buffer<-left_join(x=selected_passport_ids, y=ind_human_maize_20km_buffer.pca, by=c("gwas_sample_id"="ind_id")) %>% 
  right_join(x=., y=ind_human_maize_20km_buffer, by=c("gwas_sample_id"="Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize"))

ggplot(ind_human_maize_20km_buffer, aes(x=PC1, y=PC2))+
  geom_point(aes(color=locations_elevation.maize))+theme_bw()+
  scale_color_viridis_c(name="Elevation")+
  xlab(paste0("PC1 (", signif(ind_human_maize_20km_buffer.eigenval$PVE[1]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(ind_human_maize_20km_buffer.eigenval$PVE[2]*100,3),"%)"))+
  ggtitle("21km radius maize PCA")

#Goal is to find the centroid of PCs for each set of maize samples for a given human sample
centroidPCs_ind_human<-ind_human_maize_20km_buffer %>% 
  group_by(sample_id.human,Region.human,Ethnicity.human,latitude.human,longitude.human) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)),
            mean_elevation.maize=mean(locations_elevation.maize))
colnames(centroidPCs_ind_human)<-c("sample_name","Region","Ethnicity","latitude","longitude",
                                   paste0("PC",1:422,".maize"), "mean_elevation.maize")
#add in the human PCs, specifically the Mexican only PCs
centroidPCs_ind_human<-inner_join(centroidPCs_ind_human, Mexhuman, by=("sample_name"))

#### 2. PCA calculations and figure quality plots ####

### Human data full ### 

ggplot(human_noOutliers, aes(x=longitude, y=latitude))+
  geom_point(aes(color = Region), shape = 3)+
  theme_bw()+theme(text = element_text(size=20))+
  xlab("Longitude")+ylab("Latitude")+
  ggtitle("Human Samples Location")+
  scale_color_manual(
    values = c("#538fff","#ff5b58","#ecc500","#28d2ab","#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957","#ec00c5"),
    name = "Region", labels = c("Amazonia", "Andean Highland","Central South America",
                                "Chaco Amerindian", "Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West","Patagonia"))+
  guides(color = "none")
  #guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/HumanSamples_byRegion_LatLong.png",
       device = "png", height = 5.1, width = 5.4, dpi = 300)

ggplot(human_noOutliers, aes(x=PC_1, y=PC_2))+
  geom_point(aes(color = Region))+
  theme_bw()+theme(text = element_text(size=20))+
  xlab("PC 1")+ylab("PC 2")+
  ggtitle("Human PCA")+
  scale_color_manual(
    values = c("#538fff","#ff5b58","#ecc500","#28d2ab","#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957","#ec00c5"),
    name = "Region", labels = c("Amazonia", "Andean Highland","Central South America",
                                "Chaco Amerindian", "Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West","Patagonia"))+
  guides(color = "none")
#guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/HumanSamples_byRegion_PCA.png",
       device = "png", height = 5.1, width = 5.4, dpi = 300)


### Human Mexican only ###
#Mexican subset of humans
Mexhuman_pca<-read_csv("/group/jrigrp11/snodgras_maizePopulations/human_PCs/all-mexican.pca.csv")
Mexhuman<-inner_join(Mexhuman_pca, human_meta, by = c("sample_name" = "sample_id" ))

ggplot(Mexhuman, aes(x=PC_1, y=PC_2))+
  geom_point(aes(color = Region), alpha = 0.6)+
  theme_bw()+theme(text = element_text(size=20))+
  xlab("PC 1")+ylab("PC 2")+
  ggtitle("Human PCA")+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))#guides(color = "none")

### Maize full data ###
#plot pca
ggplot(maize, aes(x=PC1, y=PC2))+
  geom_point(aes(color = countries_country_name))+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(maize_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(maize_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Maize PCA")
#test function
plot_maize_pca(maize, "Maize PCA", maize_eigenval[1,2], maize_eigenval[2,2])

#plot pca PC2 vs. PC3
ggplot(maize, aes(x=PC2, y=PC3))+
  geom_point(aes(color = countries_country_name))+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC2 (", signif(maize_eigenval[2,2],3),"%)"))+
  ylab(paste0("PC3 (", signif(maize_eigenval[3,2],3),"%)"))+
  ggtitle("Maize PCA, PC2 vs PC3")

### Maize Mexican only ###
plot_maize_pca(Mexmaize, "Mexican Maize PCA", Mexmaize_eigenval[1,2],Mexmaize_eigenval[2,2])
ggplot(Mexmaize, aes(x=PC1, y=PC2))+
  geom_point(aes(color = locations_elevation))+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(Mexmaize_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(Mexmaize_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Mexican Maize PCA")+ scale_color_viridis_c(name = "Elevation (m)")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/Mexmaize_PCA_PC1_PC2_byElevation.png",
       device = "png", dpi= 300, width = 6, height = 6)

ggplot(Mexmaize, aes(x=PC2, y=PC3))+
  geom_point(aes(color = locations_elevation))+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC2 (", signif(Mexmaize_eigenval[2,2],3),"%)"))+
  ylab(paste0("PC3 (", signif(Mexmaize_eigenval[3,2],3),"%)"))+
  ggtitle("Mexican Maize PCA, PC2 vs PC3")+ scale_color_viridis_c()

#yucatan rough polygon: North = 22N lat, South =18N lat, East = 88W long, West = 90 long
mutate(Mexmaize, Yucatan = case_when(locations_latitude <= 22 & locations_latitude >= 18 & 
                                       locations_longitude >= -93 ~ TRUE,
                                     .default = FALSE)) %>%
  ggplot(aes(x=locations_longitude, y=locations_latitude, color = Yucatan))+geom_point()+theme_bw()
mutate(Mexmaize, Yucatan = case_when(locations_latitude <= 22 & locations_latitude >= 18 & 
                                       locations_longitude >= -93 ~ TRUE,
                                     .default = FALSE)) %>%
  ggplot(aes(x=PC1, y=PC2,color = Yucatan))+
  geom_point(aes(), alpha =0.6)+
  theme_bw()+
  xlab(paste0("PC1 (", signif(Mexmaize_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(Mexmaize_eigenval[2,2]*100,3),"%)"))

mutate(maize, Yucatan = case_when(locations_latitude <= 22 & locations_latitude >= 18 & 
                                   locations_longitude >= -93 & locations_longitude <= -86 ~ TRUE,
                                 .default = FALSE)) %>%
  ggplot(aes(x=locations_longitude, y=locations_latitude, color = Yucatan))+geom_point()+theme_bw()
mutate(maize, Yucatan = case_when(locations_latitude <= 22 & locations_latitude >= 18 & 
                                    locations_longitude >= -93 & locations_longitude <= -86 ~ TRUE,
                                  .default = FALSE)) %>%
  ggplot(aes(x=PC1, y=PC2,color = Yucatan))+
  geom_point(aes(), alpha =0.6)+
  theme_bw()+
  xlab(paste0("PC1 (", signif(maize_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(maize_eigenval[2,2]*100,3),"%)"))+ggtitle("Non-transformed, all samples")

mutate(maize_admixRemoved, Yucatan = case_when(locations_latitude <= 22 & locations_latitude >= 18 & 
                                    locations_longitude >= -93 & locations_longitude <= -86 ~ TRUE,
                                  .default = FALSE)) %>%
  ggplot(aes(x=PC1, y=PC2,color = Yucatan))+
  geom_point(aes(), alpha =0.6)+
  theme_bw()+
  xlab(paste0("PC1 (", signif(maize_admixRemoved_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(maize_admixRemoved_eigenval[2,2]*100,3),"%)"))+ggtitle("Mexicana admixture removed, all samples")

#maybe try plotting with respect to grain type?
ggplot(Mexmaize, aes(x=PC1, y=PC2))+
  geom_point(aes(color = GrainType1), alpha = 0.75)+
  #coord_equal()+
  theme_bw()+
  #xlab(paste0("PC1 (", signif(maize_eigenval[1,2]*100,3),"%)"))+
  #ylab(paste0("PC2 (", signif(maize_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Maize by Grain Type 1 PCA")+ scale_color_manual(values = c("Dent"="#3C91E6", "Flint"="#FA824C",
                                                                      "Floury"="#646F4B","Popcorn"="#CB429F",
                                                                      "Semi-Dent"="#10477E","Semi-flint"="#B33C05",
                                                                      "Sweet"="#C8AB2D","NA"="#342E37"))
#GrainType3 is all NAs
maize %>% group_by(GrainType1, GrainType2) %>% count() %>% View()
#the majority are dents with some flints and a handful of other types

#### Centroid 21 km pairings ####
#plot the maize (centroid) PCs
ggplot(centroidPCs_ind_human, aes(x=PC1.maize,y=PC2.maize))+
  geom_point(aes(color=mean_elevation.maize))+
  theme_bw()+scale_color_viridis_c(name="Mean Elevation")+
  ggtitle("21km radius maize centroids")

#plot the maize (centroid) PCs colored by human regions
ggplot(centroidPCs_ind_human, aes(x=PC1.maize,y=PC2.maize))+
  geom_point(aes(color=Region))+
  theme_bw()+ theme(legend.text=element_text(size=12))+
  xlab("maize PC 1")+ylab("maize PC 2")+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1aff",
               "#b100ea","#386955","#a7c957"),
    name = "Region", labels = c(
      "Mexico, Center","Mexico, Gulf", 
      "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
      "Mexico, Oaxaca","Mexico, West"))

#### 3. Mexicana admixture removal from maize ####

#because some maize samples don't have metadata (location info)
#write the sample IDs that do have metadata for filtering ahead of PCA
select(maize_meta, Sample_ID_of_DNA_from_single_plants_used_in_GWAS) %>% pull() %>% write_lines("/group/jrigrp11/snodgras_maizePopulations/clean_seedspassport.sampleIDsOnly.txt")

maize_admixRemoved_pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/transformed_LDprune_MAF0.01.pca")
maize_admixRemoved_eigenval<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/transformed_LDprune_MAF0.01.eigens") 

#join meta and pcs
maize_admixRemoved<-mutate(maize_admixRemoved_pca, sample_id = str_split(ind_id, "\\.",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

#plot the PCA
#plot_maize_pca(maize_admixRemoved, "Maize with mexicana anc. removed PCA", maize_admixRemoved_pve)
ggplot(maize_admixRemoved, aes(x=PC1, y=PC2))+
  geom_point(aes(color = countries_country_name), alpha=0.5)+
  theme_bw()+#scale_color_viridis_c(name = "Elevation (m)")+
  xlab(paste0("PC1 (", signif(maize_admixRemoved_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(maize_admixRemoved_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Maize with mexicana anc. removed PCA")+
  theme(text = element_text(size = 20))

ggplot(maize_admixRemoved, aes(x=PC1, y=PC2))+
  geom_point(aes(color = locations_elevation), alpha=0.5)+
  theme_bw()+scale_color_viridis_c(name = "Elevation (m)")+
  xlab(paste0("PC1 (", signif(maize_admixRemoved_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(maize_admixRemoved_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Maize with mexicana anc. removed PCA")+
  theme(text = element_text(size = 20))
#ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/maize_admixtureRemoved.PC1PC2.byElevation.png",
#       height = 7, width = 8, dpi=300, device = "png")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/maize_admixtureRemoved.PC1PC2.byElevation.allSNPs.png",
       height = 7, width = 8, dpi=300, device = "png")

ggplot(maize_admixRemoved, aes(x=PC2, y=PC3))+
  geom_point(aes(color = countries_country_name), alpha=0.5)+
  theme_bw()+#scale_color_viridis_c()+
  xlab(paste0("PC2 (", signif(maize_admixRemoved_eigenval[2,2]*100,3),"%)"))+
  ylab(paste0("PC3 (", signif(maize_admixRemoved_eigenval[3,2]*100,3),"%)"))+
  ggtitle("Maize with mexicana anc. removed PCA: PC2 vs PC3")

#double checking admixture proportions by geography
global_maize_admix<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/admixProp.GBSsamples.txt", col_names = c("ID","Admix"))
global_maize_admix<-mutate(global_maize_admix, sample_id = str_split(ID, "\\.",simplify = T)[,1])

inner_join(maize_admixRemoved, global_maize_admix, by = "sample_id") %>% 
  ggplot(aes(x=locations_longitude, y=locations_latitude, color = Admix))+
  geom_point(alpha = 0.5)+
  theme_bw()+scale_color_viridis_c(name = "Mexicana \nAdmixture \nProportion")+
  xlab("Longitude")+ylab("Latitude")+
  theme(text=element_text(size=20))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/mexicanaAdmixtureProportion_byGeography.png",
       device="png",width=7,height = 5.5, dpi = 300)

#plot PC1 by geography
ggplot(maize_admixRemoved, aes(x=locations_longitude, y=locations_latitude, color = PC2))+
  geom_point(alpha=0.5)+
  theme_bw()+scale_color_viridis_c()

#correlation between PC1 and PC2
cor.test(maize_admixRemoved$PC1, maize_admixRemoved$PC2)
#correlation between PC2 and PC3
cor.test(maize_admixRemoved$PC3, maize_admixRemoved$PC2)

#correlation between PC1 and global admixture
temp<-inner_join(maize_admixRemoved, global_maize_admix, by = "sample_id") %>% 
  select(Admix, PC1, PC2, locations_elevation) 
cor.test(temp$Admix,temp$PC1)
cor.test(temp$Admix,temp$PC2)
cor.test(temp$locations_elevation,temp$PC1)
cor.test(temp$locations_elevation,temp$PC2)

####This section would be the same as the maize object ####
#test against untransformed maize (same cutoffs)
#This is with the subset of samples with meta data for lat/long
maize_nontransformed_pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/subset_with_LatLong/nontransformed_LDprune_MAF0.01_SubsetZeaGBS.pca")
maize_nontransformed_eigenval<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/subset_with_LatLong/nontransformed_LDprune_MAF0.01_SubsetZeaGBS.eigenvec") 
#unneeded when using the importance summary writing 
#maize_nontransformed_pve<- data.frame(PC = 1:length(maize_nontransformed_eigenval), pve = (maize_nontransformed_eigenval/sum(maize_nontransformed_eigenval))*100)

#the filtering of incomplete meta data removes individuals (4845 --> 2894)
maize_nontransformed<-mutate(maize_nontransformed_pca, sample_id = str_split(rowname, "\\.",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))
maize_nontransformed<-inner_join(maize_nontransformed, global_maize_admix, by = "sample_id")
#but if we do a left join so we keep all individuals even if they don't have filtered meta data
#maize_nontransformed<-mutate(maize_nontransformed_pca, sample_id = str_split(rowname, "\\.",simplify = T)[,1]) %>%
#  left_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

ggplot(maize_nontransformed, aes(x=PC1, y=PC2))+
  geom_point(aes(color = locations_elevation), alpha=0.5)+
  theme_bw()+scale_color_viridis_c(name = "Elevation (m)")+
  xlab(paste0("PC1 (", signif(maize_nontransformed_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(maize_nontransformed_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Maize non-transformed PCA")+
  theme(text = element_text(size = 20))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/maize_nontransformed.PC1PC2.byElevation.png",
       height = 7, width = 8, dpi=300, device = "png")

ggplot(maize_nontransformed, aes(x=PC2, y=PC3))+
  geom_point(aes(color = countries_country_name), alpha=0.5)+
  theme_bw()+#scale_color_viridis_c()+
  xlab(paste0("PC2 (", signif(maize_nontransformed_eigenval[2,2]*100,3),"%)"))+
  ylab(paste0("PC3 (", signif(maize_nontransformed_eigenval[3,2]*100,3),"%)"))+
  ggtitle("Maize non-transformed PCA")

#correlation between PC1 and global admixture
temp<-inner_join(maize_nontransformed, global_maize_admix, by = "sample_id") %>% 
  select(Admix, PC1, PC2, locations_elevation) 
cor.test(temp$Admix,temp$PC1)
cor.test(temp$Admix,temp$PC2)
cor.test(temp$locations_elevation,temp$PC1)
cor.test(temp$locations_elevation,temp$PC2)

#### 5. Procrustes ####

### Human vs. geography- full set minus outliers ###

#compute procrustes using vegan 2.6-4
human_procrustes<-procrustes(Y = select(human_noOutliers, c("PC_1","PC_2")),
                             X = select(human_noOutliers, c("longitude","latitude"))
)
summary(human_procrustes)
plot(human_procrustes)
plot(human_procrustes, to.target = F, ar.col = NA)
#test the significance of two configurations
protest(X = select(human_noOutliers, c("PC_1","PC_2")),
        Y = select(human_noOutliers, c("longitude","latitude")))

#human custom plots:
#numbers specified here should be number of samples
human.translation<-matrix(rep(human_procrustes$translation,797),nrow=797,ncol=2,byrow =TRUE) 
human_noOutliers<-human_noOutliers %>% add_column(as.data.frame(human_procrustes$Yrot))
colnames(human_noOutliers)[18:19]<-c("rotated_PC1","rotated_PC2")

human_noOutliers<-add_column(human_noOutliers, 
                             transformed_PC1 = (human_procrustes$scale*t((human_procrustes$rotation) %*% t(as.matrix(human_noOutliers[,c("PC_1","PC_2")])))+human.translation)[,1],
                             transformed_PC2 = (human_procrustes$scale*t((human_procrustes$rotation) %*% t(as.matrix(human_noOutliers[,c("PC_1","PC_2")])))+human.translation)[,2])


human_noOutliers%>% 
  ggplot(aes(x=transformed_PC1, y=transformed_PC2, color = Region))+ 
  geom_point(alpha = 0.5)+
  geom_point(aes(x=longitude,y=latitude), shape = 3, color = "black")+
  theme_bw()+
  xlab("Longitude")+ylab("Latitude")+
  scale_color_manual(#values = c("#CA1551","#FB4D3D","#F38939","#EAC435",
    #           "#386641","#62747D","#1C949D","#03CEA4","#473bf0","#6A994E","#A7C957",
    #           "#7F3773"),
    values = c("#538fff","#ff5b58","#ecc500","#28d2ab","#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957","#ec00c5"),
    name = "Region", labels = c("Amazonia", "Andean Highland","Central South America",
                                "Chaco Amerindian", "Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West","Patagonia"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave('/group/jrigrp11/snodgras_maizePopulations/Plots/human_PConGeography.woAdmixedHumans.png', device = "png",dpi = 300,
       height = 6, width = 7.5)

human_noOutliers%>% 
  ggplot(aes(x=transformed_PC1, y=transformed_PC2, color = Region))+ 
  geom_point(aes(x=longitude,y=latitude), shape = 3, color = "black")+
  geom_point(alpha = 0.5)+
  theme_bw()+
  xlab("Longitude")+ylab("Latitude")+
  scale_color_manual(
    values = c("#538fff","#ff5b58","#ecc500","#28d2ab","#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957","#ec00c5"),
    name = "Region", labels = c("Amazonia", "Andean Highland","Central South America",
                                "Chaco Amerindian", "Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West","Patagonia"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  ggtitle("Human Geography-PCA Procrustes")+theme(text=element_text(size=20))

ggsave('/group/jrigrp11/snodgras_maizePopulations/Plots/HumanSamples_Procrustes_PConGeography.woAdmixedHumans.png', 
       device = "png",dpi = 300,
       height = 6, width = 8)

human_procrustes_rotation<-human_procrustes$rotation

#slope = (y2-y1)/(x2-x1)

rotated_axes_slopes<-function(rot_mat){
  #create origin axes
  x_point.1 = c(-1,0)
  x_point.2 = c(1,0)
  y_point.1 = c(0,-1)
  y_point.2 = c(0,1)
  
  #rotate points around origin
  x_prime.1 = rot_mat %*% x_point.1  
  x_prime.2 = rot_mat %*% x_point.2 
  y_prime.1 = rot_mat %*% y_point.1 
  y_prime.2 = rot_mat %*% y_point.2 
 
  #calculate slopes
  x_prime.slope = (x_prime.2[2,] - x_prime.1[2,])/(x_prime.2[1,]-x_prime.1[1,])
  y_prime.slope = (y_prime.2[2,] - y_prime.1[2,])/(y_prime.2[1,]-y_prime.1[1,])
  
  return(c(x_prime.slope, y_prime.slope))
}
rotated_axes_slopes(human_procrustes_rotation)

human_rot_PCA<-ggplot()+
  geom_abline(slope = rotated_axes_slopes(human_procrustes_rotation)[1])+
  geom_abline(slope = rotated_axes_slopes(human_procrustes_rotation)[2])+
  #geom_segment(data = as.data.frame(t(axes)*human_procrustes$scale+human_procrustes$translation) %>% mutate(axis_name = c("PC1.orig","PC2.orig")),
  #             aes(x=V1, xend=V3,y=V2,yend=V4, color = axis_name))+
  #geom_segment(data = as.data.frame(human_rot_axes*human_procrustes$scale+human_procrustes$translation) %>% mutate(axis_name = c("PC1.rot","PC2.rot")),
  #             aes(x=V1,y=V2, xend=V3, yend=V4, color = axis_name))+
  geom_point(data = human_noOutliers, aes(x=rotated_PC1, y=rotated_PC2, color = Region)) 
ggsave(plot = human_rot_PCA, "/group/jrigrp11/snodgras_maizePopulations/Plots/human_rot_PCAonGeography.png", device ="png",dpi=300, width = 6, height = 5)  

### Human Mexican only ###
Mexhuman_procrustes.proj<-procrustes(Y = select(Mexhuman, c("PC_1","PC_2")),
                                     X = select(Mexhuman, c("proj.longitude","proj.latitude"))
)
plot(Mexhuman_procrustes.proj)
protest(Y = select(Mexhuman, c("PC_1","PC_2")),
        X = select(Mexhuman, c("proj.longitude","proj.latitude")))

### Maize full data vs. geography ###
#compute procrustes using vegan 2.6-4

####Original maize procrustes vs geography####
maize_procrustes<-procrustes(Y = select(maize, c("PC1","PC2")),
                             X = select(maize, c("locations_longitude","locations_latitude"))
)
#summary(maize_procrustes)
#plot(maize_procrustes)
#plot(maize_procrustes, to.target = F, ar.col = NA)
#test the significance of two configurations
protest(Y = select(maize, c("PC1","PC2")),
        X = select(maize, c("locations_longitude","locations_latitude")))
my.translation<-matrix(rep(maize_procrustes$translation,2909),nrow=2909,ncol=2,byrow =TRUE) 

plot(maize_procrustes$scale*t((maize_procrustes$rotation) %*% t(as.matrix(maize[,c("PC1","PC2")])))+my.translation)

maize<-mutate(maize, transformed_PC1 = (maize_procrustes$scale*t((maize_procrustes$rotation) %*% t(as.matrix(maize[,c("PC1","PC2")])))+my.translation)[,1],
                  transformed_PC2 = (maize_procrustes$scale*t((maize_procrustes$rotation) %*% t(as.matrix(maize[,c("PC1","PC2")])))+my.translation)[,2])

maize<-mutate(maize, Regions = case_when(countries_country_name %in% c("BOLIVIA","CHILE","PERU") ~ "AndeanHighland",
                                                                   countries_country_name %in% c("BRAZIL") ~ "Amazonia",
                                                                   countries_country_name %in% c("COLOMBIA","ECUADOR","FRENCH GUIANA","HONDURAS","SURINAME","VENEZUELA") ~ "CentralSouthAmerica",
                                                                   countries_country_name %in% c("PARAGUAY","URUGUAY") ~ "ChacoAmerindian",
                                                                   countries_country_name %in% c("MEXICO") ~ "Mexico",
                                                                   countries_country_name %in% c("ARGENTINA") ~ "Patagonia",
                                                                   countries_country_name %in% c("COSTA RICA","EL SALVADOR","GUATEMALA",
                                                                                                 "NICARAGUA","PANAMA") ~ "Mesoamerica",
                                                                   countries_country_name %in% c("ANTIGUA AND BARBUDA","BARBADOS","CUBA",
                                                                                                 "DOMINICAN REPUBLIC","GRENADA","GUADELOUPE",
                                                                                                 "HAITI","MARTINIQUE","PUERTO RICO","SAINT VINCENT AND THE GRENADINES",
                                                                                                 "TRINIDAD AND TOBAGO","VIRGIN ISLANDS (BRITISH)","VIRGIN ISLANDS (U.S.)") ~ "Caribbean"))

ggplot(maize, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2, color = Regions), alpha = 0.5)+
  theme_bw()+scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))


protest(Y = select(maize, c("PC2","PC3")),
        X = select(maize, c("locations_longitude","locations_latitude")))
maize.2v3_procrustes<-procrustes(Y = select(maize, c("PC2","PC3")),
                             X = select(maize, c("locations_longitude","locations_latitude"))
)
my.translation.2v3<-matrix(rep(maize.2v3_procrustes$translation,2909),nrow=2909,ncol=2,byrow =TRUE) 
maize<-mutate(maize, 
                  transformed_PC2.2v3 = (maize.2v3_procrustes$scale*t((maize.2v3_procrustes$rotation) %*% t(as.matrix(maize[,c("PC2","PC3")])))+my.translation.2v3)[,1],
                  transformed_PC3.2v3 = (maize.2v3_procrustes$scale*t((maize.2v3_procrustes$rotation) %*% t(as.matrix(maize[,c("PC2","PC3")])))+my.translation.2v3)[,2])
ggplot(maize, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC2.2v3, y=transformed_PC3.2v3, color = Regions), alpha = 0.5)+
  theme_bw()+scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))

### Full maize set without the outliers vs. geography ###

maize.protest<-protest(Y = select(maize_noOutliers, c("PC1","PC2")),
                       X = select(maize_noOutliers, c("locations_longitude","locations_latitude")))

maize_noOutliers<-maize_noOutliers %>% add_column(as.data.frame(maize.protest$Yrot))
colnames(maize_noOutliers)[67:68]<-c("rotated_PC1","rotated_PC2")
maize_noOutliers<-maize_noOutliers %>%mutate(rotated_PC1 = rotated_PC1+maize.protest$translation[,1],
                                             rotated_PC2 = rotated_PC2+maize.protest$translation[,2])
maize_noOutliers<-mutate(maize_noOutliers, Regions = case_when(countries_country_name %in% c("BOLIVIA","CHILE","PERU") ~ "AndeanHighland",
                                                               countries_country_name %in% c("BRAZIL") ~ "Amazonia",
                                                               countries_country_name %in% c("COLOMBIA","ECUADOR","FRENCH GUIANA","HONDURAS","SURINAME","VENEZUELA") ~ "CentralSouthAmerica",
                                                               countries_country_name %in% c("PARAGUAY","URUGUAY") ~ "ChacoAmerindian",
                                                               countries_country_name %in% c("MEXICO") ~ "Mexico",
                                                               countries_country_name %in% c("ARGENTINA") ~ "Patagonia",
                                                               countries_country_name %in% c("COSTA RICA","EL SALVADOR","GUATEMALA",
                                                                                             "NICARAGUA","PANAMA") ~ "Mesoamerica",
                                                               countries_country_name %in% c("ANTIGUA AND BARBUDA","BARBADOS","CUBA",
                                                                                             "DOMINICAN REPUBLIC","GRENADA","GUADELOUPE",
                                                                                             "HAITI","MARTINIQUE","PUERTO RICO","SAINT VINCENT AND THE GRENADINES",
                                                                                             "TRINIDAD AND TOBAGO","VIRGIN ISLANDS (BRITISH)","VIRGIN ISLANDS (U.S.)") ~ "Caribbean"))

maize_noOutliers%>% 
  ggplot(aes(x=rotated_PC1, y=rotated_PC2, color = Regions))+ 
  geom_point(aes(x=locations_longitude,y=locations_latitude), shape = 3, color = "black")+
  geom_point(alpha = 0.5)+
  theme_bw()+
  xlab("Longitude")+ylab("Latitude")+
  scale_color_manual(name = "Regions by Country",
                     values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                     labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/maize_PConGeography.png", device = "png",dpi=300,
       height = 6, width=7.5)

### Maize vs. geography mexico only ###
Mexmaize_procrustes<-procrustes(Y = select(Mexmaize, c("PC1","PC2")),
                                X = select(Mexmaize, c("longitude","latitude"))
)
plot(Mexmaize_procrustes)

Mexmaize.protest<-protest(Y = select(Mexmaize, c("PC1","PC2")),
                          X = select(Mexmaize, c("longitude","latitude")))

ggplot(Mexmaize, aes(x=PC1, y=PC2))+
  geom_point(aes(color = locations_elevation, alpha = longitude))+
  theme_bw()+
  scale_color_gradient(low="orange",high = "blue")+
  ggtitle("Mexican maize PCs by elevation")

####Transformed maize procrustes vs. geography####
#### Transformed procrustes####

transformed_maize.procrustes <- procrustes(X = select(maize_admixRemoved, c(locations_longitude, locations_latitude)),
                                           Y = select(maize_admixRemoved, c(PC1, PC2)))
transformed_maize.translation<-matrix(rep(transformed_maize.procrustes$translation,2909),nrow=2909,ncol=2,byrow =TRUE) 
maize_admixRemoved<-add_column(maize_admixRemoved, 
                               transformed_PC1 = (transformed_maize.procrustes$scale*t((transformed_maize.procrustes$rotation) %*% t(as.matrix(maize_admixRemoved[,c("PC1","PC2")])))+transformed_maize.translation)[,1],
                               transformed_PC2 = (transformed_maize.procrustes$scale*t((transformed_maize.procrustes$rotation) %*% t(as.matrix(maize_admixRemoved[,c("PC1","PC2")])))+transformed_maize.translation)[,2])
maize_admixRemoved<-mutate(maize_admixRemoved, Regions = case_when(countries_country_name %in% c("BOLIVIA","CHILE","PERU") ~ "AndeanHighland",
                                                                   countries_country_name %in% c("BRAZIL") ~ "Amazonia",
                                                                   countries_country_name %in% c("COLOMBIA","ECUADOR","FRENCH GUIANA","HONDURAS","SURINAME","VENEZUELA") ~ "CentralSouthAmerica",
                                                                   countries_country_name %in% c("PARAGUAY","URUGUAY") ~ "ChacoAmerindian",
                                                                   countries_country_name %in% c("MEXICO") ~ "Mexico",
                                                                   countries_country_name %in% c("ARGENTINA") ~ "Patagonia",
                                                                   countries_country_name %in% c("COSTA RICA","EL SALVADOR","GUATEMALA",
                                                                                                 "NICARAGUA","PANAMA") ~ "Mesoamerica",
                                                                   countries_country_name %in% c("ANTIGUA AND BARBUDA","BARBADOS","CUBA",
                                                                                                 "DOMINICAN REPUBLIC","GRENADA","GUADELOUPE",
                                                                                                 "HAITI","MARTINIQUE","PUERTO RICO","SAINT VINCENT AND THE GRENADINES",
                                                                                                 "TRINIDAD AND TOBAGO","VIRGIN ISLANDS (BRITISH)","VIRGIN ISLANDS (U.S.)") ~ "Caribbean"))

ggplot(maize_admixRemoved, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2, color = Regions), alpha = 0.5)+
  theme_bw()+scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  xlab("Longitude")+ylab("Latitude")


ggplot(maize_admixRemoved, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2, color = locations_elevation))+
  theme_bw()+scale_color_viridis_c()

ggplot(test, aes(x=PC1,y=PC2))+
  geom_point(aes(color = admix))+
  scale_color_viridis_c()+
  theme_bw()

protest(X = select(maize_admixRemoved, locations_latitude, locations_longitude),
        Y = select(maize_admixRemoved, c(PC1, PC2)))
protest(X = select(maize_admixRemoved, locations_latitude, locations_longitude),
        Y = select(maize_admixRemoved, PC2, PC3))

ggplot(maize_admixRemoved, aes(x=PC2, y=PC3, color = countries_country_name))+
  geom_point( alpha = 0.5)+
  theme_bw()

transformed_maize.2v3.procrustes <- procrustes(X = select(maize_admixRemoved, c(locations_longitude, locations_latitude)),
                                               Y = select(maize_admixRemoved, c(PC2, PC3)))
transformed_maize.2v3.translation<-matrix(rep(transformed_maize.procrustes$translation,2909),nrow=2909,ncol=2,byrow =TRUE) 
maize_admixRemoved<-add_column(maize_admixRemoved, 
                               transformed_PC2.2v3 = (transformed_maize.2v3.procrustes$scale*t((transformed_maize.2v3.procrustes$rotation) %*% t(as.matrix(maize_admixRemoved[,c("PC2","PC3")])))+transformed_maize.2v3.translation)[,1],
                               transformed_PC3.2v3 = (transformed_maize.2v3.procrustes$scale*t((transformed_maize.2v3.procrustes$rotation) %*% t(as.matrix(maize_admixRemoved[,c("PC2","PC3")])))+transformed_maize.2v3.translation)[,2])

ggplot(maize_admixRemoved, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC2.2v3, y=transformed_PC3.2v3, color = Regions), alpha = 0.5)+
  theme_bw()+scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))

maize_admixRemoved<-add_column(maize_admixRemoved, 
                               rotated_PC1 = transformed_maize.procrustes$Yrot[,1],
                               rotated_PC2 = transformed_maize.procrustes$Yrot[,2])

test_axes<- tibble(x1 = c(x_prime.1[1,1],y_prime.1[1,1]), 
                   x2 = c(x_prime.2[1,1],y_prime.2[1,1]), 
                   y1 = c(x_prime.1[2,1],y_prime.1[2,1]), 
                   y2 = c(x_prime.2[2,1],y_prime.2[2,1]))

#maize_admixremoved_rot<-
ggplot(maize_admixRemoved, aes(x=rotated_PC1, y=rotated_PC2))+
  geom_point(aes(color = Regions), alpha = 0.5)+
  geom_abline(slope = rotated_axes_slopes(maize_admixRemoved.procrustes$rotation)[1])+
  geom_abline(slope = rotated_axes_slopes(maize_admixRemoved.procrustes$rotation)[2])+
  geom_segment(data= test_axes*10, aes(x=x1, xend=x2, y=y1, yend=y2),color =c("red","yellow"))+
  
  theme_bw()+scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave(plot=maize_admixremoved_rot, filename="/group/jrigrp11/snodgras_maizePopulations/Plots/maize_admixRemoved_rot_PCAonGeography.png",device = "png",dpi=300, width = 6, height = 6)

maize_admixremoved_pca.plot<-ggplot(maize_admixRemoved, aes(x=PC1, y=PC2))+
  geom_vline(xintercept = 0)+geom_hline(yintercept = 0)+
  geom_point(aes(color = Regions), alpha = 0.5)+
  theme_bw()+scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave(plot=maize_admixremoved_pca.plot, filename="/group/jrigrp11/snodgras_maizePopulations/Plots/maize_admixRemoved_PCAonGeography.png",device = "png",dpi=300, width = 6, height = 6)

###Mexican maize procrustes###
Mexmaize_procrustes<-procrustes(X = select(Mexmaize, c(locations_longitude, locations_latitude)),
                                Y = select(Mexmaize, c(PC1, PC2)))
Mexmaize.translation<-matrix(rep(Mexmaize_procrustes$translation,1727),nrow=1727,ncol=2,byrow =TRUE) 
Mexmaize<-add_column(Mexmaize, 
                     transformed_PC1 = (Mexmaize_procrustes$scale*t((Mexmaize_procrustes$rotation) %*% t(as.matrix(Mexmaize[,c("PC1","PC2")])))+Mexmaize.translation)[,1],
                     transformed_PC2 = (Mexmaize_procrustes$scale*t((Mexmaize_procrustes$rotation) %*% t(as.matrix(Mexmaize[,c("PC1","PC2")])))+Mexmaize.translation)[,2])
ggplot(Mexmaize, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2, color = locations_elevation), alpha = 0.5)+
  theme_bw()+scale_color_viridis_c(name = "Elevation (m)")+xlab("Longitude")+ylab("Latitude")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/Mexmaize_non-transformed_procrustes_withGeography.png",
       dpi=300, device = "png", width = 8, height = 6)

protest(X = select(Mexmaize, locations_latitude, locations_longitude),
        Y = select(Mexmaize, c(PC1, PC2)))

Mexmaize.2v3_procrustes<-procrustes(X = select(Mexmaize, c(locations_longitude, locations_latitude)),
                                    Y = select(Mexmaize, c(PC2, PC3)))
Mexmaize.2v3.translation<-matrix(rep(Mexmaize.2v3_procrustes$translation,1727),nrow=1727,ncol=2,byrow =TRUE) 
Mexmaize<-add_column(Mexmaize, 
                     transformed_PC2.2v3 = (Mexmaize.2v3_procrustes$scale*t((Mexmaize.2v3_procrustes$rotation) %*% t(as.matrix(Mexmaize[,c("PC2","PC3")])))+Mexmaize.2v3.translation)[,1],
                     transformed_PC3.2v3 = (Mexmaize.2v3_procrustes$scale*t((Mexmaize.2v3_procrustes$rotation) %*% t(as.matrix(Mexmaize[,c("PC2","PC3")])))+Mexmaize.2v3.translation)[,2])
ggplot(Mexmaize, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC2.2v3, y=transformed_PC3.2v3, color = locations_elevation), alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()
protest(X = select(Mexmaize, locations_latitude, locations_longitude),
        Y = select(Mexmaize, c(PC2, PC3)))

###Mexican maize transformed procrustes###
Mexmaize_admixRemoved_pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/Mex_subset/Mex_transformed_LDprune_MAF0.01.pca")
Mexmaize_admixRemoved_eigenval<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/Mex_subset/Mex_transformed_LDprune_MAF0.01.eigens") 

#join meta and pcs
Mexmaize_admixRemoved<-mutate(Mexmaize_admixRemoved_pca, sample_id = str_split(ind_id, "\\.",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

ggplot(Mexmaize_admixRemoved, aes(x=PC1, y=PC2))+
  geom_point(aes(color = locations_elevation), alpha=0.5)+
  theme_bw()+scale_color_viridis_c(name = "Elevation (m)")+
  xlab(paste0("PC1 (", signif(Mexmaize_admixRemoved_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(Mexmaize_admixRemoved_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Mex maize with mexicana anc. removed PCA")+
  theme(text = element_text(size = 18))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/Mexmaize_transformed_PCA_PC1_PC2_byElevation.png", 
       device = "png", dpi = 300, width = 9, height = 6)
ggplot(Mexmaize_admixRemoved, aes(x=PC2, y=PC3))+
  geom_point(aes(color = locations_elevation), alpha=0.5)+
  theme_bw()+scale_color_viridis_c(name = "Elevation (m)")+
  xlab(paste0("PC2 (", signif(Mexmaize_admixRemoved_eigenval[2,2]*100,3),"%)"))+
  ylab(paste0("PC3 (", signif(Mexmaize_admixRemoved_eigenval[3,2]*100,3),"%)"))+
  ggtitle("Mex maize with mexicana anc. removed PCA")+
  theme(text = element_text(size = 20))

Mexmaize_admixRemoved_procrustes<-procrustes(X = select(Mexmaize_admixRemoved, c(locations_longitude, locations_latitude)),
                                             Y = select(Mexmaize_admixRemoved, c(PC1, PC2)))
Mexmaize_admixRemoved.translation<-matrix(rep(Mexmaize_admixRemoved_procrustes$translation,1727),nrow=1727,ncol=2,byrow =TRUE) 
Mexmaize_admixRemoved<-add_column(Mexmaize_admixRemoved, 
                                  transformed_PC1 = (Mexmaize_admixRemoved_procrustes$scale*t((Mexmaize_admixRemoved_procrustes$rotation) %*% t(as.matrix(Mexmaize_admixRemoved[,c("PC1","PC2")])))+Mexmaize_admixRemoved.translation)[,1],
                                  transformed_PC2 = (Mexmaize_admixRemoved_procrustes$scale*t((Mexmaize_admixRemoved_procrustes$rotation) %*% t(as.matrix(Mexmaize_admixRemoved[,c("PC1","PC2")])))+Mexmaize_admixRemoved.translation)[,2])
ggplot(Mexmaize_admixRemoved, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2, color = locations_elevation), alpha = 0.5)+
  theme_bw()+scale_color_viridis_c(name = "Elevation (m)")+xlab("Longitude")+ylab("Latitude")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/Mexmaize_transformed_procrustes_withGeography.png",
       dpi = 300, device = "png", width = 8, height = 6)
protest(X = select(Mexmaize_admixRemoved, locations_latitude, locations_longitude),
        Y = select(Mexmaize_admixRemoved, c(PC1, PC2)))

Mexmaize_admixRemoved.2v3_procrustes<-procrustes(X = select(Mexmaize_admixRemoved, c(locations_longitude, locations_latitude)),
                                                 Y = select(Mexmaize_admixRemoved, c(PC2, PC3)))
Mexmaize_admixRemoved.2v3.translation<-matrix(rep(Mexmaize_admixRemoved.2v3_procrustes$translation,1727),nrow=1727,ncol=2,byrow =TRUE) 
Mexmaize_admixRemoved<-add_column(Mexmaize_admixRemoved, 
                                  transformed_PC2.2v3 = (Mexmaize_admixRemoved.2v3_procrustes$scale*t((Mexmaize_admixRemoved.2v3_procrustes$rotation) %*% t(as.matrix(Mexmaize_admixRemoved[,c("PC2","PC3")])))+Mexmaize_admixRemoved.2v3.translation)[,1],
                                  transformed_PC3.2v3 = (Mexmaize_admixRemoved.2v3_procrustes$scale*t((Mexmaize_admixRemoved.2v3_procrustes$rotation) %*% t(as.matrix(Mexmaize_admixRemoved[,c("PC2","PC3")])))+Mexmaize_admixRemoved.2v3.translation)[,2])
ggplot(Mexmaize_admixRemoved, aes(x=locations_longitude,y=locations_latitude))+
  geom_point(color = "black", shape = 3)+
  geom_point(aes(x=transformed_PC2.2v3, y=transformed_PC3.2v3, color = locations_elevation), alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()
protest(X = select(Mexmaize_admixRemoved, locations_latitude, locations_longitude),
        Y = select(Mexmaize_admixRemoved, c(PC2, PC3)))

####FOR SOME REASON THE SCALE IS WAY OFF FOR THE MAIZE PCS COMPARED TO HUMAN PCS####
#procrustes of the pair to get the maize (centroid) PCs onto the human PC coordinates
centroidPCs_ind_human.procrustes<-procrustes(X=select(ungroup(centroidPCs_ind_human), PC_1,PC_2),
                                       Y=select(ungroup(centroidPCs_ind_human), PC1.maize,PC2.maize))
#ss=0.7044, correlation=0.5437, significance=0.001
#test.protest<-protest(Y=select(ungroup(centroidPCs_ind_human), PC_1,PC_2),
#                      X=select(ungroup(centroidPCs_ind_human), PC1.maize,PC2.maize))
#number of rows in centroidPCs_ind_human
centroid.translation<-matrix(rep(centroidPCs_ind_human.procrustes$translation,345),nrow=345,ncol=2,byrow =TRUE)

centroidPCs_ind_human<-centroidPCs_ind_human %>% add_column(as.data.frame(centroidPCs_ind_human.procrustes$Yrot))
colnames(centroidPCs_ind_human)[445:446]<-c("rotated_PC1","rotated_PC2")
centroidPCs_ind_human<-centroidPCs_ind_human %>% ungroup() %>%
  mutate(transformed_rotated_PC1 = (centroidPCs_ind_human.procrustes$scale*t((centroidPCs_ind_human.procrustes$rotation) %*% t(as.matrix(centroidPCs_ind_human[,c("PC1.maize","PC2.maize")])))+centroid.translation)[,1],
         transformed_rotated_PC2 = (centroidPCs_ind_human.procrustes$scale*t((centroidPCs_ind_human.procrustes$rotation) %*% t(as.matrix(centroidPCs_ind_human[,c("PC1.maize","PC2.maize")])))+centroid.translation)[,2])

#plot the maize (centroid, procrustes) PCs and human PCs on the human PC coordinates
centroids_human_maize<-centroidPCs_ind_human %>% group_by(Region.x) %>%
  summarize(mean_PC1.maize = mean(transformed_rotated_PC1, na.rm=T),
            mean_PC2.maize = mean(transformed_rotated_PC2, na.rm=T),
            mean_PC1.human = mean(PC_1, na.rm=T),
            mean_PC2.human = mean(PC_2, na.rm=T))

ggplot(centroidPCs_ind_human )+
  geom_point(shape=17, aes(x=PC_1,y=PC_2, color = Region.x))+
  geom_point(aes(x=transformed_rotated_PC1,y=transformed_rotated_PC2, color = Region.x), alpha = 0.3)+
  geom_segment(data = na.omit(centroids_human_maize),aes(x=mean_PC1.human,y=mean_PC2.human,xend=mean_PC1.maize,yend=mean_PC2.maize))+
  theme_bw()+xlab("human PC1")+ylab("human PC2")+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1aff",
               "#b100ea","#386955","#a7c957"),
    name = "Region", labels = c(
      "Mexico, Center","Mexico, Gulf", 
      "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
      "Mexico, Oaxaca","Mexico, West"))+
  ggtitle("Humans Emphasized")

ggplot(centroidPCs_ind_human )+
  geom_point(shape=17, aes(x=PC_1,y=PC_2, color = Region.x), alpha = 0.3)+
  geom_point(aes(x=transformed_rotated_PC1,y=transformed_rotated_PC2, color = Region.x))+
  geom_segment(data = na.omit(centroids_human_maize),aes(x=mean_PC1.human,y=mean_PC2.human,xend=mean_PC1.maize,yend=mean_PC2.maize))+
  theme_bw()+xlab("human PC1")+ylab("human PC2")+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1aff",
               "#b100ea","#386955","#a7c957"),
    name = "Region", labels = c(
      "Mexico, Center","Mexico, Gulf", 
      "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
      "Mexico, Oaxaca","Mexico, West"))+
  ggtitle("Maize (20km) Emphasized")


#### 6. Correlation of PCs ####

#ind_human_maize_20km_buffer = humans at the center of 21 km radius regions

ind_human_maize_20km_buffer<-select(human, starts_with("PC_"),"sample_name") %>% 
  left_join(x=ind_human_maize_20km_buffer, y= . , by = c("sample_id.human"="sample_name"))

ind_human_maize_20km_buffer<-left_join(x=ind_human_maize_20km_buffer,
          y=maize_pca, 
          by=c("Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize"="ind_id"))

centroidPCs_ind_human<-ind_human_maize_20km_buffer %>% 
  group_by(sample_id.human,Region.human,Ethnicity.human,latitude.human,longitude.human) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)),
            mean_elevation.maize=mean(locations_elevation.maize))


colnames(centroidPCs_ind_human)[6:15]<-paste0("PC",1:10,".human")
colnames(centroidPCs_ind_human)[16:(ncol(centroidPCs_ind_human)-1)]<-paste0("PC",1:2909,".maize")

## All by All Correlation ## 
cor(centroidPCs_ind_human$PC1.human,centroidPCs_ind_human$PC1.maize)
cor(centroidPCs_ind_human[,"PC1.human"],centroidPCs_ind_human[,"PC1.maize"]) 
all_by_all_correlations<-cor(centroidPCs_ind_human[,6:15],centroidPCs_ind_human[,16:(ncol(centroidPCs_ind_human)-1)])
all_by_all_correlations<-all_by_all_correlations %>% as_tibble() %>% rownames_to_column(var = "human_PC")
all_by_all_correlations<-pivot_longer(data = all_by_all_correlations, cols = ends_with(".maize"), names_to = "maize_PC", values_to = "correlation")
all_by_all_correlations<-mutate(all_by_all_correlations, maize_PC = str_remove(maize_PC,".maize"))
all_by_all_correlations<-mutate(all_by_all_correlations, maize_PC = str_remove(maize_PC,"PC"))

all_by_all_correlations$human_PC<-factor(all_by_all_correlations$human_PC, levels = 1:10) 
all_by_all_correlations$maize_PC<-factor(all_by_all_correlations$maize_PC, levels = 1:2909) 

ggplot(all_by_all_correlations)+
  geom_tile(aes(x=human_PC, y=maize_PC, fill = correlation))+
  theme_bw()+theme(axis.text.y = element_blank())+xlab("Human PCs")+ylab("Maize PCs")+
  colorspace::scale_fill_continuous_diverging()

filter(all_by_all_correlations, maize_PC %in% (1:10)) %>%
  ggplot(aes(x=human_PC, y=maize_PC))+
  geom_tile(aes( fill = correlation))+
  geom_text(aes(label = signif(correlation, digits = 3)),size=12)+
  theme_bw()+xlab("Human PCs")+ylab("Maize PCs")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size = 20, face="bold"))+
  colorspace::scale_fill_continuous_diverging()

###mexican only and mexicana admix removed ###
human_21km_pairs<-read_csv("/group/jrigrp11/snodgras_maizePopulations/Individual-human-20km-maize-samples-joined.csv")
colnames(human_21km_pairs)<-c("sample_id.human","Ethnicity.human","Country.human","Region.human","latitude.human","longitude.human","Indigenous_Proportion.human","id.maize","general_identifier.maize",
                                         "name.maize","bank_number.maize","taxonomy_id.maize","collnumb.maize","colldate.maize","location_id.maize","created_on.maize","locations_id.maize","locations_region.maize",
                                         "locations_site_name.maize","locations_elevation.maize","locations_latitude.maize","locations_longitude.maize","countries_id.maize","countries_country_code2.maize",
                                         "countries_country_code3.maize","countries_country_name.maize","taxonomies_id.maize","taxonomies_genus.maize","taxonomies_species.maize","taxonomies_crop_name.maize",
                                         "taxonomies_ploidy.maize","Sample_ID_of_DNA_from_composite_samples.maize","Sample_ID_of_DNA_from_most_recent_CML_regenerations.maize",
                                         "Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize","Tester_GID.maize","Tester_pedigree.maize","Testcross_GID.maize","PrimaryRace.maize","PrimaryPurity.maize",
                                         "SecondaryRace.maize","Pedigree.maize","GrainType1.maize","GrainType2.maize","GrainType3.maize","GrainColor1.maize","GrainColor2.maize","GrainColor3.maize","PopulationType.maize")

human_21km_pairs_PCs<-select(Mexhuman, starts_with("PC_"),"sample_name") %>% 
  left_join(x=human_21km_pairs, y= . , by = c("sample_id.human"="sample_name"))

human_21km_pairs_PCs<-left_join(x=human_21km_pairs_PCs,
                                       y=Mexmaize_admixRemoved, 
                                       by=c("Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize"="ind_id"))

centroidPCs_Mexhuman21km<-human_21km_pairs_PCs %>% 
  group_by(sample_id.human,Region.human,Ethnicity.human,latitude.human,longitude.human) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)),
            mean_elevation.maize=mean(locations_elevation.maize))
colnames(centroidPCs_Mexhuman21km)[6:15]<-paste0("PC",1:10,".human")
colnames(centroidPCs_Mexhuman21km)[16:(ncol(centroidPCs_Mexhuman21km)-1)]<-paste0("PC",1:1731,".maize")

Mexhuman_Mexmaize_procrustes<-procrustes(X=select(ungroup(centroidPCs_Mexhuman21km), PC1.human, PC2.human),
                                         Y=select(ungroup(centroidPCs_Mexhuman21km), PC1.maize, PC2.maize))
protest(X=select(ungroup(centroidPCs_Mexhuman21km), PC1.human, PC2.human),
           Y=select(ungroup(centroidPCs_Mexhuman21km), PC1.maize, PC2.maize))
Mexhuman_Mexmaize.translation<-matrix(rep(Mexhuman_Mexmaize_procrustes$translation,345),nrow=345,ncol=2,byrow =TRUE) 

centroidPCs_Mexhuman21km<-ungroup(centroidPCs_Mexhuman21km)

centroidPCs_Mexhuman21km<-centroidPCs_Mexhuman21km %>% add_column(as.data.frame(Mexhuman_Mexmaize_procrustes$Yrot))
colnames(centroidPCs_Mexhuman21km)[1748:1749]<-c("rotated_PC1","rotated_PC2")
centroidPCs_Mexhuman21km<-centroidPCs_Mexhuman21km %>%
  mutate(transformed_PC1 = (Mexhuman_Mexmaize_procrustes$scale*t((Mexhuman_Mexmaize_procrustes$rotation) %*% t(as.matrix(centroidPCs_Mexhuman21km[,c("PC1.maize","PC2.maize")])))+Mexhuman_Mexmaize.translation)[,1],
         transformed_PC2 = (Mexhuman_Mexmaize_procrustes$scale*t((Mexhuman_Mexmaize_procrustes$rotation) %*% t(as.matrix(centroidPCs_Mexhuman21km[,c("PC1.maize","PC2.maize")])))+Mexhuman_Mexmaize.translation)[,2])

Mexhuman_mexmaize_procrustes_plot<-ggplot(centroidPCs_Mexhuman21km, aes(color = Region.human))+
  geom_point(aes(x=PC1.human, y=PC2.human), shape = 3)+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2), alpha =0.5)+
  theme_bw()+xlab("human PC1")+ylab("human PC2")+
  scale_color_manual(
    values = c("Mexico-Center_of_Mexico"="#fca207",
               "Mexico-Gulf_of_Mexico"="#cb4d8e",
               "Mexico-Mayan_region"="#268189",
               "Mexico-North_of_Mesoamerica"="#2d1aff",
               "Mexico-North_of_Mexico"="#b100ea",
               "Mexico-Oaxaca"="#386955",
               "Mexico-West_of_Mexico"="#a7c957"),
    name = "Region", labels = c(
      "Mexico, Center","Mexico, Gulf", 
      "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
      "Mexico, Oaxaca","Mexico, West"))
ggsave(plot = Mexhuman_mexmaize_procrustes_plot, 
       filename = "/group/jrigrp11/snodgras_maizePopulations/Plots/21km_Mexhuman_Mexmaize_procrustes.png", dpi = 300, 
       device = "png", width = 5, height = 3)  
  

Mex21_all_by_all_correlations<-cor(centroidPCs_Mexhuman21km[,6:15],centroidPCs_Mexhuman21km[,16:(ncol(centroidPCs_Mexhuman21km)-1)])
Mex21_all_by_all_correlations<-Mex21_all_by_all_correlations %>% as_tibble() %>% rownames_to_column(var = "human_PC")
Mex21_all_by_all_correlations<-pivot_longer(data = Mex21_all_by_all_correlations, cols = ends_with(".maize"), names_to = "maize_PC", values_to = "correlation")
Mex21_all_by_all_correlations<-mutate(Mex21_all_by_all_correlations, maize_PC = str_remove(maize_PC,".maize") %>% str_remove("PC"))
Mex21_all_by_all_correlations$human_PC <- factor(Mex21_all_by_all_correlations$human_PC, levels = 1:10)
Mex21_all_by_all_correlations$maize_PC <- factor(Mex21_all_by_all_correlations$maize_PC, levels = 1:1731)

ggplot(Mex21_all_by_all_correlations)+
  geom_tile(aes(x=human_PC, y=maize_PC, fill = correlation))+
  theme_bw()+theme(axis.text.y = element_blank())+xlab("Mexican Human PCs")+ylab("Mexican Maize Admix Removed PCs")+
  colorspace::scale_fill_continuous_diverging()
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_MexicanSamples_PCcorrelationheatmap.all_v_all.png",
       device = "png", dpi = 300, width = 11,height = 8)

filter(Mex21_all_by_all_correlations, maize_PC %in% (1:10)) %>%
  ggplot(aes(x=human_PC, y=maize_PC))+
  geom_tile(aes( fill = correlation))+
  geom_text(aes(label = signif(correlation, digits = 3)),size=5)+
  theme_bw()+xlab("Mexican Human PCs")+ylab("Mexican Maize Admix Removed PCs")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size = 20, face="bold"))+
  colorspace::scale_fill_continuous_diverging()
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_MexicanSamples_PCcorrelationheatmap.PCs1_10.png",
       device = "png", dpi = 300, width = 11,height = 8)

filter(Mex21_all_by_all_correlations, correlation >= 0.5 | correlation <= -0.5) # 427 PC pairs met this criteria

sum(Mexmaize_admixRemoved_eigenval$value[1:731]) #number of maize pcs to get to 0.5
sum(Mexmaize_admixRemoved_eigenval$value[1:50]) #this only gets to the first 5%

filter(Mex21_all_by_all_correlations, maize_PC %in% 1:731) %>% 
  filter(correlation >= 0.5 | correlation <= -0.5) #254 PC pairs meet these criteria

filter(Mex21_all_by_all_correlations, maize_PC %in% 1:20) %>% 
  filter(correlation >= 0.5 | correlation <= -0.5) #6 PC pairs meet these criteria

ggplot(centroidPCs_Mexhuman21km)+
  geom_point(aes(x=PC2.human, y=PC7.maize, color = Region.human))+
  theme_bw()

plotting_pairs<-filter(Mex21_all_by_all_correlations, maize_PC %in% 1:20) %>% 
  filter(correlation >= 0.5 | correlation <= -0.5)
pdf("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_MexicanSamples_PCcorrelationPoints.pdf")

ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human,y=latitude.human, color=Region.human))+
  geom_point()+theme_bw()+ggtitle("Locations of each human 'Region' in Mexico")

for(i in 1:nrow(plotting_pairs)){
  plt<-ggplot(centroidPCs_Mexhuman21km, aes(color=Region.human))+
    geom_point(aes_string(x=paste0("PC",plotting_pairs$human_PC[i],".human"), y=paste0("PC",plotting_pairs$maize_PC[i],".maize")))+
    theme_bw()+
    ggtitle(paste0("Correlation = ",signif(plotting_pairs$correlation[i], 3)))
  print(plt)
}

dev.off()

#what is human PC3 and is any human PC correlated with elevation
cor(centroidPCs_Mexhuman21km[,6:15],centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])
#PC 3 -0.15, PC6 0.13, PC10 -0.213

cor.test(pull(centroidPCs_Mexhuman21km[,8]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])) #PC3, significant p=0.005
cor.test(pull(centroidPCs_Mexhuman21km[,11]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])) #PC6 significant p=0.014
cor.test(pull(centroidPCs_Mexhuman21km[,15]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])) #PC10 significant p=6.5e-5
cor.test(pull(centroidPCs_Mexhuman21km[,6]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])) #PC1 not significant p=0.1257

pdf("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_MexicanSamples_HumanElevationPCcorrPoints.pdf",
    height = 4, width = 6)

ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human, y=latitude.human, color = mean_elevation.maize))+
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()+ggtitle("Locations by elevation")
ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human, y=latitude.human, color = PC3.human))+
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()+ggtitle("Correlation with elevation = -0.15")
ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human, y=latitude.human, color = PC6.human))+
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()+ggtitle("Correlation with elevation = 0.13")
ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human, y=latitude.human, color = PC10.human))+
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()+ggtitle("Correlation with elevation = -0.21")


dev.off()

pdf("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_MexicanSamples_MaizePCsGeoCorrWithHumanPCs.pdf",
    height = 4, width = 6)

ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human, y=latitude.human, color = PC1.maize))+
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()+ggtitle("maize PC1, human PC1/PC3")
ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human, y=latitude.human, color = PC7.maize))+
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()+ggtitle("maize PC7, human PC1/PC2/PC3")
ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human, y=latitude.human, color = PC11.maize))+
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.5)+
  theme_bw()+scale_color_viridis_c()+ggtitle("maize PC11, human PC4")

dev.off()

## after MaizePVEbyHumanPCs.R where each Maize snp is lm'ed with the first 10PCs of humans ##
maize21km_human_lm_Rsq<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/lm_results.tsv")

ggplot(maize21km_human_lm_Rsq)+
  geom_histogram(aes(x=r.squared), fill = "navy", binwidth = 0.01)+
  #geom_histogram(aes(x=adj.r.squared), fill="red", alpha=0.5, binwidth = 0.01)+
  geom_vline(xintercept = 0)+theme_bw()

max(maize21km_human_lm_Rsq$r.squared)
max(maize21km_human_lm_Rsq$adj.r.squared)

filter(maize21km_human_lm_Rsq, adj.r.squared >= 0.3) %>% write_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/snps_atleast_0.3r.sq_lm.tsv")

# on average across all SNPs, how much does human PCs explain maize genotypes
summarize(maize21km_human_lm_Rsq, 
          mean.r.squared = mean(r.squared, na.rm = TRUE),
          mean.adj.r.squared = mean(adj.r.squared, na.rm = TRUE),
          sd.r.squared = sd(r.squared, na.rm = TRUE),
          sd.adj.r.squared = sd(adj.r.squared, na.rm = TRUE))
c(0.0310 - (3*0.0203), 0.0310 + (3*0.0203))
c(0.00715 - (3*0.0208), 0.00715 + (3*0.0208))

# so on average human pcs explain ~1-3% of the maize genotypes?
#Is this comparable to what we'd get if we ran this on the maize PCs? 
#ball park should be close to:
Mexmaize_eigenval[1:10,2] %>% sum() #0.03418

#map genotypes of top SNPs onto geographic positions
#in unix:
#grep "^#CHROM" LDprune_MAF0.01_21km_MaizeGBS.tsv > snps_atleast_0.3r.sq_lm.genotypes.tsv
# cut -f 1 snps_atleast_0.3r.sq_lm.tsv | grep -f - LDprune_MAF0.01_21km_MaizeGBS.tsv >> snps_atleast_0.3r.sq_lm.genotypes.tsv

highLM_snps<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/snps_atleast_0.3r.sq_lm.genotypes.tsv")
highLM_snps<-pivot_longer(highLM_snps, cols = starts_with("SEEDGWAS"), names_to = "ind_id", values_to = "genotypes")
highLM_snps<-mutate(highLM_snps, sample_id = str_split(ind_id, ":", simplify=T)[,1])

highLM_snps<-left_join(human_21km_pairs_PCs, highLM_snps, by =c("Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize" = "sample_id")) %>% 
  select('genotypes',"ID","Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize", ends_with(".human"), contains("elevation.maize"), "locations_latitude.maize","locations_longitude.maize") 

filter(highLM_snps, ID == "S3_97146428") %>% mutate(genotypes = as.character(genotypes)) %>% 
  ggplot(aes(x=locations_longitude.maize, y=locations_latitude.maize, color = genotypes, shape = Region.human))+
  geom_point(alpha = 0.6)+ #scale_color_manual(values = c(NA = "#555555", 0 = "#b30000", 1 = "#600080", 2 = "#003399"))+
  theme_bw()+ggtitle("S3_97146428, r.sq = 0.334")

snps_of_interest<-filter(maize21km_human_lm_Rsq, adj.r.squared >= 0.3)

pdf("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_MexicanSamples_highSNPlmGenotypes.pdf")

for(i in 1:nrow(snps_of_interest)){
  
  plt<-filter(highLM_snps, ID == snps_of_interest$snp_ID[i]) %>% mutate(genotypes = as.character(genotypes)) %>% 
    ggplot(aes(x=locations_longitude.maize, y=locations_latitude.maize, color = genotypes, shape = Region.human))+
    geom_point(alpha = 0.8)+ #scale_color_manual(values = c(NA = "#555555", 0 = "#b30000", 1 = "#600080", 2 = "#003399"))+
    theme_bw()+ggtitle(paste("SNP",snps_of_interest$snp_ID[i],", R.sq =",snps_of_interest$r.squared[i]))
  print(plt)
}

dev.off()

maize21km_human_lm_Rsq<-maize21km_human_lm_Rsq %>% 
  mutate(Chromosome = str_split(snp_ID, "_", simplify = T)[,1] %>% str_remove_all("S"),
         bp = str_split(snp_ID, "_", simplify = T)[,2] %>% as.numeric())
#code from https://r-graph-gallery.com/101_Manhattan_plot.html
manhattan_plot_rsq<- maize21km_human_lm_Rsq %>% 
  group_by(Chromosome) %>% summarize(chr_len = max(bp)) %>% # compute chromosome size
  mutate(tot = cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% # calculate cumulative position of each chr
  left_join(maize21km_human_lm_Rsq, ., by= c("Chromosome"="Chromosome")) %>% #add to initial dataset
  arrange(Chromosome, bp) %>% mutate(BPcum = bp+tot) #add cumulative position of each SNP
axisdf = manhattan_plot_rsq %>% group_by(Chromosome) %>% summarize(center = (max(BPcum) + min(BPcum) ) / 2)

manhattan_plot_rsq$Chromosome<-factor(manhattan_plot_rsq$Chromosome, levels = 
                                        c("1","2","3","4","5","6","7","8","9","10","0"))

manhattan_plot_rsq.plot<-ggplot(manhattan_plot_rsq, aes(x=BPcum, y= r.squared)) +
  geom_point(aes(color = Chromosome, alpha = 0.8, size = 1.3))+
  geom_hline(yintercept = quantile(manhattan_plot_rsq$r.squared,prob=1-5/100))+ #line at top 5%
  scale_color_manual(values = c(rep(c("navy","goldenrod"), 5),"darkgrey"))+
  scale_x_continuous(label = axisdf$Chromosome, breaks = axisdf$center)+
  scale_y_continuous(expand = c(0,0) )+
  theme_bw()+theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
ggsave(plot = manhattan_plot_rsq.plot, 
       filename = "/group/jrigrp11/snodgras_maizePopulations/Plots/LMr-squared_manhattan_plot.png",
       device = "png", dpi = 300, width = 6, height = 4)

