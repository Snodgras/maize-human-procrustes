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

### maize Mexican only ###

Mexmaize_pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/Mex_subset/Mex_non_transformed_LDprune_MAF0.01.pca")
Mexmaize_eigenval<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/admix_removal/Mex_subset/Mex_non_transformed_LDprune_MAF.01.eigens") 

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

transformed_21km_maize.pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.pca")
transformed_21km_maize.eigenval<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.eigens")
colnames(transformed_21km_maize.eigenval)<-c("PC","PVE")
transformed_21km_maize<-right_join(x=selected_passport_ids, y=transformed_21km_maize.pca, by=c("gwas_sample_id"="ind_id")) %>% 
  left_join(x=., y=maize_meta, by=c("gwas_sample_id"="Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

ggplot(transformed_21km_maize, aes(x=PC1, y=PC2))+
  geom_point(aes(color=locations_elevation))+theme_bw()+
  scale_color_viridis_c(name="Elevation")+
  xlab(paste0("PC1 (", signif(transformed_21km_maize.eigenval$PVE[1]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(transformed_21km_maize.eigenval$PVE[2]*100,3),"%)"))+
  ggtitle("Transformed 21km radius maize PCA")

#Goal is to find the centroid of PCs for each set of maize samples for a given human sample
centroidPCs_ind_human<-ind_human_maize_20km_buffer %>% 
  group_by(sample_id.human,Region.human,Ethnicity.human,latitude.human,longitude.human) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)),
            mean_elevation.maize=mean(locations_elevation.maize))
colnames(centroidPCs_ind_human)<-c("sample_name","Region","Ethnicity","latitude","longitude",
                                   paste0("PC",1:422,".maize"), "mean_elevation.maize")
#add in the human PCs, specifically the Mexican only PCs
centroidPCs_ind_human<-inner_join(centroidPCs_ind_human, Mexhuman, by=("sample_name"))
centroidPCs_ind_human_Americas<-inner_join(centroidPCs_ind_human, human, by = "sample_name")

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
  ggtitle("Mexican only, Human PCA")+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))#guides(color = "none")

### Maize full data ###
#plot pca
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

ggplot(maize, aes(x=PC1, y=PC2))+
  geom_point(aes(color = Regions), alpha = 0.6)+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(maize_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(maize_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Maize PCA")+scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  theme(text = element_text(size = 20))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/maize_PC1PC2.byRegion.png",device="png",dpi=300)

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
  xlab(paste0("PC2 (", signif(Mexmaize_eigenval[2,2]*100,3),"%)"))+
  ylab(paste0("PC3 (", signif(Mexmaize_eigenval[3,2]*100,3),"%)"))+
  ggtitle("Mexican Maize PCA, PC2 vs PC3")+ scale_color_viridis_c()

maize_Mex_states<-read_csv("/group/jrigrp11/snodgras_maizePopulations/Intersection-MexicanMaize-StatesofMexico.csv")
maize_Mex_states<-maize_Mex_states %>% 
  mutate(Regional_Name = case_when(
    shapeName %in% c("Baja California","Baja California Sur") ~ "Baja California",
    shapeName %in% c("Chihuahua","Durango","Zacatecas","Aguascalientes","Coahuila de Zaragoza") ~ "Mexican Plateau",
    shapeName %in% c("Nuevo Leon","San Luis Potosi") ~ "Sierra Madre Oriental",
    shapeName %in% c("Guanajuato","Hidalgo","Queretaro de Arteaga") ~ "Mesa Central",
    shapeName %in% c("Sinaloa","Sonora","Nayarit") ~ "Sierra Madre Occidental and Pacific Coastal Lowlands",
    shapeName %in% c("Jalisco","Colima","Michoacan de Ocampo","Mexico","Distrito Federal","Morelos","Puebla","Tlaxcala") ~ "Cordillera Neovolcanica",
    shapeName %in% c("Tamaulipas","Veracruz de Ignacio de la Llave","Tabasco") ~ "Gulf Coastal Plain",
    shapeName %in% c("Guerrero","Oaxaca","Chiapas") ~ "Southern Highlands",
    shapeName %in% c("Yucatan","Campeche","Quintana Roo") ~ "Yucatan Peninsula"
  ))
colnames(maize_Mex_states)[1:41]<-colnames(maize_meta)

Mexmaize<-inner_join(Mexmaize, maize_Mex_states,by = c("id", "general_identifier", "name", "bank_number", "taxonomy_id", "collnumb", "colldate", "location_id",
                                                  "locations_id", "locations_region", "locations_site_name", "locations_elevation", "locations_latitude",
                                                  "locations_longitude", "countries_id", "countries_country_code2", "countries_country_code3", "countries_country_name",
                                                  "taxonomies_id", "taxonomies_genus", "taxonomies_species", "taxonomies_crop_name", "taxonomies_ploidy",
                                                  "Sample_ID_of_DNA_from_composite_samples", "Sample_ID_of_DNA_from_most_recent_CML_regenerations", "Tester_GID",
                                                  "Tester_pedigree", "Testcross_GID", "PrimaryRace", "PrimaryPurity", "SecondaryRace", "Pedigree", "GrainType1", "GrainType2", "GrainType3",
                                                  "GrainColor1", "GrainColor2", "GrainColor3", "PopulationType"))
ggplot(Mexmaize, aes(x=PC1, y=PC2))+
  geom_point(aes(color = Regional_Name), alpha = 0.8)+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(Mexmaize_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(Mexmaize_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Mexican Maize PCA, PC1 vs PC2")+ 
scale_color_manual(values = c("Baja California"="#E28AFF",
                             "Mexican Plateau"="#2d1a77",
                             "Sierra Madre Oriental"="#6E54D9",
                             "Sierra Madre Occidental and Pacific Coastal Lowlands"="#b100ea",
                             "Mesa Central"="#fca207",
                             "Gulf Coastal Plain"="#cb4d8e",
                             "Cordillera Neovolcanica"="#a7c957",
                             "Southern Highlands"="#386651",
                             "Yucatan Peninsula"="#268189"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  theme(text = element_text(size = 20))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/Mexmaize_PCA_PC1_PC2_byRegion.png")

#make screeplot
Mexmaize_eigenval[1:50,] %>%
  mutate(name = str_remove(name, pattern = "PC") %>% as.numeric()) %>%
  ggplot(aes(x=name, y=(value*100)))+
  geom_bar(stat = "identity")+
  geom_line(color = "red")+
  geom_point(color = "red")+
  xlab("PC")+ylab("Percent Variance Explained")+
  ggtitle("Mexican Maize PC Scree Plot")+
  theme_bw()
ggsave("/group/jrigrp11/snodgras/snodgras_maizePopulations/Plots/2026-06-05-MexMaizeScreePlot.png",
       device = "png", dpi = 300, 
       width = 5, height = 4, units = "in")

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

maize_admixRemoved<-maize_admixRemoved %>% mutate(Regions = case_when(countries_country_name %in% c("BOLIVIA","CHILE","PERU") ~ "AndeanHighland",
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

ggplot(maize_admixRemoved, aes(x=PC1, y=PC2))+
  geom_point(aes(color = Regions), alpha=0.5)+
  theme_bw()+
  xlab(paste0("PC1 (", signif(maize_admixRemoved_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(maize_admixRemoved_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Maize with mexicana anc. removed PCA")+
  theme(text = element_text(size = 20))+
  scale_color_manual(name = "Regions by Country",
                     values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                     labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/maize_admixtureRemoved.PC1PC2.byRegion.allSNPs.png",
       height = 7, width = 8, dpi=300, device = "png")

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
temp<-inner_join(maize, global_maize_admix, by = "sample_id") %>% 
  select(Admix, PC1, PC2, locations_elevation) 
cor.test(temp$Admix,temp$PC1)
cor.test(temp$Admix,temp$PC2)
cor.test(temp$locations_elevation,temp$PC1)
cor.test(temp$locations_elevation,temp$PC2)
cor.test(temp$Admix, temp$locations_elevation)

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
temp<-maize_nontransformed %>% 
  select(Admix, PC1, PC2, locations_elevation, locations_longitude,locations_latitude) 
cor.test(temp$Admix,temp$PC1)
cor.test(temp$Admix,temp$PC2)
cor.test(temp$locations_elevation,temp$PC1)
cor.test(temp$locations_elevation,temp$PC2)
cor.test(temp$locations_latitude,temp$PC1)
cor.test(temp$locations_latitude,temp$PC2)
cor.test(temp$locations_longitude,temp$PC1)
cor.test(temp$locations_longitude,temp$PC2)

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

#slope = (y2-y1)/(x2-x1)

###rotated procrustes figures ###

#ref
ggplot(human_noOutliers, aes(x=PC_1, y=PC_2, color=Region))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)

axes<-tibble(x=c(-100,0),
             xend=c(100,0),
             y=c(0,-100),
             yend=c(0,100))
#rotated axes
human_rotAxes<-as.matrix(human_procrustes$rotation)%*%as.matrix(axes) %>%
  as_tibble() %>%
  mutate(humanPC = c("PC1","PC2"),
         slope = (y-yend)/(x-xend),
         intercept = human_procrustes$translation[1,2] - slope * (human_procrustes$translation[1,1]-0))
#plot
ggplot()+
  geom_abline(data = human_rotAxes, aes(slope = slope, intercept = intercept, linetype = humanPC))+
   geom_point(data = human_noOutliers, aes(x=transformed_PC1, y=transformed_PC2, color = Region), alpha = 0.5)+
  scale_color_manual(
    values = c("#538fff","#ff5b58","#ecc500","#28d2ab","#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957","#ec00c5"),
    name = "Region", labels = c("Amazonia", "Andean Highland","Central South America",
                                "Chaco Amerindian", "Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West","Patagonia"),
    guide = guide_legend(override.aes = list(alpha = 1) ))+
  theme_bw()+xlab("Longitude")+ylab("Latitude")+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/human_rot_PCAonGeography.png",
       device ="png",dpi=300, width = 6, height = 5)  

### Human Mexican only ###
Mexhuman_procrustes<-procrustes(Y = select(Mexhuman, c("PC_1","PC_2")),
                                     X = select(Mexhuman, c("longitude","latitude"))
)
plot(Mexhuman_procrustes)
protest(Y = select(Mexhuman, c("PC_1","PC_2")),
        X = select(Mexhuman, c("longitude","latitude")))

### Maize full data vs. geography ###
#compute procrustes using vegan 2.6-4

####Original all maize procrustes vs geography####
maize_procrustes<-procrustes(Y = select(maize, c("PC1","PC2")),
                             X = select(maize, c("locations_longitude","locations_latitude"))
)

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

ggplot(maize_admixRemoved, aes(x=rotated_PC1, y=rotated_PC2))+
  geom_point(aes(color = Regions), alpha = 0.5)+
  theme_bw()+scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  xlab("")+ylab("")+ theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
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

Mexmaize_rotAxes<-as.matrix(Mexmaize_procrustes$rotation)%*%as.matrix(axes) %>%
  as_tibble() %>%
  mutate(maizePC = c("PC1","PC2"),
         slope = (y-yend)/(x-xend),
         intercept = Mexmaize_procrustes$translation[1,2] - slope * (Mexmaize_procrustes$translation[1,1]-0))


ggplot(Mexmaize, aes(x=locations_longitude,y=locations_latitude))+
  geom_abline(data = Mexmaize_rotAxes, aes(slope = slope, intercept = intercept, linetype = maizePC))+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2, color = Regional_Name), alpha = 0.5)+
  theme_bw()+xlab("Longitude")+ylab("Latitude")+coord_fixed()+
  scale_color_manual(values = c("Baja California"="#E28AFF",
                                "Mexican Plateau"="#2d1a77",
                                "Sierra Madre Oriental"="#6E54D9",
                                "Sierra Madre Occidental and Pacific Coastal Lowlands"="#b100ea",
                                "Mesa Central"="#fca207",
                                "Gulf Coastal Plain"="#cb4d8e",
                                "Cordillera Neovolcanica"="#a7c957",
                                "Southern Highlands"="#386651",
                                "Yucatan Peninsula"="#268189"))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/Mexmaize_procrustes_raw_latlong_byRegionalName.png", 
       device="png",dpi=300,width = 8, height=7)


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
Mexmaize_admixRemoved<-inner_join(Mexmaize_admixRemoved, 
           select(maize_Mex_states, "Sample_ID_of_DNA_from_single_plants_used_in_GWAS","Regional_Name"),
           by=c("sample_id"="Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

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

ggplot(Mexmaize_admixRemoved, aes(x=PC1, y=PC2))+
  geom_point(aes(color = Regional_Name), alpha=0.8)+
  theme_bw()+
  xlab(paste0("PC1 (", signif(Mexmaize_admixRemoved_eigenval[1,2]*100,3),"%)"))+
  ylab(paste0("PC2 (", signif(Mexmaize_admixRemoved_eigenval[2,2]*100,3),"%)"))+
  ggtitle("Mex maize with mexicana anc. removed PCA")+
  scale_color_manual(values = c("Baja California"="#E28AFF",
                                "Mexican Plateau"="#2d1a77",
                                "Sierra Madre Oriental"="#6E54D9",
                                "Sierra Madre Occidental and Pacific Coastal Lowlands"="#b100ea",
                                "Mesa Central"="#fca207",
                                "Gulf Coastal Plain"="#cb4d8e",
                                "Cordillera Neovolcanica"="#a7c957",
                                "Southern Highlands"="#386651",
                                "Yucatan Peninsula"="#268189"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  theme(text = element_text(size = 20))
  
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/Mexmaize_transformed_PCA_PC1_PC2_byRegions.png", 
       device = "png", dpi = 300, width = 9, height = 6)

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


Mexmaize_admixRemoved_rotAxes<-as.matrix(Mexmaize_admixRemoved_procrustes$rotation)%*%as.matrix(axes) %>%
  as_tibble() %>%
  mutate(maizePC = c("PC1","PC2"),
         slope = (y-yend)/(x-xend),
         intercept = Mexmaize_admixRemoved_procrustes$translation[1,2] - slope * (Mexmaize_admixRemoved_procrustes$translation[1,1]-0))

Mexmaize_admixRemoved_rotAxes<-Mexmaize_admixRemoved_rotAxes %>% 
  mutate(label.x = c(-95,-101),
         label.y = slope*label.x+intercept)
  
ggplot(Mexmaize_admixRemoved, aes(x=locations_longitude,y=locations_latitude))+
  geom_abline(data = Mexmaize_admixRemoved_rotAxes, aes(slope = slope, intercept = intercept, linetype = maizePC))+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2, color = Regional_Name), alpha = 0.5)+
  geom_text(data=Mexmaize_admixRemoved_rotAxes, aes(x = label.x-0.5,
                                                    y=label.y-0.5, 
                                                    label=paste0(signif(Mexmaize_admixRemoved_eigenval$value[1:2]*100,3),"%")))+
  theme_bw()+xlab("Longitude")+ylab("Latitude")+coord_fixed()+
  scale_color_manual(values = c("Baja California"="#CCCCCC",
                               "Mexican Plateau"="#2d1a77",
                               "Sierra Madre Oriental"="#6E54D9",
                               "Sierra Madre Occidental and Pacific Coastal Lowlands"="#b100ea",
                               "Mesa Central"="#fca207",
                               "Gulf Coastal Plain"="#cb4d8e",
                               "Cordillera Neovolcanica"="#a7c957",
                               "Southern Highlands"="#386651",
                               "Yucatan Peninsula"="#268189"))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/Mexmaize_procrustes_transformed_latlong_byRegionalName.png", device="png",dpi=300,width = 8, height=7)

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

#### Maize Human Procrustes ####
#procrustes of the pair to get the maize (centroid) PCs onto the human PC coordinates (Mexhuman)
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

#procrustes of the pair to get the maize (centroid) PCs onto the human PC coordinates (Mexhuman)
centroidPCs_ind_human_Americas.procrustes<-procrustes(X=select(ungroup(centroidPCs_ind_human_Americas), PC_1.y,PC_2.y),
                                             Y=select(ungroup(centroidPCs_ind_human_Americas), PC1.maize,PC2.maize))
protest(X=select(ungroup(centroidPCs_ind_human_Americas), PC_1.y,PC_2.y),
           Y=select(ungroup(centroidPCs_ind_human_Americas), PC1.maize,PC2.maize))

centroid_Americas.translation<-matrix(rep(centroidPCs_ind_human_Americas.procrustes$translation,345),nrow=345,ncol=2,byrow =TRUE)

#centroidPCs_ind_human_Americas<-centroidPCs_ind_human_Americas %>% add_column(as.data.frame(centroidPCs_ind_human_Americas.procrustes$Yrot))
#colnames(centroidPCs_ind_human_Americas)[445:446]<-c("rotated_PC1","rotated_PC2")
centroidPCs_ind_human_Americas<-centroidPCs_ind_human_Americas %>% ungroup() %>%
  mutate(transformed_rotated_PC1 = (centroidPCs_ind_human_Americas.procrustes$scale*t((centroidPCs_ind_human_Americas.procrustes$rotation) %*% t(as.matrix(centroidPCs_ind_human_Americas[,c("PC1.maize","PC2.maize")])))+centroid_Americas.translation)[,1],
         transformed_rotated_PC2 = (centroidPCs_ind_human_Americas.procrustes$scale*t((centroidPCs_ind_human_Americas.procrustes$rotation) %*% t(as.matrix(centroidPCs_ind_human_Americas[,c("PC1.maize","PC2.maize")])))+centroid_Americas.translation)[,2])

Mexmaize.eigens<-read_tsv("lm_results/Mex_non_transformed_LDprune_MAF.01.eigens")
Mexmaize.pca<-read_tsv("lm_results/Mex_non_transformed_LDprune_MAF0.01.pca")
Mexmaize<-mutate(Mexmaize.pca, sample_id = str_split(ind_id, ":",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

Mexmaize_by_states<-read_csv("Intersection-MexicanMaize-StatesofMexico.csv")
Mexmaize<-inner_join(x=Mexmaize, y=select(Mexmaize_by_states, id, shapeName), join_by(id))

Mexmaize_PC2v3_procrustes<-procrustes(Y = select(Mexmaize, c("PC2","PC3")),
                                      X = select(Mexmaize, c("locations_longitude","locations_latitude"))
)
protest(Y = select(Mexmaize, c("PC2","PC3")),
        X = select(Mexmaize, c("locations_longitude","locations_latitude"))
) #r-squared =0.5414


Mexmaize_PC2v3_rotation<-tibble(MX_State = Mexmaize$shapeName,
                                latitude = Mexmaize$locations_latitude,
                                longitude = Mexmaize$locations_longitude,
                                rotated_PC2.maize = Mexmaize_PC2v3_procrustes$Yrot[,1],
                                rotated_PC3.maize = Mexmaize_PC2v3_procrustes$Yrot[,2],
                                PC2.maize = Mexmaize$PC2,
                                PC3.maize = Mexmaize$PC3) 

Mexmaize_PC2v3_rotation<-Mexmaize_PC2v3_rotation %>% mutate(Regional_Name = case_when(
  MX_State %in% c("Baja California","Baja California Sur") ~ "Baja California",
  MX_State %in% c("Chihuahua","Durango","Zacatecas","Aguascalientes","Coahuila de Zaragoza") ~ "Mexican Plateau",
  MX_State %in% c("Nuevo Leon","San Luis Potosi") ~ "Sierra Madre Oriental",
  MX_State %in% c("Guanajuato","Hidalgo","Queretaro de Arteaga") ~ "Mesa Central",
  MX_State %in% c("Sinaloa","Sonora","Nayarit") ~ "Sierra Madre Occidental and Pacific Coastal Lowlands",
  MX_State %in% c("Jalisco","Colima","Michoacan de Ocampo","Mexico","Distrito Federal","Morelos","Puebla","Tlaxcala") ~ "Cordillera Neovolcanica",
  MX_State %in% c("Tamaulipas","Veracruz de Ignacio de la Llave","Tabasco") ~ "Gulf Coastal Plain",
  MX_State %in% c("Guerrero","Oaxaca","Chiapas") ~ "Southern Highlands",
  MX_State %in% c("Yucatan","Campeche","Quintana Roo") ~ "Yucatan Peninsula"
))

Mexmaize_PC2v3.translation<-matrix(rep(Mexmaize_PC2v3_procrustes$translation,1726),nrow=1726,ncol=2,byrow =TRUE) 

Mexmaize_PC2v3_rotation<-Mexmaize_PC2v3_rotation %>% 
  add_column(transformed_PC2 = (Mexmaize_PC2v3_procrustes$scale*t((Mexmaize_PC2v3_procrustes$rotation) %*% t(as.matrix(Mexmaize[,c("PC2","PC3")])))+Mexmaize_PC2v3.translation)[,1],
             transformed_PC3 = (Mexmaize_PC2v3_procrustes$scale*t((Mexmaize_PC2v3_procrustes$rotation) %*% t(as.matrix(Mexmaize[,c("PC2","PC3")])))+Mexmaize_PC2v3.translation)[,2])

Mexmaize_PC2v3_axes<-as.matrix(Mexmaize_PC2v3_procrustes$rotation)%*%as.matrix(axes) %>%
  as_tibble() %>% 
  mutate(maizePC = c("PC2","PC3"),
         slope = (y-yend)/(x-xend))

#to find intercept: y=mx+b
x=Mexmaize_PC2v3_procrustes$translation[1,1]
y=Mexmaize_PC2v3_procrustes$translation[1,2]
y - Mexmaize_PC2v3_axes$slope[1]*x

Mexmaize_PC2v3_axes<-Mexmaize_PC2v3_axes %>% 
  mutate(x_intercept = Mexmaize_PC2v3_procrustes$translation[1,2] - slope*Mexmaize_PC2v3_procrustes$translation[1,1])

ggplot(Mexmaize_PC2v3_rotation)+
  geom_point(aes(x=transformed_PC2, 
                 y=transformed_PC3, 
                 color = Regional_Name), alpha = 0.67)+
  geom_abline(data = Mexmaize_PC2v3_axes, 
              aes(slope = slope, intercept = x_intercept, linetype = maizePC))+
  theme_bw()+
  scale_color_manual(values = c("Baja California"="#CCCCCC",
                                "Mexican Plateau"="#2d1a77",
                                "Sierra Madre Oriental"="#6E54D9",
                                "Sierra Madre Occidental and Pacific Coastal Lowlands"="#b100ea",
                                "Mesa Central"="#fca207",
                                "Gulf Coastal Plain"="#cb4d8e",
                                "Cordillera Neovolcanica"="#a7c957",
                                "Southern Highlands"="#386651",
                                "Yucatan Peninsula"="#268189"),
                     name = "MX State Regions")+xlab("Longitude")+ylab("Latitude")+
  guides(color = "none")+theme(legend.position = "bottom")

ggsave("Plots/Mexmaize_procrustes_raw_latlong_byRegionalName.PC2v3.png",
       dpi = 300, device = "png", width = 3.5, height = 3.75, units = "in")

#now doing PC2 and 3 with human 21km data instead of lat/long
centroid_21km_PC2v3_procrustes<-procrustes(X=select(centroidPCs_Mexhuman21km, PC1.human,PC2.human), 
                                           Y=select(centroidPCs_Mexhuman21km, PC2.maize, PC3.maize))
protest(X=select(centroidPCs_Mexhuman21km, PC1.human,PC2.human), 
        Y=select(centroidPCs_Mexhuman21km, PC2.maize, PC3.maize))
#correlation = 0.4302

centroid_21km_PC2v3_rotation<-tibble(Region = centroidPCs_Mexhuman21km$Region.human,
                                     PC1.human = centroidPCs_Mexhuman21km$PC1.human,
                                     PC2.human = centroidPCs_Mexhuman21km$PC2.human,
                                     PC2.maize = centroidPCs_Mexhuman21km$PC2.maize,
                                     PC3.maize = centroidPCs_Mexhuman21km$PC3.maize) 

centroid_21km_PC2v3.translation<-matrix(rep(centroid_21km_PC2v3_procrustes$translation,345),nrow=345,ncol=2,byrow =TRUE) 

centroid_21km_PC2v3_rotation<-centroid_21km_PC2v3_rotation %>% 
  add_column(transformed_PC2 = (centroid_21km_PC2v3_procrustes$scale*t((centroid_21km_PC2v3_procrustes$rotation) %*% t(as.matrix(centroidPCs_Mexhuman21km[,c("PC2.maize","PC3.maize")])))+centroid_21km_PC2v3.translation)[,1],
             transformed_PC3 = (centroid_21km_PC2v3_procrustes$scale*t((centroid_21km_PC2v3_procrustes$rotation) %*% t(as.matrix(centroidPCs_Mexhuman21km[,c("PC2.maize","PC3.maize")])))+centroid_21km_PC2v3.translation)[,2])

centroid_21km_PC2v3_axes<-as.matrix(centroid_21km_PC2v3_procrustes$rotation)%*%as.matrix(axes) %>%
  as_tibble() %>% 
  mutate(maizePC = c("PC2","PC3"),
         slope = (y-yend)/(x-xend))

#to find intercept: y=mx+b

centroid_21km_PC2v3_axes<-centroid_21km_PC2v3_axes %>% 
  mutate(x_intercept = centroid_21km_PC2v3_procrustes$translation[1,2] - slope*centroid_21km_PC2v3_procrustes$translation[1,1])

centroid_21km_PC2v3_rotation %>%
  select(Region, PC1.human, PC2.human, transformed_PC2, transformed_PC3) %>%
  mutate(ID = 1:nrow(centroid_21km_PC2v3_rotation)) %>% #dummy variable to pair them back up on the other end of the pivots
  pivot_longer(cols = contains("PC"), names_to = "PC_type", values_to = "PC") %>%
  mutate(axes = case_when(PC_type %in% c("PC1.human","transformed_PC2") ~ "x",
                          PC_type %in% c("PC2.human","transformed_PC3") ~ "y"),
         species = case_when(PC_type %in% c("PC1.human","PC2.human") ~ "human",
                             PC_type %in% c("transformed_PC2","transformed_PC3") ~ "maize")) %>%
  pivot_wider(names_from = "axes", values_from = "PC", id_cols = c("Region","species","ID")) %>%
  ggplot(aes(x = x, y=y, shape = species))+
  geom_abline(data = centroid_21km_PC2v3_axes, 
              aes(slope = slope, intercept = x_intercept, linetype = maizePC))+
  geom_point(aes(color = Region), alpha = 0.75)+
  theme_bw()+scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77","#b100ea","#386651","#a7c957"), #All colored
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  coord_fixed()+xlab("Human PC 1")+ylab("Human PC 2")+
  theme(panel.grid = element_blank())+guides(color = "none")+
  scale_shape_manual(values = c("human"=3,"maize"=16))
ggsave("Plots/21km_Mexhuman_Mexmaize_procrustes.PC2v3.raw.png",
       dpi = 300, width = 7, height = 3.4, units = "in")

# Now do it with the admixture removed data
Mexmaize.transf.eigens<-read_tsv("lm_results/Mex_transformed_LDprune_MAF0.01.eigens")
Mexmaize.transf.pca<-read_tsv("lm_results/Mex_transformed_LDprune_MAF0.01.pca")
Mexmaize.transf<-mutate(Mexmaize.transf.pca, sample_id = str_split(ind_id, ":",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

Mexmaize.transf<-inner_join(x=Mexmaize.transf, y=select(Mexmaize_by_states, id, shapeName), join_by(id))

Mexmaize_transf_procrustes<-procrustes(Y = select(Mexmaize.transf, c("PC1","PC2")),
                                       X = select(Mexmaize.transf, c("locations_longitude","locations_latitude"))
)
protest(Y = select(Mexmaize.transf, c("PC1","PC2")),
        X = select(Mexmaize.transf, c("locations_longitude","locations_latitude"))
) #r-squared =0.5387


Mexmaize_transf_rotation<-tibble(MX_State = Mexmaize.transf$shapeName,
                                 latitude = Mexmaize.transf$locations_latitude,
                                 longitude = Mexmaize.transf$locations_longitude,
                                 rotated_PC2.maize = Mexmaize_transf_procrustes$Yrot[,1],
                                 rotated_PC3.maize = Mexmaize_transf_procrustes$Yrot[,2],
                                 PC1.maize = Mexmaize.transf$PC1,
                                 PC2.maize = Mexmaize.transf$PC2) 

Mexmaize_transf_rotation<-Mexmaize_transf_rotation %>% mutate(Regional_Name = case_when(
  MX_State %in% c("Baja California","Baja California Sur") ~ "Baja California",
  MX_State %in% c("Chihuahua","Durango","Zacatecas","Aguascalientes","Coahuila de Zaragoza") ~ "Mexican Plateau",
  MX_State %in% c("Nuevo Leon","San Luis Potosi") ~ "Sierra Madre Oriental",
  MX_State %in% c("Guanajuato","Hidalgo","Queretaro de Arteaga") ~ "Mesa Central",
  MX_State %in% c("Sinaloa","Sonora","Nayarit") ~ "Sierra Madre Occidental and Pacific Coastal Lowlands",
  MX_State %in% c("Jalisco","Colima","Michoacan de Ocampo","Mexico","Distrito Federal","Morelos","Puebla","Tlaxcala") ~ "Cordillera Neovolcanica",
  MX_State %in% c("Tamaulipas","Veracruz de Ignacio de la Llave","Tabasco") ~ "Gulf Coastal Plain",
  MX_State %in% c("Guerrero","Oaxaca","Chiapas") ~ "Southern Highlands",
  MX_State %in% c("Yucatan","Campeche","Quintana Roo") ~ "Yucatan Peninsula"
))

Mexmaize_transf.translation<-matrix(rep(Mexmaize_transf_procrustes$translation,1726),nrow=1726,ncol=2,byrow =TRUE) 

Mexmaize_transf_rotation<-Mexmaize_transf_rotation %>% 
  add_column(transformed_PC1 = (Mexmaize_transf_procrustes$scale*t((Mexmaize_transf_procrustes$rotation) %*% t(as.matrix(Mexmaize.transf[,c("PC1","PC2")])))+Mexmaize_transf.translation)[,1],
             transformed_PC2 = (Mexmaize_transf_procrustes$scale*t((Mexmaize_transf_procrustes$rotation) %*% t(as.matrix(Mexmaize.transf[,c("PC1","PC2")])))+Mexmaize_transf.translation)[,2])

Mexmaize_transf_axes<-as.matrix(Mexmaize_transf_procrustes$rotation)%*%as.matrix(axes) %>%
  as_tibble() %>% 
  mutate(maizePC = c("PC1","PC2"),
         slope = (y-yend)/(x-xend))

#to find intercept: y=mx+b
x=Mexmaize_transf_procrustes$translation[1,1]
y=Mexmaize_transf_procrustes$translation[1,2]
y - Mexmaize_transf_axes$slope[1]*x

Mexmaize_transf_axes<-Mexmaize_transf_axes %>% 
  mutate(x_intercept = Mexmaize_transf_procrustes$translation[1,2] - slope*Mexmaize_transf_procrustes$translation[1,1])

ggplot(Mexmaize_transf_rotation)+
  geom_point(aes(x=transformed_PC1, 
                 y=transformed_PC2, 
                 color = Regional_Name), alpha = 0.67)+
  geom_abline(data = Mexmaize_transf_axes, 
              aes(slope = slope, intercept = x_intercept, linetype = maizePC))+
  theme_bw()+
  scale_color_manual(values = c("Baja California"="#CCCCCC",
                                "Mexican Plateau"="#2d1a77",
                                "Sierra Madre Oriental"="#6E54D9",
                                "Sierra Madre Occidental and Pacific Coastal Lowlands"="#b100ea",
                                "Mesa Central"="#fca207",
                                "Gulf Coastal Plain"="#cb4d8e",
                                "Cordillera Neovolcanica"="#a7c957",
                                "Southern Highlands"="#386651",
                                "Yucatan Peninsula"="#268189"),
                     name = "MX State Regions")+xlab("Longitude")+ylab("Latitude")+
  guides(color = "none")+theme(legend.position = "bottom")+
  coord_equal()

ggsave("Plots/Mexmaize_procrustes_transf_latlong_byRegionalName.PC1v2.png",
       dpi = 300, device = "png", width = 3.5, height = 3.75, units = "in")

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
  geom_tile(aes(x=human_PC, y=maize_PC, fill = abs(correlation)))+
  theme_bw()+theme(axis.text.y = element_blank())+xlab("Human PCs")+ylab("Maize PCs")+
  colorspace::scale_fill_continuous_diverging()

filter(all_by_all_correlations, maize_PC %in% (1:10)) %>%
  ggplot(aes(x=human_PC, y=maize_PC))+
  geom_tile(aes( fill = abs(correlation)))+
  geom_text(aes(label = signif(correlation, digits = 2)),size=12)+
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
#note that the human PCs come from the Mexican only calculated PCs
human_21km_pairs_PCs<-select(Mexhuman, starts_with("PC_"),"sample_name") %>% 
  left_join(x=human_21km_pairs, y= . , by = c("sample_id.human"="sample_name"))

#human_21km_pairs_PCs_transf<-left_join(x=human_21km_pairs_PCs,
#                                       #y=Mexmaize_admixRemoved, 
#                                y=transformed_21km_maize,
#                                      by=c("Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize"="gwas_sample_id"))

raw_21km_maize.pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_non_transformed_LDprune_MAF0.01.pca")
raw_21km_maize.eigen<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_non_transformed_LDprune_MAF.01.eigens")
raw_21km_maize<-left_join(raw_21km_maize.pca, maize_meta, by=c("ind_id"="Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

human_21km_pairs_PCs<-left_join(x=human_21km_pairs_PCs,
                                y=raw_21km_maize,
                                by=c("Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize"="ind_id"))

centroidPCs_Mexhuman21km<-human_21km_pairs_PCs %>% 
  group_by(sample_id.human,Region.human,Ethnicity.human,latitude.human,longitude.human) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)),
            mean_elevation.maize=mean(locations_elevation.maize))
colnames(centroidPCs_Mexhuman21km)[6:15]<-paste0("PC",1:10,".human")
colnames(centroidPCs_Mexhuman21km)[16:(ncol(centroidPCs_Mexhuman21km)-1)]<-paste0("PC",1:422,".maize")

Mexhuman_Mexmaize_procrustes<-procrustes(X=select(ungroup(centroidPCs_Mexhuman21km), PC1.human, PC2.human),
                                         Y=select(ungroup(centroidPCs_Mexhuman21km), PC1.maize, PC2.maize))
protest(X=select(ungroup(centroidPCs_Mexhuman21km), PC1.human, PC2.human),
           Y=select(ungroup(centroidPCs_Mexhuman21km), PC1.maize, PC2.maize))
#0.4304 #transformedd
#0.388 #raw

Mexhuman_Mexmaize.translation<-matrix(rep(Mexhuman_Mexmaize_procrustes$translation,345),nrow=345,ncol=2,byrow =TRUE) 

centroidPCs_Mexhuman21km<-ungroup(centroidPCs_Mexhuman21km)

centroidPCs_Mexhuman21km<-centroidPCs_Mexhuman21km %>% add_column(as.data.frame(Mexhuman_Mexmaize_procrustes$Yrot))
colnames(centroidPCs_Mexhuman21km)[439:440]<-c("rotated_PC1","rotated_PC2")
centroidPCs_Mexhuman21km<-centroidPCs_Mexhuman21km %>%
  mutate(transformed_PC1 = (Mexhuman_Mexmaize_procrustes$scale*t((Mexhuman_Mexmaize_procrustes$rotation) %*% t(as.matrix(centroidPCs_Mexhuman21km[,c("PC1.maize","PC2.maize")])))+Mexhuman_Mexmaize.translation)[,1],
         transformed_PC2 = (Mexhuman_Mexmaize_procrustes$scale*t((Mexhuman_Mexmaize_procrustes$rotation) %*% t(as.matrix(centroidPCs_Mexhuman21km[,c("PC1.maize","PC2.maize")])))+Mexhuman_Mexmaize.translation)[,2])

axes<-tibble(x=c(-100,0),
             xend=c(100,0),
             y=c(0,-100),
             yend=c(0,100))
rotated_axes<-as.matrix(Mexhuman_Mexmaize_procrustes$rotation)%*%as.matrix(axes) %>%
  as_tibble() %>% 
  mutate(maizePC = c("PC1","PC2"),
         slope = (y-yend)/(x-xend))

Mexhuman_mexmaize_procrustes_plot<-
  ggplot(centroidPCs_Mexhuman21km, aes(color = Region.human))+
  geom_point(aes(x=PC1.human, y=PC2.human), shape = 3)+
  geom_point(aes(x=transformed_PC1, y=transformed_PC2), alpha =0.5)+
  geom_abline(data = rotated_axes, 
              aes(slope = slope, intercept = 0, linetype = maizePC))+
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
      "Mexico, Mayan","Mexico, North of \nMesoamerica","Mexico, North",
      "Mexico, Oaxaca","Mexico, West"))+
  coord_fixed()+
    scale_x_continuous(breaks = seq(-200, 800, by=100))+
    scale_y_continuous(breaks = seq(-200,300, by=100))#+
  #theme(legend.position = "bottom")

  ggsave(plot = Mexhuman_mexmaize_procrustes_plot, 
       filename = "/group/jrigrp11/snodgras_maizePopulations/Plots/21km_Mexhuman_Mexmaize_procrustes.raw.png", dpi = 300, 
       device = "png", width = 7, height = 5)  

Mex21_all_by_all_correlations<-cor(centroidPCs_Mexhuman21km[,6:15],centroidPCs_Mexhuman21km[,16:(ncol(centroidPCs_Mexhuman21km)-1)])
Mex21_all_by_all_correlations<-Mex21_all_by_all_correlations %>% as_tibble() %>% rownames_to_column(var = "human_PC")
Mex21_all_by_all_correlations<-pivot_longer(data = Mex21_all_by_all_correlations, cols = ends_with(".maize"), names_to = "maize_PC", values_to = "correlation")
Mex21_all_by_all_correlations<-mutate(Mex21_all_by_all_correlations, maize_PC = str_remove(maize_PC,".maize") %>% str_remove("PC"))
Mex21_all_by_all_correlations$human_PC <- factor(Mex21_all_by_all_correlations$human_PC, levels = 1:10)
Mex21_all_by_all_correlations$maize_PC <- factor(Mex21_all_by_all_correlations$maize_PC, levels = 1:422)
Mex21_all_by_all_correlations<-mutate(Mex21_all_by_all_correlations, r_squared = correlation^2)

ggplot(Mex21_all_by_all_correlations)+
  geom_tile(aes(x=human_PC, y=maize_PC, fill = abs(correlation)))+
  theme_bw()+xlab("Mexican Human PCs")+ylab("Mexican Maize PCs")+
  colorspace::scale_fill_continuous_diverging()+labs(fill = "Correlation")+
  scale_y_discrete(labels = c("1",rep("",420),"422"))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/2026-01-23-21km_MexicanSamples_PCcorrelationheatmap.all_v_all.raw.png",
       device = "png", dpi = 300, width = 4,height = 3.5, units = "in")

#Dan's test
Dan_cor_test<-human_21km_pairs_PCs %>% 
  group_by(Region.human,latitude.human,longitude.human) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)),
            mean_elevation.maize=mean(locations_elevation.maize))
colnames(Dan_cor_test)[4:13]<-paste0("PC",1:10,".human")
colnames(Dan_cor_test)[14:(ncol(Dan_cor_test)-1)]<-paste0("PC",1:422,".maize")

Dan_Mex21_all_by_all_correlations<-cor(Dan_cor_test[,4:13],Dan_cor_test[,14:(ncol(Dan_cor_test)-1)])
Dan_Mex21_all_by_all_correlations<-Dan_Mex21_all_by_all_correlations %>% as_tibble() %>% rownames_to_column(var = "human_PC")
Dan_Mex21_all_by_all_correlations<-pivot_longer(data = Dan_Mex21_all_by_all_correlations, cols = ends_with(".maize"), names_to = "maize_PC", values_to = "correlation")
Dan_Mex21_all_by_all_correlations<-mutate(Dan_Mex21_all_by_all_correlations, maize_PC = str_remove(maize_PC,".maize") %>% str_remove("PC"))
Dan_Mex21_all_by_all_correlations$human_PC <- factor(Dan_Mex21_all_by_all_correlations$human_PC, levels = 1:10)
Dan_Mex21_all_by_all_correlations$maize_PC <- factor(Dan_Mex21_all_by_all_correlations$maize_PC, levels = 1:422)
Dan_Mex21_all_by_all_correlations<-mutate(Dan_Mex21_all_by_all_correlations, r_squared = correlation^2)

ggplot(Dan_Mex21_all_by_all_correlations)+
  geom_tile(aes(x=human_PC, y=maize_PC, fill = r_squared))+
  theme_bw()+xlab("Mexican Human PCs")+ylab("Mexican Maize PCs")+
  colorspace::scale_fill_continuous_diverging()+labs(fill = "Correlation")+
  scale_y_discrete(labels = c("1",rep("",420),"422"))
ggsave("/group/jrigrp11/snodgras/snodgras_maizePopulations/Plots/2026-01-24-21km_MexicanSamples_PCcorrelation-heatmap.all_v_all.raw.PCsCentroid.png",
       width = 1200, height = 1050, units = "px", dpi = 300)

summary(Dan_Mex21_all_by_all_correlations$r_squared)
arrange(Dan_Mex21_all_by_all_correlations, -r_squared)

make_scatterplot<-function(HPC, MPC){
  plt<-ggplot(Dan_cor_test, aes_string(HPC, MPC))+
    geom_smooth(method = "lm", color = "black")+
    geom_point(aes(color = Region.human), alpha = 0.75)+
    theme_bw()+scale_color_manual(
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
  return(plt)
}
make_scatterplot("PC4.human", "PC3.maize") #r_squared =0.217
make_scatterplot("PC1.human", "PC3.maize") #r_squared =0.213
make_scatterplot("PC4.human", "PC6.maize") #r_squared =0.206
make_scatterplot("PC1.human", "PC128.maize") #r_squared =0.205
make_scatterplot("PC1.human", "PC2.maize") #r_squared =0.204
# End of Dan's test

ggplot(Mex21_all_by_all_correlations)+
  geom_tile(aes(x=human_PC, y=maize_PC, fill = r_squared))+
  theme_bw()+theme(axis.text.y = element_blank())+xlab("Mexican Human PCs")+ylab("Mexican Maize Raw PCs")+
  colorspace::scale_fill_continuous_diverging()

filter(Mex21_all_by_all_correlations, maize_PC %in% (1:10)) %>%
  ggplot(aes(x=human_PC, y=maize_PC))+
  geom_tile(aes( fill = correlation))+
  geom_text(aes(label = signif(correlation, digits = 3)),size=5)+
  theme_bw()+xlab("Mexican Human PCs")+ylab("Mexican Maize Raw PCs")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size = 20, face="bold"))+
  colorspace::scale_fill_continuous_diverging()
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/2025-10-22-21km_MexicanSamples_PCcorrelationheatmap.PCs1_10.raw.png",
       device = "png", dpi = 300, width = 11,height = 8)

filter(Mex21_all_by_all_correlations, maize_PC %in% (1:10)) %>%
  ggplot(aes(x=human_PC, y=maize_PC))+
  geom_tile(aes( fill = r_squared))+
  geom_text(aes(label = signif(r_squared, digits = 2)),size=5)+
  theme_bw()+xlab("Mexican Human PCs")+ylab("Mexican Maize Raw PCs")+
  theme(axis.text = element_text(size=20), axis.title = element_text(size = 20, face="bold"))+
  colorspace::scale_fill_continuous_diverging()

filter(Mex21_all_by_all_correlations, correlation >= 0.5 | correlation <= -0.5) # 427 PC pairs met this criteria

ggplot(centroidPCs_Mexhuman21km)+
  geom_point(aes(x=PC2.human, y=PC7.maize, color = Region.human))+
  theme_bw()

ggplot(centroidPCs_Mexhuman21km, aes(x=PC4.human, y=PC6.maize))+
  geom_smooth(method = "lm", color = "black")+
  geom_point(aes(color = Region.human), alpha = 0.75)+
  theme_bw()+scale_color_manual(
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

#originally had only looked at the first 20 PCs

plotting_pairs<-filter(Mex21_all_by_all_correlations, correlation >= 0.6)
pdf("/group/jrigrp11/snodgras_maizePopulations/Plots/2025-12-08-21km_MexicanSamples_PCcorrelationPoints.raw.pdf",
    height = 4, width = 4)

ggplot(centroidPCs_Mexhuman21km, aes(x=longitude.human,y=latitude.human, color=Region.human))+
  geom_point()+theme_bw()+ggtitle("Locations of each human 'Region' in Mexico")+
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
      "Mexico, Mayan","Mexico, North \nof Mesoamerica","Mexico, North",
      "Mexico, Oaxaca","Mexico, West"))+
  ggtitle("Map of Human Samples\nColor Guide")+
  theme(plot.title = element_text(size=12), 
        legend.text = element_text(size=6), 
        legend.title = element_text(size=8))+
  xlab("Long.")+ylab("Lat.")
  

for(i in 1:nrow(plotting_pairs)){
  plt<-ggplot(centroidPCs_Mexhuman21km, aes_string(x=paste0("PC",plotting_pairs$human_PC[i],".human"), y=paste0("PC",plotting_pairs$maize_PC[i],".maize")))+
    geom_smooth(method = "lm", color = "black")+
    geom_point(aes(color=Region.human), alpha = 0.75)+
    theme_bw()+scale_color_manual(
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
        "Mexico, Oaxaca","Mexico, West"))+
    guides(color="none")+theme(plot.title = element_text(size=12))+
    ggtitle(paste0("Cor. = ",signif(plotting_pairs$correlation[i], 3),"; R-sq. = ", signif((plotting_pairs$correlation[i])^2, 3)))
  print(plt)
}

dev.off()

ggplot(centroidPCs_Mexhuman21km, aes(color=Region.human,x = PC1.human, y=PC2.maize))+
  geom_smooth(method = "lm", col="black")+
  geom_point(alpha =0.75)+
  theme_bw()+ 
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  guides(color = "none")
#To try and figure out what's going on with PC 3
Mexhuman %>% ggplot(aes(x=PC_2, y=PC_3, color=Ethnicity,shape = Region))+geom_point()+scale_shape_manual(values=c(0,1,2,3,4,5,10))
#are there any correlations that continue if the North Mexico Samples are removed?
Mex21_woNMex<-filter(centroidPCs_Mexhuman21km, Region.human != "Mexico-North_of_Mexico")
Mex21_woNMex_correlations<-cor(Mex21_woNMex[,6:15],Mex21_woNMex[,16:437])
Mex21_woNMex_correlations<-Mex21_woNMex_correlations %>% as_tibble() %>% rownames_to_column(var = "human_PC")
Mex21_woNMex_correlations<-pivot_longer(data = Mex21_woNMex_correlations, cols = ends_with(".maize"), names_to = "maize_PC", values_to = "correlation")
Mex21_woNMex_correlations<-mutate(Mex21_woNMex_correlations, maize_PC = str_remove(maize_PC,".maize") %>% str_remove("PC"))
Mex21_woNMex_correlations$human_PC <- factor(Mex21_woNMex_correlations$human_PC, levels = 1:10)
Mex21_woNMex_correlations$maize_PC <- factor(Mex21_woNMex_correlations$maize_PC, levels = 1:422)
Mex21_woNMex_correlations<-mutate(Mex21_woNMex_correlations, r_squared = correlation^2)

ggplot(Mex21_woNMex_correlations)+
  geom_tile(aes(x=human_PC, y=maize_PC, fill = correlation))+
  theme_bw()+theme(axis.text.y = element_blank())+xlab("Mexican Human PCs")+ylab("Mexican Maize Raw PCs")+
  colorspace::scale_fill_continuous_diverging()

ggplot(Mex21_woNMex_correlations)+
  geom_tile(aes(x=human_PC, y=maize_PC, fill = r_squared))+
  theme_bw()+theme(axis.text.y = element_blank())+xlab("Mexican Human PCs")+ylab("Mexican Maize Raw PCs")+
  colorspace::scale_fill_continuous_diverging()

plotting_pairs_woNMex<-filter(Mex21_woNMex_correlations, r_squared >= 0.4)
pdf("/group/jrigrp11/snodgras_maizePopulations/Plots/2025-10-22-21km_MexicanSamples_woNMex_PCcorrelationPoints.raw.pdf")

for(i in 1:nrow(plotting_pairs_woNMex)){
  plt<-ggplot(Mex21_woNMex, aes_string(x=paste0("PC",plotting_pairs_woNMex$human_PC[i],".human"), y=paste0("PC",plotting_pairs_woNMex$maize_PC[i],".maize")))+
    geom_smooth(method = "lm", color = "black")+
    geom_point(aes(color=Region.human), alpha = 0.75)+
    theme_bw()+scale_color_manual(
      values = c("Mexico-Center_of_Mexico"="#fca207",
                 "Mexico-Gulf_of_Mexico"="#cb4d8e",
                 "Mexico-Mayan_region"="#268189",
                 "Mexico-North_of_Mesoamerica"="#2d1aff",
                 "Mexico-Oaxaca"="#386955",
                 "Mexico-West_of_Mexico"="#a7c957"),
      name = "Region", labels = c(
        "Mexico, Center","Mexico, Gulf", 
        "Mexico, Mayan","Mexico, North of Mesoamerica",
        "Mexico, Oaxaca","Mexico, West"))+
    ggtitle(paste0("Correlation = ",signif(plotting_pairs_woNMex$correlation[i], 3),"; R-squared = ", signif((plotting_pairs_woNMex$correlation[i])^2, 3)))
  print(plt)
}

dev.off()

#which sets of PCs have high correlations no matter if NMex group is removed?
pair_set1<-unique(select(plotting_pairs, human_PC, maize_PC)) %>% mutate(pair = str_c(human_PC,",",maize_PC))
pair_set2<-unique(select(plotting_pairs_woNMex, human_PC, maize_PC))%>% mutate(pair = str_c(human_PC,",",maize_PC))

filter(pair_set1, pair %in% pair_set2$pair) #kept the same between with North of Mexico and without that group
filter(pair_set1, !pair %in% pair_set2$pair) #only found with North of Mexico
filter(pair_set2, !pair %in% pair_set1$pair) #only found without North of Mexico

#what is human PC3 and is any human PC correlated with elevation
cor(centroidPCs_Mexhuman21km[,6:15],centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])
#PC 3 -0.15, PC6 0.13, PC10 -0.213

cor.test(pull(centroidPCs_Mexhuman21km[,8]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])) #PC3, significant p=0.005
cor.test(pull(centroidPCs_Mexhuman21km[,11]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])) #PC6 significant p=0.014
cor.test(pull(centroidPCs_Mexhuman21km[,15]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])) #PC10 significant p=6.5e-5
cor.test(pull(centroidPCs_Mexhuman21km[,6]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)])) #PC1 not significant p=0.1257

for(i in 6:15){
  cor.test.result<-cor.test(pull(centroidPCs_Mexhuman21km[,i]),pull(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)]))
  if(cor.test.result$p.value <= 0.05/10){
    print(paste(colnames(centroidPCs_Mexhuman21km)[i], 
                "is significantly correlated with", 
                colnames(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)]),";",
                "estimated correlation =", signif(cor.test.result$estimate, 3),
                "p value =", signif(cor.test.result$p.value, 3)))
  }else{
    print(paste(colnames(centroidPCs_Mexhuman21km)[i], 
                "is NOT significantly correlated with", 
                colnames(centroidPCs_Mexhuman21km[,ncol(centroidPCs_Mexhuman21km)]),";",
                "estimated correlation =", signif(cor.test.result$estimate, 3),
                "p value =", signif(cor.test.result$p.value, 3)))
  }
}

pdf("/group/jrigrp11/snodgras_maizePopulations/Plots/2025-06-06-21km_MexicanSamples_HumanElevationPCcorrPoints.pdf",
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

pdf("/group/jrigrp11/snodgras_maizePopulations/Plots/2025-06-06-21km_MexicanSamples_MaizePCsGeoCorrWithHumanPCs.pdf",
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

#### after MaizePVEbyHumanPCs.R where each Maize snp is lm'ed with the first 10PCs of humans ####
maize21km_human_lm_Rsq<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01.genotypes.all_models.lm_results.tsv")
transformed_maize21km_human_lm_Rsq<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_transformed_LDprune_MAF0.01.genotypes.all_models.lm_results.tsv")

#ggplot(maize21km_human_lm_Rsq)+
#  geom_histogram(aes(x=r.squared, fill = model_formula), binwidth = 0.01)+
#  geom_vline(xintercept = 0)+theme_bw()

#clean up model names
maize21km_human_lm_Rsq <- maize21km_human_lm_Rsq %>% 
  mutate(model_formula = str_remove_all(model_formula, ".maize"),
         model_formula = str_remove_all(model_formula, "locations_"),
         model_formula = str_replace_all(model_formula, "latitude","lat"),
         model_formula = str_replace_all(model_formula, "longitude", "long"),
         model_formula = str_replace_all(model_formula, "elevation", "elev")) 

transformed_maize21km_human_lm_Rsq<-transformed_maize21km_human_lm_Rsq %>% 
  mutate(model_formula = str_remove_all(model_formula, ".maize"),
         model_formula = str_remove_all(model_formula, "locations_"),
         model_formula = str_replace_all(model_formula, "latitude","lat"),
         model_formula = str_replace_all(model_formula, "longitude", "long")) 

#What's the max/min r-squared value for each model?
group_by(maize21km_human_lm_Rsq, model_formula) %>% 
  summarize(max.value = max(r.squared), 
            min.value = min(r.squared),
            mean.value = mean(r.squared, na.rm = T)) %>%
  ggplot(aes(x=model_formula))+
  geom_point(aes(y=max.value))+
  geom_point(aes(y=min.value))+
  geom_point(aes(y=mean.value), shape = 23, color = "red", fill = "red")+
  geom_text(aes(y=max.value+0.05,label = signif(max.value, 3)))+
  geom_text(aes(y=mean.value+0.05,label = signif(mean.value, 2)), color = "darkred")+
  geom_text(aes(y=min.value-0.05,label = signif(min.value, 1)))+
  coord_flip()+xlab("")+ylab("R-squared")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Min,Mean, Max R-squared, Raw Genotypes")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_raw_points_MinMeanMax.png", 
       width = 11, height = 8.5, dpi =300, device = "png")

#### MAKING SUPPLEMENTAL FIGURE HERE ####
group_by(maize21km_human_lm_Rsq, model_formula) %>% 
  summarize(max.value = max(adj.r.squared), 
            min.value = min(adj.r.squared),
            mean.value = mean(adj.r.squared, na.rm = T)) %>%
  ungroup() %>%
  filter(!model_formula %in% c("lat+long+PC_1","lat+long+PC_1+PC_2","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10",
                               "lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8",
                               "lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6",
                               "lat+long+PC_1+PC_2+PC_3+PC_4+PC_5","lat+long+PC_1+PC_2+PC_3+PC_4",
                               "lat+long+PC_1+PC_2+PC_3"))%>%
  mutate(PCs = case_when(model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10","lat+long+elev+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10") ~ 10,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9","lat+long+elev+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9") ~ 9,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8","lat+long+elev+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8") ~ 8,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7","lat+long+elev+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7") ~ 7,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6","lat+long+elev+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6") ~ 6,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5","lat+long+elev+PC_1+PC_2+PC_3+PC_4+PC_5") ~ 5,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4","lat+long+elev+PC_1+PC_2+PC_3+PC_4") ~ 4,
                         model_formula %in% c("PC_1+PC_2+PC_3","lat+long+elev+PC_1+PC_2+PC_3") ~ 3,
                         model_formula %in% c("PC_1+PC_2","lat+long+elev+PC_1+PC_2") ~ 2,
                         model_formula %in% c("PC_1","lat+long+elev+PC_1") ~ 1,
                         model_formula %in% c("lat+long+elev","lat+long","elev") ~ 0,
  ),
  model_type = case_when(str_detect(model_formula, "lat|elev")~"Geography (LLE)",
                         str_detect(model_formula, "lat|elev", negate = TRUE) ~ "Human Only")) %>%
  filter(!model_formula %in% c("lat+long", "elev")) %>%
  ggplot(aes(x=PCs))+
  geom_point(aes(y=max.value))+
  geom_point(aes(y=min.value))+
  geom_point(aes(y=mean.value), shape = 23, color = "red", fill = "red")+
  geom_text(aes(y=max.value+0.1,label = signif(max.value, 3)))+
  geom_text(aes(y=mean.value+0.1,label = signif(mean.value, 2)), color = "darkred")+
  geom_text(aes(y=min.value-0.1,label = signif(min.value, 1)))+
  coord_flip()+xlab("")+ylab("Adj. R-squared")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Min,Mean, Max Adj.R-squared, Raw Genotypes")+
  facet_wrap(vars(model_type))+
  scale_x_continuous(breaks = c(0:10), 
                     labels = c("LLE + 0 PC","1 PC","2 PCs","3 PCs","4 PCs","5 PCs","6 PCs","7 PCs","8 PCs","9 PCs","10 PCs"))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_raw_adjR2_points_MinMeanMax.png",
       dpi = 300, device = "png", width = 11, height = 3.75)  

group_by(transformed_maize21km_human_lm_Rsq, model_formula) %>% 
  summarize(max.value = max(r.squared), 
            min.value = min(r.squared),
            mean.value = mean(r.squared, na.rm = T)) %>%
  ggplot(aes(x=model_formula))+
  geom_point(aes(y=max.value))+
  geom_point(aes(y=min.value))+
  geom_point(aes(y=mean.value), shape = 23, color = "red", fill = "red")+
  geom_text(aes(y=max.value+0.05,label = signif(max.value, 3)))+
  geom_text(aes(y=mean.value+0.05,label = signif(mean.value, 2)), color = "darkred")+
  geom_text(aes(y=min.value-0.05,label = signif(min.value, 1)))+
  coord_flip()+xlab("")+ylab("R-squared")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Min,Mean, Max R-squared, Transformed Genotypes")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_transformed_points_MinMeanMax.png", 
       width = 11, height = 8.5, dpi =300, device = "png")

group_by(transformed_maize21km_human_lm_Rsq, model_formula) %>% 
  summarize(max.value = max(adj.r.squared), 
            min.value = min(adj.r.squared),
            mean.value = mean(adj.r.squared, na.rm = T)) %>%
  ungroup() %>%
  mutate(PCs = case_when(model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10") ~ 10,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9") ~ 9,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8") ~ 8,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7") ~ 7,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5+PC_6","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6") ~ 6,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4+PC_5","lat+long+PC_1+PC_2+PC_3+PC_4+PC_5") ~ 5,
                         model_formula %in% c("PC_1+PC_2+PC_3+PC_4","lat+long+PC_1+PC_2+PC_3+PC_4") ~ 4,
                         model_formula %in% c("PC_1+PC_2+PC_3","lat+long+PC_1+PC_2+PC_3") ~ 3,
                         model_formula %in% c("PC_1+PC_2","lat+long+PC_1+PC_2") ~ 2,
                         model_formula %in% c("PC_1","lat+long+PC_1") ~ 1,
                         model_formula %in% c("lat+long") ~ 0,
  ),
  model_type = case_when(str_detect(model_formula, "lat")~"Geography (Latitude, Longitude)",
                         str_detect(model_formula, "lat", negate = TRUE) ~ "Human Only")) %>%
  ggplot(aes(x=PCs))+
  geom_point(aes(y=max.value))+
  geom_point(aes(y=min.value))+
  geom_point(aes(y=mean.value), shape = 23, color = "red", fill = "red")+
  geom_text(aes(y=max.value+0.1,label = signif(max.value, 3)))+
  geom_text(aes(y=mean.value+0.1,label = signif(mean.value, 2)), color = "darkred")+
  geom_text(aes(y=min.value-0.1,label = signif(min.value, 1)))+
  coord_flip()+xlab("")+ylab("Adj. R-squared")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  ggtitle("Min,Mean, Max Adj.R-squared, Admix Removed Genotypes")+
  facet_wrap(vars(model_type))+
  scale_x_continuous(breaks = c(0:10), 
                     labels = c("Lat,Long + 0 PC","1 PC","2 PCs","3 PCs","4 PCs","5 PCs","6 PCs","7 PCs","8 PCs","9 PCs","10 PCs"))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_transf_adjR2_points_MinMeanMax.png",
       dpi = 300, device = "png", width = 11, height = 3.75)  
  
# on average across all SNPs, how much does human PCs explain maize genotypes
maize21km_human_lm_Rsq %>% group_by(model_formula) %>%
  summarize(mean.r.squared = mean(r.squared, na.rm = TRUE),
            mean.adj.r.squared = mean(adj.r.squared, na.rm = TRUE),
            sd.r.squared = sd(r.squared, na.rm = TRUE),
            sd.adj.r.squared = sd(adj.r.squared, na.rm = TRUE))
# so on average human pcs explain ~1-3% of the maize genotypes?

transformed_maize21km_human_lm_Rsq %>% group_by(model_formula) %>%
  summarize(mean.r.squared = mean(r.squared, na.rm = TRUE),
            mean.adj.r.squared = mean(adj.r.squared, na.rm = TRUE),
            sd.r.squared = sd(r.squared, na.rm = TRUE),
            sd.adj.r.squared = sd(adj.r.squared, na.rm = TRUE))
#c(0.0263 - (3*0.0184), 0.0263 + (3*0.0184))

all_lm_Rsq<-left_join(maize21km_human_lm_Rsq, 
           transformed_maize21km_human_lm_Rsq, 
           by=c("snp_ID", "model_formula"), 
           suffix = c(".raw",".transf"))

#### Line graph across all linear models ####
summary_all_lm_Rsq<-all_lm_Rsq %>% group_by(model_formula) %>%
  summarize(mean.raw = mean(r.squared.raw, na.rm = T),
            mean.transf = mean(r.squared.transf, na.rm = T),
            sd.raw = sd(r.squared.raw, na.rm = T),
            sd.transf = sd(r.squared.transf, na.rm = T),
            min.raw = min(r.squared.raw),
            min.transf = min(r.squared.transf),
            max.raw = max(r.squared.raw),
            max.transf = max(r.squared.transf)
            )

summary_all_lm_Rsq<-summary_all_lm_Rsq %>% 
  mutate(model_type = case_when(str_detect(model_formula,fixed("+elev")) ~ "Lat, Long, Elev.", 
                                str_equal(model_formula, "elev") ~ "Elevation Only",
                                str_detect(model_formula, "lat") & str_detect(model_formula, "elev", negate = T) ~ "Lat, Long",
                                str_detect(model_formula, "lat", negate = T) & str_detect(model_formula, "elev", negate = T) ~ "Human PC only"),
         model_type = factor(model_type, levels = c("Human PC only", "Lat, Long","Lat, Long, Elev.", "Elevation Only")),
         PC_number = case_when(str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2", negate = T)
                               & str_detect(model_formula, "PC_3", negate = T)
                               & str_detect(model_formula, "PC_4", negate = T)
                               & str_detect(model_formula, "PC_5", negate = T)
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 1,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3", negate = T)
                               & str_detect(model_formula, "PC_4", negate = T)
                               & str_detect(model_formula, "PC_5", negate = T)
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 2,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4", negate = T)
                               & str_detect(model_formula, "PC_5", negate = T)
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 3,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5", negate = T)
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 4,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 5,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 6,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7")
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 7,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7")
                               & str_detect(model_formula, "PC_8")
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 8,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7")
                               & str_detect(model_formula, "PC_8")
                               & str_detect(model_formula, "PC_9")
                               & str_detect(model_formula, "PC_10", negate = T) ~ 9,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7")
                               & str_detect(model_formula, "PC_8")
                               & str_detect(model_formula, "PC_9")
                               & str_detect(model_formula, "PC_10") ~ 10,
                               .default = 0)
         )

ggplot(summary_all_lm_Rsq, aes(x=PC_number))+
  geom_line(aes(y=mean.raw, color = model_type))+
  geom_point(aes(y=mean.raw, color = model_type, shape = model_type))+
  geom_line(aes(y=mean.raw+sd.raw, color = model_type),linetype = 2)+
  geom_line(aes(y=mean.raw-sd.raw, color = model_type), linetype = 2)+
  #geom_segment(aes(x=PC_number, xend = PC_number, y=mean.raw-sd.raw, yend=mean.raw+sd.raw, color = model_type))+
  geom_point(aes(y=mean.raw+sd.raw, color = model_type, shape = model_type))+
  geom_point(aes(y=mean.raw-sd.raw, color = model_type, shape = model_type))+
  ggtitle("Raw maize genotype, r-squared")+ylab("R-squared")+xlab("Human PCs Included in Model")+
  scale_x_continuous(limits = c(0,10),breaks = c(0:10))+
  theme_bw()
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_raw_linegraph_mean_1sd.png", 
       width = 11, height = 8.5, dpi =300, device = "png")

filter(summary_all_lm_Rsq, model_type %in% c("Human PC only","Lat, Long")) %>% 
  ggplot(aes(x=PC_number))+
  geom_line(aes(y=mean.transf, color = model_type))+
  geom_point(aes(y=mean.transf, color = model_type, shape = model_type))+
  geom_line(aes(y=mean.transf+sd.transf, color = model_type),linetype = 2)+
  geom_line(aes(y=mean.transf-sd.transf, color = model_type), linetype = 2)+
  #geom_segment(aes(x=PC_number, xend = PC_number, y=mean.transf-sd.transf, yend=mean.transf+sd.transf, color = model_type))+
  geom_point(aes(y=mean.transf+sd.transf, color = model_type, shape = model_type))+
  geom_point(aes(y=mean.transf-sd.transf, color = model_type, shape = model_type))+
  ggtitle("Transformed maize genotype, r-squared")+ylab("R-squared")+xlab("Human PCs Included in Model")+
  scale_x_continuous(limits = c(0,10),breaks = c(0:10))+
  theme_bw()
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_transf_linegraph_mean_1sd.png", 
       width = 11, height = 8.5, dpi =300, device = "png")

#### LM for adjusted r-squared ####
summary_all_lm_Rsq_adj<-all_lm_Rsq %>% group_by(model_formula) %>%
  summarize(mean.raw = mean(adj.r.squared.raw, na.rm = T),
            mean.transf = mean(adj.r.squared.transf, na.rm = T),
            sd.raw = sd(adj.r.squared.raw, na.rm = T),
            sd.transf = sd(adj.r.squared.transf, na.rm = T),
            min.raw = min(adj.r.squared.raw),
            min.transf = min(adj.r.squared.transf),
            max.raw = max(adj.r.squared.raw),
            max.transf = max(adj.r.squared.transf),
            stderr.raw = sd.raw / sqrt(nrow(all_lm_Rsq)),
            stderr.transf = sd.transf / sqrt(nrow(all_lm_Rsq))
  )

summary_all_lm_Rsq_adj<-summary_all_lm_Rsq_adj %>% 
  mutate(model_type = case_when(str_detect(model_formula,fixed("+elev")) ~ "Lat, Long, Elev.", 
                                str_equal(model_formula, "elev") ~ "Elevation Only",
                                str_detect(model_formula, "lat") & str_detect(model_formula, "elev", negate = T) ~ "Lat, Long",
                                str_detect(model_formula, "lat", negate = T) & str_detect(model_formula, "elev", negate = T) ~ "Human PC only"),
         model_type = factor(model_type, levels = c("Human PC only", "Lat, Long","Lat, Long, Elev.", "Elevation Only")),
         PC_number = case_when(str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2", negate = T)
                               & str_detect(model_formula, "PC_3", negate = T)
                               & str_detect(model_formula, "PC_4", negate = T)
                               & str_detect(model_formula, "PC_5", negate = T)
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 1,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3", negate = T)
                               & str_detect(model_formula, "PC_4", negate = T)
                               & str_detect(model_formula, "PC_5", negate = T)
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 2,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4", negate = T)
                               & str_detect(model_formula, "PC_5", negate = T)
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 3,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5", negate = T)
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 4,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6", negate = T)
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 5,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7", negate = T)
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 6,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7")
                               & str_detect(model_formula, "PC_8", negate = T)
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 7,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7")
                               & str_detect(model_formula, "PC_8")
                               & str_detect(model_formula, "PC_9", negate = T)
                               & str_detect(model_formula, "PC_10", negate = T) ~ 8,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7")
                               & str_detect(model_formula, "PC_8")
                               & str_detect(model_formula, "PC_9")
                               & str_detect(model_formula, "PC_10", negate = T) ~ 9,
                               str_detect(model_formula, "PC_1") 
                               & str_detect(model_formula, "PC_2")
                               & str_detect(model_formula, "PC_3")
                               & str_detect(model_formula, "PC_4")
                               & str_detect(model_formula, "PC_5")
                               & str_detect(model_formula, "PC_6")
                               & str_detect(model_formula, "PC_7")
                               & str_detect(model_formula, "PC_8")
                               & str_detect(model_formula, "PC_9")
                               & str_detect(model_formula, "PC_10") ~ 10,
                               .default = 0)
  )

ggplot(summary_all_lm_Rsq_adj, aes(x=PC_number))+
  geom_line(aes(y=mean.raw, color = model_type))+
  geom_point(aes(y=mean.raw, color = model_type, shape = model_type))+
  geom_line(aes(y=mean.raw+stderr.raw, color = model_type),linetype = 2)+
  geom_line(aes(y=mean.raw-stderr.raw, color = model_type), linetype = 2)+
  #geom_segment(aes(x=PC_number, xend = PC_number, y=mean.raw-sd.raw, yend=mean.raw+sd.raw, color = model_type))+
  geom_point(aes(y=mean.raw+stderr.raw, color = model_type, shape = model_type))+
  geom_point(aes(y=mean.raw-stderr.raw, color = model_type, shape = model_type))+
  ggtitle("Raw maize genotype, adjusted r-squared")+ylab("Adj. R-squared")+xlab("Human PCs Included in Model")+
  scale_x_continuous(limits = c(0,10),breaks = c(0:10))+
  theme_bw()
#ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_raw_adjRsq_linegraph_mean_1sd.png", 
#       width = 11, height = 8.5, dpi =300, device = "png")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_raw_adjRsq_linegraph_mean_1stderr.png", 
              width = 11, height = 8.5, dpi =300, device = "png")

filter(summary_all_lm_Rsq_adj, model_type %in% c("Human PC only","Lat, Long")) %>% 
  ggplot(aes(x=PC_number))+
  geom_line(aes(y=mean.transf, color = model_type))+
  geom_point(aes(y=mean.transf, color = model_type, shape = model_type))+
  geom_line(aes(y=mean.transf+sd.transf, color = model_type),linetype = 2)+
  geom_line(aes(y=mean.transf-sd.transf, color = model_type), linetype = 2)+
  #geom_segment(aes(x=PC_number, xend = PC_number, y=mean.transf-sd.transf, yend=mean.transf+sd.transf, color = model_type))+
  geom_point(aes(y=mean.transf+sd.transf, color = model_type, shape = model_type))+
  geom_point(aes(y=mean.transf-sd.transf, color = model_type, shape = model_type))+
  ggtitle("Transformed maize genotype, adjusted r-squared")+ylab("Adj. R-squared")+xlab("Human PCs Included in Model")+
  scale_x_continuous(limits = c(0,10),breaks = c(0:10))+
  theme_bw()
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/lm_results_transf_adjRsq_linegraph_mean_1sd.png", 
       width = 11, height = 8.5, dpi =300, device = "png")

#### LM with admixture ####
LLEA_lm_results<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmix_lm_results.tsv")
LLEA_lm_results<-add_column(LLEA_lm_results, model_formula = "LLEA")

for(i in c("PC1","PC1_2","PC1_3","PC1_4","PC1_5","PC1_6","PC1_7","PC1_8","PC1_9","PC1_10")){
  LLEA_lm_results<-add_row(LLEA_lm_results, 
                           read_tsv(paste0("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmix",i,"_lm_results.tsv")),
                           model_formula = paste0("LLEA_",i))
}
summary_LLEA_lm<-LLEA_lm_results %>% group_by(model_formula) %>%
  summarize(mean = mean(adj.r.squared, na.rm = T),
            sd = sd(adj.r.squared, na.rm = T),
            min = min(adj.r.squared),
            max = max(adj.r.squared),
            stderr = sd / sqrt(nrow(filter(LLEA_lm_results, model_formula == "LLEA")))
  )

summary_LLEA_lm<-mutate(summary_LLEA_lm, 
                        model_formula = factor(model_formula, levels = c("LLEA","LLEA_PC1",paste0("LLEA_PC1_",2:10))))

ggplot(summary_LLEA_lm, aes(x=model_formula,y=mean))+
  geom_line(group=1)+
  geom_line(aes(x=c(1:11),
                y=c(Mexmaize_eigenval$value[1],cumsum(Mexmaize_eigenval$value[1:10]))),
            group=1,
            linetype=2)+
  geom_point()+
  geom_point(aes(x=c("LLEA","LLEA_PC1",paste0("LLEA_PC1_",2:10)),
                 y=c(Mexmaize_eigenval$value[1],cumsum(Mexmaize_eigenval$value[1:10]))),
            shape=2)+
  geom_segment(aes(x=model_formula, xend=model_formula,
                   y=mean-3*stderr, yend=mean+3*stderr))+
  theme_bw()+
  ylab("Adj. R-squared")+xlab("Model Parameters")+
  theme(axis.text.x = element_text(angle = 90))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_LLEA_humanPCs_adjRsquared_linegraph.png", 
       dpi=300, device="png", height = 4, width = 5)

summary_LLEA_lm %>% filter(model_formula %in% c("LLEA","LLEA_PC1_10")) %>%
  ggplot(aes(x=model_formula, y=mean))+
  geom_bar(aes(fill=model_formula), stat = "identity")+
  geom_segment(aes(yend = mean+stderr,xend = model_formula))+
  coord_flip()+theme_bw()+
  guides(fill="none")+#guides(fill=guide_legend("model"))+
  xlab("")+ylab("mean adj. R-squared")+
  scale_x_discrete(labels = c("LLEA"="LLEA", "LLEA_PC1_10"="LLEA + \n10 human PCs"))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_LLEA_vs_LLEA10PCs_adjRsquared_barplot.png",
       dpi = 300, device="png")

#### 21km FDR ####
LLEA_lm_results<-mutate(LLEA_lm_results, model_formula = factor(model_formula, levels = c("LLEA","LLEA_PC1",paste0("LLEA_PC1_",2:10))))
LLEA_lm_results %>%
 group_by(model_formula) %>%
  count()
#n = 138268
snp_CHROM<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv") %>% select("#CHROM")
snp_bp<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv") %>% select("POS")
#the SNPs are reported in the same order as they were listed in the genotype file, so it should just be the meta info for each snp from the genotype file repeated 11 times
LLEA_lm_results <- add_column(LLEA_lm_results, Chromosome = rep(pull(snp_CHROM),11))
LLEA_lm_results$Chromosome<-factor(LLEA_lm_results$Chromosome, levels = c("1","2","3","4","5","6","7","8","9","10","0"))
LLEA_lm_results <- add_column(LLEA_lm_results, bp = rep(pull(snp_bp),11))

filter(LLEA_lm_results, p.value <= 0.05/138268)
filter(LLEA_lm_results, p.value <= 0.05/nrow(LLEA_lm_results))

filter(LLEA_lm_results, p.value <= 0.05/nrow(LLEA_lm_results)) %>%
  group_by(model_formula) %>%
  count() %>% 
  ggplot(aes(x=model_formula, y=n))+
  geom_line(group=1)+
  geom_point()+
  theme_bw()+xlab("Model Formula")+ylab("Count of Signif. SNPs post-FDR correction")

filter(LLEA_lm_results, p.value <= 0.05/nrow(LLEA_lm_results)) %>%
  group_by(model_formula, Chromosome) %>% count() %>% 
  ggplot(aes(x=model_formula, y=n))+
  geom_line(aes(group = Chromosome, color=Chromosome))+
  geom_point(aes(color = Chromosome))+
  scale_color_viridis_d(option = "turbo")+
  theme_bw()+xlab("Model Formula")+ylab("Count of Signif. SNPs post-FDR correction")

manhattan_objects<-function(df,model_formula_string){
  manhattan_plot_rsq<-filter(df, model_formula == model_formula_string) %>%
    group_by(Chromosome) %>% summarize(chr_len = max(bp)) %>% # compute chromosome size
    mutate(tot = cumsum(chr_len)-chr_len) %>% select(-chr_len) %>% # calculate cumulative position of each chr
    left_join(filter(df,model_formula == model_formula_string), ., by= c("Chromosome"="Chromosome")) %>% #add to initial dataset
    arrange(Chromosome, bp) %>% mutate(BPcum = bp+tot) #add cumulative position of each SNP
  manhattan_plot_rsq$Chromosome<-factor(manhattan_plot_rsq$Chromosome, levels = 
                                          c("1","2","3","4","5","6","7","8","9","10","0"))
  return(manhattan_plot_rsq)
}
bon_21km_top_manhattan<-manhattan_objects(LLEA_lm_results, "LLEA")
axisdf = bon_21km_top_manhattan %>% group_by(Chromosome) %>% summarize(center = (max(BPcum) + min(BPcum) ) / 2)

ggplot(bon_21km_top_manhattan, aes(x=BPcum, y= -log(p.value))) +
  geom_point(aes(color = Chromosome, alpha = 0.8))+
  #geom_point(data=filter(bon_21km_top_manhattan, p.value <= 0.05/138268),
  #           aes(x=BPcum, y=-log(p.value)),
  #           color="red")+ #this highlights bonferroni corrected snps
  geom_hline(aes(yintercept = -log(0.05/138268)), color = "red", linetype = 2)+
  geom_hline(aes(yintercept = -log(0.05/(138268*11))), color = "blue", linetype = 2)+
  scale_color_manual(values = c(rep(c("#000000","#999999"), 5),"darkgrey"))+
  scale_x_continuous(label = axisdf$Chromosome, breaks = axisdf$center)+
  #scale_y_continuous(expand = c(0,0),limits = c(-0.03,1))+ #this will remove points with -adj r-squared values (43846 of them)
  theme_bw()+xlab("")+ylab("-log(p.value)")+
  ggtitle("LLEA p-values, 21km")+
  theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

bon_21km_LLEA10_manhattan<-manhattan_objects(LLEA_lm_results, "LLEA_PC1_10")
axisdf = bon_21km_LLEA10_manhattan %>% group_by(Chromosome) %>% summarize(center = (max(BPcum) + min(BPcum) ) / 2)

ggplot(bon_21km_LLEA10_manhattan, aes(x=BPcum, y= -log(p.value))) +
  geom_point(aes(color = Chromosome, alpha = 0.8))+
  #geom_point(data=filter(bon_21km_top_manhattan, p.value <= 0.05/138268),
  #           aes(x=BPcum, y=-log(p.value)),
  #           color="red")+ #this highlights bonferroni corrected snps
  geom_hline(aes(yintercept = -log(0.05/138268)), color = "red", linetype = 2)+
  geom_hline(aes(yintercept = -log(0.05/(138268*11))), color = "blue", linetype = 1)+
  scale_color_manual(values = c(rep(c("#000000","#999999"), 5),"darkgrey"))+
  scale_x_continuous(label = axisdf$Chromosome, breaks = axisdf$center)+
  #scale_y_continuous(expand = c(0,0),limits = c(-0.03,1))+ #this will remove points with -adj r-squared values (43846 of them)
  theme_bw()+xlab("")+ylab("-log(p.value)")+
  ggtitle("LLEA+10 PCs p-values, 21km")+
  theme(legend.position = "none", panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())


#### LM with admixture 100Km####
LLEA_100km_lm_results<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/100km_maize/100km_raw_LDprune_MAF0.01_LatLongElevAdmix_lm_results.tsv")
LLEA_100km_lm_results<-add_column(LLEA_100km_lm_results, model_formula = "LLEA")

for(i in c("PC1","PC1_2","PC1_3","PC1_4","PC1_5","PC1_6","PC1_7","PC1_8","PC1_9","PC1_10")){
  LLEA_100km_lm_results<-add_row(LLEA_100km_lm_results, 
                           read_tsv(paste0("/group/jrigrp11/snodgras_maizePopulations/100km_maize/100km_raw_LDprune_MAF0.01_LatLongElevAdmix",i,"_lm_results.tsv")),
                           model_formula = paste0("LLEA_",i))
}
summary_LLEA_100km_lm<-LLEA_100km_lm_results %>% group_by(model_formula) %>%
  summarize(mean = mean(adj.r.squared, na.rm = T),
            sd = sd(adj.r.squared, na.rm = T),
            min = min(adj.r.squared),
            max = max(adj.r.squared),
            stderr = sd / sqrt(nrow(filter(LLEA_100km_lm_results, model_formula == "LLEA")))
  )

summary_LLEA_100km_lm<-mutate(summary_LLEA_100km_lm, 
                        model_formula = factor(model_formula, levels = c("LLEA","LLEA_PC1",paste0("LLEA_PC1_",2:10))))

ggplot(summary_LLEA_100km_lm, aes(x=model_formula,y=mean))+
  geom_line(group=1)+
  geom_line(aes(x=c(1:11),
                y=c(Mexmaize_eigenval$value[1],cumsum(Mexmaize_eigenval$value[1:10]))),
            group=1,
            linetype=2)+
  geom_point()+
  geom_point(aes(x=c("LLEA","LLEA_PC1",paste0("LLEA_PC1_",2:10)),
                 y=c(Mexmaize_eigenval$value[1],cumsum(Mexmaize_eigenval$value[1:10]))),
             shape=2)+
  geom_segment(aes(x=model_formula, xend=model_formula,
                   y=mean-3*stderr, yend=mean+3*stderr))+
  theme_bw()+
  ylab("Adj. R-squared")+xlab("Model Parameters")

summary_LLEA_100km_lm %>% filter(model_formula %in% c("LLEA","LLEA_PC1_10")) %>%
  ggplot(aes(x=model_formula, y=mean))+
  geom_bar(aes(fill=model_formula), stat = "identity")+
  geom_segment(aes(yend = mean+stderr,xend = model_formula))+
  coord_flip()+theme_bw()+
  guides(fill="none")+#guides(fill=guide_legend("model"))+
  xlab("")+ylab("mean adj. R-squared")+
  scale_x_discrete(labels = c("LLEA"="LLEA", "LLEA_PC1_10"="LLEA + \n10 human PCs"))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/100km_LLEA_vs_LLEA10PCs_adjRsquared_barplot.png",
       dpi = 300, device="png")

### Trying Jeff's plot

PVE_LLEA_100km<-0.0144*100
PVE_LLEAPC1_10_100km<-0.0170*100
PVE_difference_100km<- PVE_LLEAPC1_10_100km-PVE_LLEA_100km 
PVE_maize10PCs<-sum(Mexmaize_eigenval$value[1:10])*100

PVE_LLEA<-0.0133*100 #for the 21km results
PVE_LLEAPC1_10<-0.0162*100
PVE_difference<- PVE_LLEAPC1_10-PVE_LLEA 


tibble(lm_PVEs = c(PVE_LLEA_100km, PVE_difference_100km), 
       maize_PVEs = c(PVE_maize10PCs, PVE_maize10PCs),
       lm_type = c("LLEA","Human PCs")) %>%
  ggplot(aes(x=c("Maize PCs","lm PVEs")))+
  geom_bar(stat="identity",position="stack",
           fill="#999999",
           data = tibble(maize_PVE=Mexmaize_eigenval$value[1:10]*100,
                         PC= factor(1:10, levels = c(10:1))),
           aes(y=maize_PVE, x="Maize PCs", color = PC))+
  scale_color_manual(values = c("1"="#000000","2"="#000000","3"="#000000","4"="#000000","5"="#000000","6"="#000000","7"="#000000","8"="#000000","9"="#000000","10"="#000000"))+
  scale_fill_manual(values = c("LLEA"="#2C8C99",
                               "Human PCs"="#42D9C8"))+
  geom_bar(stat = "identity", position = "stack",
           aes(y=lm_PVEs, x="lm PVEs", fill=lm_type))+
  theme_bw()+
  guides(fill = guide_legend(title = "Model: ", position = "bottom"),
         color="none")+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold"))+
  xlab("")+ylab("Percent Variance Explained")+
  ggtitle("100km maize samples")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/100km_lm_pve_vs_maizePCs_barplot.png",
       width = 5, height = 5)

#7C9EB2, #FFBA49
#2C8C99, #42D9C8
#F0A202, #748E54
#028090, #00A896

tibble(lm_PVEs = c(PVE_LLEA, PVE_difference), 
       maize_PVEs = c(PVE_maize10PCs, PVE_maize10PCs),
       lm_type = c("LLEA","Human PCs")) %>%
  ggplot(aes(x=c("top 10 maize PCs","model")))+
  geom_bar(stat="identity",position="stack",
           fill="#999999",
           data = tibble(maize_PVE=Mexmaize_eigenval$value[1:10]*100,
                         PC= factor(1:10, levels = c(10:1))),
           aes(y=maize_PVE, x="top 10 maize PCs", color = PC))+
  scale_color_manual(values = c("1"="#000000","2"="#000000","3"="#000000","4"="#000000","5"="#000000","6"="#000000","7"="#000000","8"="#000000","9"="#000000","10"="#000000"))+
  scale_fill_manual(values = c("LLEA"="#2C8C99",
                               "Human PCs"="#42D9C8"))+
  geom_bar(stat = "identity", position = "stack",
           aes(y=lm_PVEs, x="model", fill=lm_type))+
  theme_bw()+
  guides(fill = guide_legend(title = "Model: ", position = "bottom"),
         color="none")+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold"))+
  xlab("")+ylab("Percent Variance Explained")+
  ggtitle("21km maize samples")+
  guides(color="none")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/21km_lm_pve_vs_maizePCs_barplot.png",
       height=5,width=5, units = "in",dpi=300)

#plotting both 21km and 100km results on same plot
tibble(lm_PVEs = c(PVE_LLEA_100km, PVE_difference_100km,
                   PVE_LLEA, PVE_difference), 
       maize_PVEs = c(PVE_maize10PCs, PVE_maize10PCs,
                      PVE_maize10PCs, PVE_maize10PCs),
       lm_type = c("LLEA","Human PCs",
                   "LLEA","Human PCs"),
       model_type = c("100km","100km",
                      "21km","21km")) %>%
  ggplot(aes(x=c("Maize PCs","100km", "21km")))+
  geom_bar(stat="identity",position="stack",
           fill="#999999",
           data = tibble(maize_PVE=Mexmaize_eigenval$value[1:10]*100,
                         PC= factor(1:10, levels = c(10:1))),
           aes(y=maize_PVE, x="Maize PCs", color = PC))+
  scale_color_manual(values = c("1"="#000000","2"="#000000","3"="#000000","4"="#000000","5"="#000000","6"="#000000","7"="#000000","8"="#000000","9"="#000000","10"="#000000"))+
  scale_fill_manual(values = c("LLEA"="#2C8C99",
                               "Human PCs"="#42D9C8"))+
  geom_bar(stat = "identity", position = "stack",
           aes(y=lm_PVEs, x=model_type, fill=lm_type))+
  theme_bw()+
  guides(fill = guide_legend(title = "Model: ", position = "bottom"),
         color="none")+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.title = element_text(size=16,face="bold"))+
  xlab("")+ylab("Percent Variance Explained")
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/100km_vs_21km_lm_pve_vs_maizePCs_barplot.pdf",
       dpi=300, height=5, width = 6, device = "pdf")
