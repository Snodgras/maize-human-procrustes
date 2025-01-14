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
#plot pca
ggplot(human, aes(x=PC_1,y=PC_2))+
  geom_point(aes(color = Region))+
  theme_bw()+
  xlab("PC1")+ylab("PC2")+
  ggtitle("Indigenous human PCA")

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

#try projecting the longitude and latitude using the Gall-Peters equal area cylindrical projections
#to convert to radians, multiple the degrees by pi/180
#if in radians
#x=radius of earth*longitude
#y=2*radius of earth*sin(latitude)
#using radius of globe = 3,958.8 miles
human<-mutate(human, 
              proj.longitude = longitude*(pi/180),
              proj.longitude = (3958.8*proj.longitude),
              proj.latitude = latitude*(pi/180),
              proj.latitude = 2*3958.8*sin(proj.latitude))
ggplot(human, aes(x=proj.longitude,y=proj.latitude))+
  geom_point(aes(color = Region))+
  theme_bw()+
  xlab("Longitude (radians)")+ylab("Latitude (radians)")+
  ggtitle("Indigenous human location, Gall-Peters projection")

#try procrustes on the projections
human_procrustes.proj<-procrustes(Y = select(human_noOutliers, c("PC_1","PC_2")),
                             X = select(human_noOutliers, c("proj.longitude","proj.latitude"))
)
plot(human_procrustes.proj)
protest(X = select(human_noOutliers, c("PC_1","PC_2")),
        Y = select(human_noOutliers, c("proj.longitude","proj.latitude")))

#human custom plots:
#numbers specified here should be number of samples
human.translation<-matrix(rep(human_procrustes.proj$translation,813),nrow=813,ncol=2,byrow =TRUE) 

human_noOutliers<-add_column(human_noOutliers, 
                             transformed_PC1 = (human_procrustes.proj$scale*t((human_procrustes.proj$rotation) %*% t(as.matrix(human_noOutliers[,c("PC_1","PC_2")])))+human.translation)[,1],
           transformed_PC2 = (human_procrustes.proj$scale*t((human_procrustes.proj$rotation) %*% t(as.matrix(human_noOutliers[,c("PC_1","PC_2")])))+human.translation)[,2])
ggplot(human_noOutliers, aes(x=transformed_PC1,y=transformed_PC2))+
  geom_point(aes(color=Region))+
  theme_bw()

human_procrustes.proj$Yrot %>% cbind(.,human_noOutliers$Region) %>% as.data.frame() %>% 
  ggplot(aes(x=as.numeric(V1), y=as.numeric(V2), color = V3))+ 
  geom_point()

human_noOutliers<-human_noOutliers %>% add_column(as.data.frame(human_procrustes$Yrot))
colnames(human_noOutliers)[24:25]<-c("rotated_PC1","rotated_PC2")
human_noOutliers<-human_noOutliers %>%mutate(rotated_PC1 = rotated_PC1+human_procrustes$translation[,1],
                                             rotated_PC2 = rotated_PC2+human_procrustes$translation[,2])
human_noOutliers%>% 
  ggplot(aes(x=rotated_PC1, y=rotated_PC2, color = Region))+ 
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

#Mexican subset of humans
Mexhuman_pca<-read_csv("/group/jrigrp11/snodgras_maizePopulations/human_PCs/all-mexican.pca.csv")
Mexhuman<-inner_join(Mexhuman_pca, human_meta, by = c("sample_name" = "sample_id" ))

#plot pca
ggplot(Mexhuman, aes(x=PC_1,y=PC_2))+
  geom_point(aes(color = Region))+
  theme_bw()+
  xlab("PC1")+ylab("PC2")+
  ggtitle("Mexican Indigenous human PCA")

Mexhuman<-mutate(Mexhuman, 
              proj.longitude = longitude*(pi/180),
              proj.longitude = (3958.8*proj.longitude),
              proj.latitude = latitude*(pi/180),
              proj.latitude = 2*3958.8*sin(proj.latitude))
ggplot(Mexhuman, aes(x=proj.longitude,y=proj.latitude))+
  geom_point(aes(color = Region))+
  theme_bw()+
  xlab("Longitude (radians)")+ylab("Latitude (radians)")+
  ggtitle("Indigenous Mexican human location, Gall-Peters projection")

Mexhuman_procrustes.proj<-procrustes(Y = select(Mexhuman, c("PC_1","PC_2")),
                                  X = select(Mexhuman, c("proj.longitude","proj.latitude"))
)
plot(Mexhuman_procrustes.proj)
protest(Y = select(Mexhuman, c("PC_1","PC_2")),
           X = select(Mexhuman, c("proj.longitude","proj.latitude")))

####Maize data####
#read in data from plink
maize_pca<-read_table2("/group/jrigrp11/snodgras_maizePopulations/plink/maize.eigenvec",col_names = FALSE)
maize_eigenval<-scan("/group/jrigrp11/snodgras_maizePopulations/plink/maize.eigenval")
#remove nuisance column for the double ID from plink
maize_pca<-maize_pca[,-1]
#set names
names(maize_pca)[1]<-"ind"
names(maize_pca)[2:ncol(maize_pca)]<- paste0("PC",1:(ncol(maize_pca)-1))

#convert eigenvalues to % variance explained
pve<- data.frame(PC = 1:20, pve = maize_eigenval/sum(maize_eigenval)*100)
#plot percent variance explained
ggplot(pve, aes(PC, pve)) +geom_bar(stat = "identity")+ylab("Percentage variance explained")+theme_bw()
#to caluclate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

#plot pca
ggplot(maize_pca, aes(x=PC1, y=PC2))+
  geom_point()+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(pve$pve[1],3),"%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2],3),"%)"))+
  ggtitle("Maize PCA")

#PCA from mean centered
maize_pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/filtered_centeredMaizeGBS.pca")
maize_eigenval<-read_lines("/group/jrigrp11/snodgras_maizePopulations/filtered_centeredMaizeGBS.eigenval") %>% as.numeric()

pve<- data.frame(PC = 1:length(maize_eigenval), pve = (maize_eigenval/sum(maize_eigenval))*100)

ggplot(maize_pca, aes(x=PC1, y=PC2))+
  geom_point()+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(pve$pve[1],3),"%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2],3),"%)"))+
  ggtitle("Maize PCA (after mean centering)")

#load in the maize location data
maize_meta<-read_csv("/group/jrigrp11/snodgras_maizePopulations/clean_seedspassport.csv")
ggplot(maize_meta, aes(x=locations_longitude,y=locations_latitude))+
  geom_point()+
  theme_bw()+
  xlab("Longitude")+ylab("Latitude")+
  ggtitle("Maize Location")

colnames(maize_pca)[1]<-"ind"

#join meta and pcs
maize<-mutate(maize_pca, sample_id = str_split(ind, ":",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

#procrustes
#compute procrustes using vegan 2.6-4
maize_procrustes<-procrustes(X = select(maize_noOutliers, c("PC1","PC2")),
                             Y = select(maize_noOutliers, c("locations_longitude","locations_latitude"))
)
summary(maize_procrustes)
plot(maize_procrustes)
plot(maize_procrustes, to.target = F, ar.col = NA)
#test the significance of two configurations
protest(X = select(maize, c("PC1","PC2")),
        Y = select(maize, c("locations_longitude","locations_latitude")))

#project the maize lat and long using gall-peters
maize<-mutate(maize, 
              proj.longitude = locations_longitude*(pi/180),
              proj.longitude = (3958.8*proj.longitude),
              proj.latitude = locations_latitude*(pi/180),
              proj.latitude = 2*3958.8*sin(proj.latitude))

ggplot(maize, aes(x=proj.longitude,y=proj.latitude))+
  geom_point(aes(color = countries_country_name))+
  theme_bw()+
  xlab("Longitude (radians)")+ylab("Latitude (radians)")+
  ggtitle("Maize Location, Gall-Peters Projection")

ggplot(maize, aes(x=PC1, y=PC2))+
  geom_point(aes(color = countries_country_name))+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(pve$pve[1],3),"%)"))+
  ylab(paste0("PC2 (", signif(pve$pve[2],3),"%)"))+
  ggtitle("Maize PCA (after mean centering)")

#compute a projected procrustes for maize
maize_procrustes.proj<-procrustes(Y = select(maize, c("PC1","PC2")),
                             X = select(maize, c("proj.longitude","proj.latitude"))
)
plot(maize_procrustes.proj)

maize.protest<-protest(Y = select(maize, c("PC1","PC2")),
        X = select(maize, c("proj.longitude","proj.latitude")))

my.translation<-matrix(rep(maize_procrustes.proj$translation,2899),nrow=2899,ncol=2,byrow =TRUE) 

plot(maize_procrustes.proj$scale*t((maize_procrustes.proj$rotation) %*% t(as.matrix(maize[,c("PC1","PC2")])))+my.translation)

maize<-add_column(maize, transformed_PC1 = (maize_procrustes.proj$scale*t((maize_procrustes.proj$rotation) %*% t(as.matrix(maize[,c("PC1","PC2")])))+my.translation)[,1],
                  transformed_PC2 = (maize_procrustes.proj$scale*t((maize_procrustes.proj$rotation) %*% t(as.matrix(maize[,c("PC1","PC2")])))+my.translation)[,2])
ggplot(maize, aes(x=transformed_PC1,y=transformed_PC2))+
  geom_point(aes(color = countries_country_name))+
  theme_bw()

#Maize subset of Mexico

Mexmaize_pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/Mex_subset/filtered_centeredMaizeGBS.pca")
colnames(Mexmaize_pca)[1]<-"ind"
Mexmaize_eigenval<-read_lines("/group/jrigrp11/snodgras_maizePopulations/Mex_subset/filtered_centeredMaizeGBS.eigenval") %>% as.numeric()

Mexpve<- data.frame(PC = 1:length(Mexmaize_eigenval), pve = (Mexmaize_eigenval/sum(Mexmaize_eigenval))*100)

ggplot(Mexmaize_pca, aes(x=PC1, y=PC2))+
  geom_point()+
  coord_equal()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(Mexpve$pve[1],3),"%)"))+
  ylab(paste0("PC2 (", signif(Mexpve$pve[2],3),"%)"))+
  ggtitle("Mexican Maize PCA (after mean centering)")

filter(maize, countries_country_name == "MEXICO") %>%
ggplot(aes(x=proj.longitude,y=proj.latitude))+
  geom_point(aes(color = locations_elevation))+
  theme_bw()+
  xlab("Longitude (radians)")+ylab("Latitude (radians)")+
  scale_color_gradient(low="orange", high = "blue")
  ggtitle("Mexican Maize Location, Gall-Peters Projection")

#join meta and pcs
Mexmaize<-mutate(Mexmaize_pca, sample_id = str_split(ind, ":",simplify = T)[,1]) %>%
  inner_join(y=maize_meta, by = c("sample_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

Mexmaize<-mutate(Mexmaize, 
              proj.longitude = locations_longitude*(pi/180),
              proj.longitude = (3958.8*proj.longitude),
              proj.latitude = locations_latitude*(pi/180),
              proj.latitude = 2*3958.8*sin(proj.latitude))

Mexmaize_procrustes.proj<-procrustes(Y = select(Mexmaize, c("PC1","PC2")),
                                  X = select(Mexmaize, c("proj.longitude","proj.latitude"))
)
plot(Mexmaize_procrustes.proj)

Mexmaize.protest<-protest(Y = select(Mexmaize, c("PC1","PC2")),
                       X = select(Mexmaize, c("proj.longitude","proj.latitude")))

ggplot(Mexmaize, aes(x=PC1, y=PC2))+
  geom_point(aes(color = locations_elevation, alpha = proj.longitude))+
  theme_bw()+
  scale_color_gradient(low="orange",high = "blue")+
  ggtitle("Mexican maize PCs by elevation")

####Maize map####
maize.outlier.stats<-maize %>% group_by(countries_country_name) %>% summarize(mean.PC1 = mean(PC1, na.rm = T),
                                                              mean.PC2 = mean(PC2, na.rm = T),
                                                              sd.PC1 = sd(PC1, na.rm = T),
                                                              sd.PC2 = sd(PC2, na.rm = T))
#filter out any samples where they're beyond 3 sd from mean for their regional group
maize.outlier<-inner_join(maize, maize.outlier.stats, by = "countries_country_name") %>% 
  mutate(PC1_lower.cutoff = mean.PC1 - 3*sd.PC1, 
         PC1_upper.cutoff = mean.PC1 + 3*sd.PC1,
         PC2_lower.cutoff = mean.PC2 - 3*sd.PC2, 
         PC2_upper.cutoff = mean.PC2 + 3*sd.PC2) %>%
  filter(PC1 < PC1_lower.cutoff | PC1 > PC1_upper.cutoff | PC2 < PC2_lower.cutoff | PC2 > PC2_upper.cutoff)

maize.outlier %>% ggplot(aes(x=locations_longitude, y=locations_latitude, color = countries_country_name))+geom_point()
maize.outlier %>% ggplot(aes(x=PC1, y=PC2, color = countries_country_name))+geom_point()

maize_noOutliers<-filter(maize, !ind %in% maize.outlier$ind)
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

#### Pair humans with maize PC centroids based on 20km buffer ####
ind_human_maize_20km_buffer<-read_csv("/group/jrigrp11/snodgras_maizePopulations/Individual-human-20km-maize-samples-joined.csv")
colnames(ind_human_maize_20km_buffer)<-c("sample_id","Ethnicity","Country.human","Region","latitude","longitude","Indigenous_Proportion","id","general_identifier",
                                         "name","bank_number","taxonomy_id","collnumb","colldate","location_id","created_on","locations_id","locations_region",
                                         "locations_site_name","locations_elevation","locations_latitude","locations_longitude","countries_id","countries_country_code2",
                                         "countries_country_code3","countries_country_name","taxonomies_id","taxonomies_genus","taxonomies_species","taxonomies_crop_name",
                                         "taxonomies_ploidy","Sample_ID_of_DNA_from_composite_samples","Sample_ID_of_DNA_from_most_recent_CML_regenerations",
                                         "Sample_ID_of_DNA_from_single_plants_used_in_GWAS","Tester_GID","Tester_pedigree","Testcross_GID","PrimaryRace","PrimaryPurity",
                                         "SecondaryRace","Pedigree","GrainType1","GrainType2","GrainType3","GrainColor1","GrainColor2","GrainColor3","PopulationType")
ind_human_maize_20km_buffer %>% group_by(sample_id) %>% count() %>% 
  ggplot(aes(x=n))+
  geom_histogram(binwidth = 1)

#write out the ids so that pca can be rerun on just those samples
selected_passport_ids<-read_csv("/group/jrigrp11/snodgras_maizePopulations/selected_genotypeIDs.csv", col_names = c("vcf_sample_id","gwas_sample_id"))
filter(selected_passport_ids, gwas_sample_id %in% ind_human_maize_20km_buffer$Sample_ID_of_DNA_from_single_plants_used_in_GWAS) %>% 
  write_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_maize_samples_ids.tsv")

#after running pca
ind_human_maize_20km_buffer.pca<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/21km_maize/filtered_centeredMaizeGBS.pca")
colnames(ind_human_maize_20km_buffer.pca)[1]<-"ind"
ind_human_maize_20km_buffer.eigenval<-read_lines("/group/jrigrp11/snodgras_maizePopulations/21km_maize/filtered_centeredMaizeGBS.eigenval") %>% as.numeric()

km21_pve<- data.frame(PC = 1:length(ind_human_maize_20km_buffer.eigenval), pve = (ind_human_maize_20km_buffer.eigenval/sum(ind_human_maize_20km_buffer.eigenval))*100)

ggplot(ind_human_maize_20km_buffer.pca, aes(x=PC1, y=PC2))+
  geom_point()+
  theme_bw()+
  xlab(paste0("PC1 (", signif(km21_pve$pve[1],3),"%)"))+
  ylab(paste0("PC2 (", signif(km21_pve$pve[2],3),"%)"))+
  ggtitle("21km radius maize PCA")

ind_human_maize_20km_buffer<-left_join(x=selected_passport_ids, y=ind_human_maize_20km_buffer.pca, by=c("vcf_sample_id"="ind")) %>% 
  right_join(x=., y=ind_human_maize_20km_buffer, by=c("gwas_sample_id"="Sample_ID_of_DNA_from_single_plants_used_in_GWAS"))

ggplot(ind_human_maize_20km_buffer, aes(x=PC1, y=PC2))+
  geom_point(aes(color=locations_elevation))+theme_bw()+
  scale_color_viridis_c(name="Elevation")+
  xlab(paste0("PC1 (", signif(km21_pve$pve[1],3),"%)"))+
  ylab(paste0("PC2 (", signif(km21_pve$pve[2],3),"%)"))+
  ggtitle("21km radius maize PCA")

#Goal is to find the centroid of PCs for each set of maize samples for a given human sample
centroidPCs_ind_human<-ind_human_maize_20km_buffer %>% 
  group_by(sample_id,Region,Ethnicity,latitude,longitude) %>% 
  summarize(across(contains("PC"), ~ mean(.x,na.rm = T)),
            mean_elevation.maize=mean(locations_elevation))
colnames(centroidPCs_ind_human)<-c("sample_name","Region","Ethnicity","latitude","longitude",
                                   paste0("PC",1:20,".maize"), "mean_elevation.maize")
#add in the human PCs
centroidPCs_ind_human<-inner_join(centroidPCs_ind_human, human_pcs, by=("sample_name"))

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

####FOR SOME REASON THE SCALE IS WAY OFF FOR THE MAIZE PCS COMPARED TO HUMAN PCS####
#procrustes of the pair to get the maize (centroid) PCs onto the human PC coordinates
centroidPCs_ind_human.protest<-protest(X=select(ungroup(centroidPCs_ind_human), PC_1,PC_2),
                                          Y=select(ungroup(centroidPCs_ind_human), PC1.maize,PC2.maize))
#ss=0.7044, correlation=0.5437, significance=0.001
test.protest<-protest(Y=select(ungroup(centroidPCs_ind_human), PC_1,PC_2),
                      X=select(ungroup(centroidPCs_ind_human), PC1.maize,PC2.maize))

centroidPCs_ind_human<-centroidPCs_ind_human %>% add_column(as.data.frame(centroidPCs_ind_human.protest$Yrot))
colnames(centroidPCs_ind_human)[37:38]<-c("rotated_PC1","rotated_PC2")
centroidPCs_ind_human<-centroidPCs_ind_human %>%mutate(rotated_PC1 = (rotated_PC1*centroidPCs_ind_human.protest$scale)+centroidPCs_ind_human.protest$translation[,1],
                                             rotated_PC2 = (rotated_PC2*centroidPCs_ind_human.protest$scale)+centroidPCs_ind_human.protest$translation[,2])

centroidPCs_ind_human<-centroidPCs_ind_human %>% add_column(as.data.frame(test.protest$Yrot))
colnames(centroidPCs_ind_human)[39:40]<-c("rotated_PC1.human","rotated_PC2.human")
centroidPCs_ind_human<-centroidPCs_ind_human %>%mutate(rotated_PC1.human = (rotated_PC1*test.protest$scale)+test.protest$translation[,1],
                                                       rotated_PC2.human = (rotated_PC2*test.protest$scale)+test.protest$translation[,2])

#plot the maize (centroid, procrustes) PCs and human PCs on the human PC coordinates
ggplot(centroidPCs_ind_human )+
  geom_point(aes(x=PC_1,y=PC_2, color=Region), shape=17)+
  geom_point(aes(x=rotated_PC1,y=rotated_PC2), shape = 1, color="goldenrod")+
  theme_bw()+xlab("human PC1")+ylab("human PC2")+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1aff",
               "#b100ea","#386955","#a7c957"),
    name = "Region", labels = c(
      "Mexico, Center","Mexico, Gulf", 
      "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
      "Mexico, Oaxaca","Mexico, West"))

ggplot(centroidPCs_ind_human,aes(color=Region) )+
  geom_point(aes(x=PC1.maize,y=PC1.maize), shape=1)+
  geom_point(aes(x=rotated_PC1.human,y=rotated_PC2.human), shape = 17)+
  theme_bw()+xlab("maize PC1")+ylab("maize PC2")

#### maybe the centroids are what's making the graphing weird
paired_procrustes<- protest(X=select(ungroup(inner_join(human_pcs, ind_human_maize_20km_buffer, by=c("sample_name"="sample_id"))),PC_1,PC_2), 
                            Y=select(ungroup(inner_join(human_pcs, ind_human_maize_20km_buffer, by=c("sample_name"="sample_id"))),PC1,PC2))
paired_procrustes #sum of squares = 0.8822, correlation = 0.3432, significance = 0.001

ungroup(inner_join(human_pcs, ind_human_maize_20km_buffer, by=c("sample_name"="sample_id"))) %>%
  add_column(as.data.frame(paired_procrustes$Yrot)) %>%
  mutate(rotated_PC1.maize = (V1*paired_procrustes$scale)+paired_procrustes$translation[,1],
         rotated_PC2.maize = (V2*paired_procrustes$scale)+paired_procrustes$translation[,2]) %>%
  ggplot(aes(color=Region))+
  geom_point(aes(x=PC_1,y=PC_2), shape=17)+
  geom_point(aes(x=rotated_PC1.maize,y=rotated_PC2.maize), shape = 1)+
  theme_bw()+xlab("human PC1")+ylab("human PC2")+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1aff",
               "#b100ea","#386955","#a7c957"),
    name = "Region", labels = c(
      "Mexico, Center","Mexico, Gulf", 
      "Mexico, Mayan","Mexico, North of Mesoamerica","Mexico, North",
      "Mexico, Oaxaca","Mexico, West"))

#writing the paired human PCs file for PCcorrelations.R
ungroup(inner_join(human_pcs, ind_human_maize_20km_buffer, by=c("sample_name"="sample_id"))) %>% 
  select(c(paste0("PC_",1:10),sample_name,vcf_sample_id)) %>% select(vcf_sample_id) %>% unique()


#### Paired points for procrustes of maize and human PCs ####
paired_output<-read_csv("/group/jrigrp11/snodgras_maizePopulations/anchored_procrustes/anchor.0.out.csv")
for(i in 1:97){
  paired_output<-rbind(paired_output, read_csv(paste0("/group/jrigrp11/snodgras_maizePopulations/anchored_procrustes/anchor.",i,".out.csv")))
}
#will need to add in 98 and 99, they're taking a while to run

ggplot(paired_output)+
  geom_histogram(aes(x=SumOfSquares), binwidth = 0.001)
ggplot(paired_output)+
  geom_histogram(aes(x = Correlation), binwidth = 0.001)

summary(paired_output)

####Correlation of PCs from paired anchors####
#create a dataframe where the first column is the iteration and second column is the dataframe (like a list?)
anchor_pairs<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/anchored_procrustes/1.anchors.tsv")

cor(anchor_pairs$pairs$PC_1.human, anchor_pairs$pairs$PC1.maize)
cor(anchor_pairs[,"PC_1.human"],anchor_pairs[,"PC1.maize"]) %>% unname()

PC_cors<-tibble(anchor_iteration = NA, 
                human_PC = NA,
                maize_PC = NA, 
                correlation = NA)
for(i in c(1:3,18:99)){ #for each iteration
  anchor_pairs<-read_tsv(paste0("/group/jrigrp11/snodgras_maizePopulations/anchored_procrustes/",
                                i,".anchors.tsv")
  )
  for(h in 1:10){ #for each human PC
    for(m in 1:20){ #for each maize PC
      PC_cors<-add_row(PC_cors,
                       anchor_iteration = i, 
                       human_PC = paste0("PC",h),
                       maize_PC = paste0("PC",m), 
                       correlation = cor(anchor_pairs[,paste0("PC_",h,".human")],
                                        anchor_pairs[,paste0("PC",m,".maize")])%>% unname()
                      )
    }
  }
}
PC_cors<-na.omit(PC_cors)

PC_cors.summary<-PC_cors %>% group_by(human_PC, maize_PC) %>% 
  summarize(cor.mean = mean(correlation,na.rm=T),
            cor.median = median(correlation,na.rm=T),
            cor.range = (max(correlation) - min(correlation)))

PC_cors.summary$human_PC<-PC_cors.summary$human_PC %>% factor(levels = c(paste0("PC",1:10)))
PC_cors.summary$maize_PC<-PC_cors.summary$maize_PC %>% factor(levels = c(paste0("PC",1:20)))

ggplot(PC_cors.summary, aes(x=human_PC, y=maize_PC))+
  #geom_tile(aes(fill = cor.mean))+
  geom_tile(aes(fill = cor.median))+
  #geom_tile(aes(fill = cor.range))+
  #scale_fill_viridis_c()+
  theme_bw()+xlab("Human PC")+ylab("Maize PC")

#try anchor points only in Mexico
ggplot(anchor_pairs, aes(x=longitude, y=latitude, color = distance))+geom_point()+theme_bw()+scale_color_viridis_c()
filter(anchor_pairs, str_detect(Region, "Mexico")) %>% 
  ggplot(aes(x=longitude, y=latitude, color = distance))+
  geom_point() +theme_bw()+scale_color_viridis_c()
filter(anchor_pairs, str_detect(Region, "Mexico")) %>% 
  ggplot(aes(x=locations_longitude, y=locations_latitude, color = distance))+
  geom_point() +theme_bw()+scale_color_viridis_c()


#example graph of maize on human PC using anchor points
anchor_pairs %>% colnames()
anchor_pairs<-inner_join(anchor_pairs, select(human_noOutliers, sample_name, latitude, longitude, Region), by = c("human_id" = "sample_name"))
anchor_pairs<-inner_join(anchor_pairs, select(maize_noOutliers, ind, locations_latitude, locations_longitude, countries_country_name, Regions))

paired_procrustes<-procrustes(X = select(anchor_pairs, "PC_1.human","PC_2.human"),
                              Y = select(anchor_pairs, "PC1.maize","PC2.maize"))
protest(X = select(anchor_pairs, "PC_1.human","PC_2.human"),
        Y = select(anchor_pairs, "PC1.maize","PC2.maize"))
#Sum of squares = 0.6505, correlation = 0.5912

anchor_pairs<-anchor_pairs %>% add_column(as.data.frame(paired_procrustes$Yrot))
colnames(anchor_pairs)[41]<-"RegionsByCountry"
colnames(anchor_pairs)[42:43]<-c("rotated_PC1","rotated_PC2")
anchor_pairs<-anchor_pairs %>%mutate(rotated_PC1 = rotated_PC1+paired_procrustes$translation[,1],
                                     rotated_PC2 = rotated_PC2+paired_procrustes$translation[,2])
paired.plot<-ggplot(anchor_pairs)+
  geom_point(aes(x=rotated_PC1, y=rotated_PC2, color = RegionsByCountry), alpha = 0.3)+
  geom_point(aes(x=PC_1.human, y=PC_2.human, color = Region),shape = 3)+
  theme_bw()+
  xlab("Human PC1")+ylab("Human PC2")+
  scale_color_manual(name = "Regions by Country",
                                values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                           "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f",
                                           "Mexico-Center_of_Mexico"="#8e79e2","Mexico-Gulf_of_Mexico"="#6447D7", 
                                           "Mexico-Mayan_region"="#4B2CC9","Mexico-North_of_Mesoamerica"="#3f25a7","Mexico-North_of_Mexico"="#1f1254",
                                           "Mexico-Oaxaca"="#190E43","Mexico-West_of_Mexico"="#130B32"),
                                labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                           "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean",
                                           "Mexico-Center_of_Mexico"="Mexico, Center","Mexico-Gulf_of_Mexico"="Mexico, Gulf", 
                                           "Mexico-Mayan_region"="Mexico, Mayan","Mexico-North_of_Mesoamerica"="Mexico, North of Mesoamerica","Mexico-North_of_Mexico"="Mexico, North",
                                           "Mexico-Oaxaca"="Mexico, Oaxaca","Mexico-West_of_Mexico"="Mexico, West"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))

centroids<-full_join(summarize(group_by(anchor_pairs, Region),mean.humanPC1 = mean(PC_1.human,na.rm=T), mean.humanPC2=mean(PC_2.human,na.rm=T)),
          summarize(group_by(anchor_pairs, RegionsByCountry),mean.maizePC1 = mean(rotated_PC1,na.rm=T), mean.maizePC2=mean(rotated_PC2,na.rm=T)),
          by=c("Region"="RegionsByCountry"))
centroids[5:11,4]<- -71.4
centroids[5:11,5]<- 17.3
centroids<-centroids %>% mutate(line.color = case_when(Region == "Amazonia"~"#0047CC",
                                                       Region == "AndeanHighland"~"#8F0200",
                                                       Region == "CentralSouthAmerica"~"#7A6600",
                                                       Region == "ChacoAmerindian"~"#177861",
                                                       Region == "Mexico-Center_of_Mexico"~"#C6BCF1",
                                                       Region == "Mexico-Gulf_of_Mexico"~"#C6BCF1",
                                                       Region == "Mexico-Mayan_region"~"#C6BCF1",
                                                       Region == "Mexico-North_of_Mesoamerica"~"#C6BCF1",
                                                       Region == "Mexico-North_of_Mexico"~"#C6BCF1",
                                                       Region == "Mexico-Oaxaca"~"#C6BCF1",
                                                       Region == "Mexico-West_of_Mexico"~"#C6BCF1",
                                                       Region == "Patagonia" ~"#7A0066",
                                                       .default = NA))
paired.plot+geom_segment(data = na.omit(centroids),aes(x=mean.humanPC1,y=mean.humanPC2,xend=mean.maizePC1,yend=mean.maizePC2))
ggplot(anchor_pairs)+
  geom_point(aes(x=rotated_PC1, y=rotated_PC2, color = RegionsByCountry), alpha = 0.1, size=1)+
  geom_point(aes(x=PC_1.human, y=PC_2.human, color = Region),shape = 3)+
  geom_segment(data = na.omit(centroids),aes(x=mean.humanPC1,y=mean.humanPC2,xend=mean.maizePC1,yend=mean.maizePC2))+
  theme_bw()+
  xlab("Human PC1")+ylab("Human PC2")+
  scale_color_manual(name = "Regions by Country",
                     values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f",
                                "Mexico-Center_of_Mexico"="#8e79e2","Mexico-Gulf_of_Mexico"="#6447D7", 
                                "Mexico-Mayan_region"="#4B2CC9","Mexico-North_of_Mesoamerica"="#3f25a7","Mexico-North_of_Mexico"="#1f1254",
                                "Mexico-Oaxaca"="#190E43","Mexico-West_of_Mexico"="#130B32"),
                     labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean",
                                "Mexico-Center_of_Mexico"="Mexico, Center","Mexico-Gulf_of_Mexico"="Mexico, Gulf", 
                                "Mexico-Mayan_region"="Mexico, Mayan","Mexico-North_of_Mesoamerica"="Mexico, North of Mesoamerica","Mexico-North_of_Mexico"="Mexico, North",
                                "Mexico-Oaxaca"="Mexico, Oaxaca","Mexico-West_of_Mexico"="Mexico, West"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/human_maize_pcProcrustes.humanemphasized.png",device="png",dpi=300,
       height = 6, width = 7.5)
ggplot(anchor_pairs)+
  geom_point(aes(x=rotated_PC1, y=rotated_PC2, color = RegionsByCountry), alpha = 0.4, size=1)+
  geom_point(aes(x=PC_1.human, y=PC_2.human, color = Region),shape = 3, alpha = 0.2)+
  geom_segment(data = na.omit(centroids),aes(x=mean.humanPC1,y=mean.humanPC2,xend=mean.maizePC1,yend=mean.maizePC2))+
  theme_bw()+
  xlab("Human PC1")+ylab("Human PC2")+
  scale_color_manual(name = "Regions by Country",
                     values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#2d1a77",
                                "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f",
                                "Mexico-Center_of_Mexico"="#8e79e2","Mexico-Gulf_of_Mexico"="#6447D7", 
                                "Mexico-Mayan_region"="#4B2CC9","Mexico-North_of_Mesoamerica"="#3f25a7","Mexico-North_of_Mexico"="#1f1254",
                                "Mexico-Oaxaca"="#190E43","Mexico-West_of_Mexico"="#130B32"),
                     labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                                "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean",
                                "Mexico-Center_of_Mexico"="Mexico, Center","Mexico-Gulf_of_Mexico"="Mexico, Gulf", 
                                "Mexico-Mayan_region"="Mexico, Mayan","Mexico-North_of_Mesoamerica"="Mexico, North of Mesoamerica","Mexico-North_of_Mexico"="Mexico, North",
                                "Mexico-Oaxaca"="Mexico, Oaxaca","Mexico-West_of_Mexico"="Mexico, West"))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave("/group/jrigrp11/snodgras_maizePopulations/Plots/human_maize_pcProcrustes.maizeemphasized.png",device="png",dpi=300,
       height = 6, width = 7.5)

#what if we assigned the maize samples the "regions" of the human sample they're paired to?
centroids.2<-summarize(group_by(anchor_pairs, Region),
                       mean.humanPC1 = mean(PC_1.human,na.rm=T), 
                       mean.humanPC2=mean(PC_2.human,na.rm=T), 
                       mean.maizePC1 = mean(rotated_PC1,na.rm=T), 
                       mean.maizePC2=mean(rotated_PC2,na.rm=T))
ggplot(anchor_pairs)+
  geom_point(aes(x=rotated_PC1, y=rotated_PC2, color = Region), alpha = 0.4, size=1)+
  geom_point(aes(x=PC_1.human, y=PC_2.human, color = Region),shape = 3, alpha = 0.2)+
  scale_color_manual(name = "Regions by Country",
                     values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab",
                                "Patagonia"="#ec00c5","Mexico-Center_of_Mexico"="#8e79e2","Mexico-Gulf_of_Mexico"="#6f2da8", 
                                "Mexico-Mayan_region"="#4B2CC9","Mexico-North_of_Mesoamerica"="#3f25a7","Mexico-North_of_Mexico"="#1f1254",
                                "Mexico-Oaxaca"="#b200ed","Mexico-West_of_Mexico"="#702963"),
                     labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian",
                                "Patagonia"="Patagonia",
                                "Mexico-Center_of_Mexico"="Mexico, Center","Mexico-Gulf_of_Mexico"="Mexico, Gulf", 
                                "Mexico-Mayan_region"="Mexico, Mayan","Mexico-North_of_Mesoamerica"="Mexico, North of Mesoamerica","Mexico-North_of_Mexico"="Mexico, North",
                                "Mexico-Oaxaca"="Mexico, Oaxaca","Mexico-West_of_Mexico"="Mexico, West"))+
  geom_segment(data = na.omit(centroids.2),aes(x=mean.humanPC1,y=mean.humanPC2,xend=mean.maizePC1,yend=mean.maizePC2, color = Region))+
  guides(color = guide_legend(override.aes = list(alpha = 1)))+
  theme_bw()+xlab("Human PC1")+ylab("Human PC2")
  

####Ploting PCs against geography####

plot_PCsByGeo<-function(df,long,lat,PC){
  plt<-ggplot(data=df, aes_string(x=long,y=lat, color=PC))+
    geom_point()+
    theme_bw()+
    scale_color_viridis_c()
  return(plt)
}
plot_PCsByGeo(human_noOutliers, "longitude", "latitude", "PC_1")

for(i in c(paste("PC",1:10,sep="_"))){
  plot_PCsByGeo(human_noOutliers, "longitude", "latitude", i)
  ggsave(
         paste0("/group/jrigrp11/snodgras_maizePopulations/Plots/human_",i,"_byGeo.png"),
         device="png",
         dpi=300
         )
}
for(i in c(paste("PC",1:20,sep=""))){
  plot_PCsByGeo(maize_noOutliers, "locations_longitude", "locations_latitude", i)
  ggsave(
    paste0("/group/jrigrp11/snodgras_maizePopulations/Plots/maize_",i,"_byGeo.png"),
    device="png",
    dpi=300
  )
}
