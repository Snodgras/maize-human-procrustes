library(tidyverse)
library(sp)
library(sf)
#library(geosphere)

#Goal: make a map of Americas where regions are colored and location data is plotted

#### Make a map with regions colored ####
world_boundaries<-st_read("world-administrative-boundaries/world-administrative-boundaries.shp")
ggplot(world_boundaries)+geom_sf()

filter(world_boundaries, continent == "Americas") %>% 
  filter(!iso3 %in% c("CAN","USA","GRL") & region != "Northern America") %>%
  ggplot()+geom_sf()

america_boundaries<-filter(world_boundaries, continent == "Americas") %>% 
  filter(!iso3 %in% c("CAN","USA","GRL") & region != "Northern America")

ggplot(america_boundaries)+geom_sf(aes(fill = region)) #shows how aes work with geom_sf

mutate(america_boundaries, country_group = case_when(iso3 %in% c("MEX","ARG") ~ "group1",
                                                     .default = "group2")) %>%
  ggplot()+geom_sf(aes(fill = country_group)) #shows how we can create new color groups

select(human, c("Country","Region")) %>% unique()

#can't color countries by region because some countries have multiple regions

#### Try adding colored points to the map ####
ggplot(america_boundaries)+
  geom_sf(fill = "white", color = "black")+
  geom_point(data = human, aes(x=longitude, y=latitude, color = Region))+
  theme(legend.position = "bottom")

#### How to show multiple humans under same lat/long? ####

#centroids?
centroid_human_locations<-select(human, c(longitude, latitude, Region)) %>% 
  group_by(Region) %>% summarize(sample_count = n(), centroid_long = mean(longitude, na.rm = T), centroid_lat = mean(latitude, na.rm = T))
ggplot(america_boundaries)+
  geom_sf(fill = "white", color = "black")+
  geom_point(data = centroid_human_locations, aes(x=centroid_long, y=centroid_lat, color = Region, size = sample_count))+
  theme(legend.position = "right")+xlab("")+ylab("")

# add size as an aes on the original map?
ggplot(america_boundaries)+
  geom_sf(fill = "white", color = "black")+
  geom_point(data = (human %>% group_by(latitude, longitude, Region) %>% count()), 
             aes(x=longitude, y=latitude, color = Region, size = n),
             alpha = 0.7)+
  theme_bw()+xlab("")+ylab("")+ggtitle("Human sampling locations")+guides(color = guide_legend(override.aes = list(alpha = 1)))+
  scale_color_manual(values = c("#538fff","#ff5b58","#ecc500","#28d2ab","#fca207","#cb4d8e","#268189","#2d1a77",
                                "#b100ea","#386651","#a7c957","#ec00c5"),
                     labels = c("Amazonia","Andean Highland","Central South America", "Chaco-Amerindian","MX, Center of Mexico","MX, Gulf of Mexico","MX, Mayan","MX, North of Mesoamerica","MX, North of Mexico","MX, Oaxaca","MX, West of Mexico","Patagonia"))
ggsave("Plots/Human_Sampling_Locations.sampleSizeHighlighted.png",device = "png",dpi=300)

# Let's try the same thing but now for maize
mutate(america_boundaries, country_fill = case_when(iso3 %in% c("BOL","CHL","PER") ~ "AndeanHighland",
                                                    iso3 %in% c("BRA") ~ "Amazonia",
                                                    iso3 %in% c("COL","ECU","GUF","SUR","VEN","GUY") ~ "CentralSouthAmerica",
                                                    iso3 %in% c("PRY") ~ "ChacoAmerindian",
                                                    iso3 %in% c("MEX") ~ "Mexico",
                                                    iso3 %in% c("ARG","URY","FLK") ~ "Patagonia",
                                                    iso3 %in% c("CRI","SLV","GTM", "NIC","PAN","HND","BLZ") ~ "Mesoamerica",
                                                    iso3 %in% c("ATG","BRB","CUB", "JAM","DOM","GRD","GLP", "HTI","MTQ","PRI","VCT","TTO","VGB","VIR") ~ "Caribbean")) %>%
  filter(!is.na(country_fill))%>%
  ggplot()+
  geom_sf(aes(fill = country_fill), color = "black")+
  geom_point(data = maize, aes(x=locations_longitude,y=locations_latitude), size =0.1, alpha = 0.8)+
  xlab("")+ylab("")+theme_bw()+ggtitle("Maize sample locations")+
  scale_fill_manual(name = "Regions by Country",
                    values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#624db7",
                               "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
                    labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
                               "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))
ggsave("Plots/Maize_sample_locations.coloredByRegions.png",device = "png",dpi=300)


mutate(america_boundaries, country_fill = case_when(iso3 %in% c("MEX") ~ "Mexico")) %>%
  filter(!is.na(country_fill))%>%
  ggplot()+
  #geom_sf(aes(fill = country_fill), color = "black", alpha = 0.5)+
  geom_sf(color = "black")+
  geom_point(data = unique(select(human_21km_pairs, Region.human, sample_id.human, latitude.human, longitude.human)),
             aes(x=longitude.human, y=latitude.human, color = Region.human), size=1)+
  geom_point(data = filter(maize_meta, Sample_ID_of_DNA_from_single_plants_used_in_GWAS %in% maize_21km_non_transformed.pca$ind_id), 
             aes(x=locations_longitude,y=locations_latitude), size =0.1)+
  xlab("")+ylab("")+theme_bw()+
  #scale_fill_manual(name = "Regions by Country",
  #                  values = c("AndeanHighland"="#ff5b58","Amazonia"="#538fff","CentralSouthAmerica"="#ecc500","ChacoAmerindian"="#28d2ab","Mexico"="#624db7",
  #                             "Patagonia"="#ec00c5","Mesoamerica"="#a7c957","Caribbean"="#802d2f"),
  #                  labels = c("AndeanHighland"="Andean Highland","Amazonia"="Amazonia","CentralSouthAmerica"="Central South America","ChacoAmerindian"="Chaco Amerindian","Mexico"="Mexico",
  #                             "Patagonia"="Patagonia","Mesoamerica"="Mesoamerica","Caribbean"="Caribbean"))+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))
ggsave("Plots/21km_sample_map.byRegion.png",device="png",dpi=300)

#### Mexico states ####
#from https://www.geoboundaries.org/countryDownloads.html
Mexico_boundaries<-st_read("geoBoundaries-MEX-ADM1-all/geoBoundaries-MEX-ADM1_simplified.shp")
Mexico_boundaries<-Mexico_boundaries %>% 
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
ggplot(Mexico_boundaries)+
  geom_sf(aes(fill = Regional_Name), color = "black")+
  #geom_point(data = filter(maize_meta, countries_country_name == "MEXICO"), 
  #           aes(x=locations_longitude,y=locations_latitude))+
  xlab("")+ylab("")+theme_bw()+
  scale_fill_manual(values = c("Baja California"="#E28AFF",
                               "Mexican Plateau"="#2d1a77",
                               "Sierra Madre Oriental"="#6E54D9",
                               "Sierra Madre Occidental and Pacific Coastal Lowlands"="#b100ea",
                               "Mesa Central"="#fca207",
                               "Gulf Coastal Plain"="#cb4d8e",
                               "Cordillera Neovolcanica"="#a7c957",
                               "Southern Highlands"="#386651",
                               "Yucatan Peninsula"="#268189"))


ggsave("Plots/Mexican_states_colored.map.png",device = "png",dpi=300, width = 11, height = 8.5)

maize_Mex_states<-read_csv("Intersection-MexicanMaize-StatesofMexico.csv")
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
ggplot(Mexico_boundaries)+
  geom_sf(aes(fill = Regional_Name), color = "black", alpha = 0.4)+
  geom_point(data = maize_Mex_states, 
             aes(x=locations_longitude,y=locations_latitude, color=Regional_Name))+
  xlab("")+ylab("")+theme_bw()

#This doesn't work well because this isn't the Mexican only PCA
inner_join(maize_admixRemoved_pca, maize_Mex_states, by=c("ind_id" = "Sample_ID_of_DNA_from_single_plants_used_in_GWAS")) %>%
  ggplot(aes(x=PC1,y=PC2, color = Regional_Name))+
  geom_point()+
  theme_bw()+
  ggtitle("Maize admix removed PCA")

