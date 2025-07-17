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


###Scratch
# PC 1: -1,0 ; 1,0
# PC 2: 0,-1 ; 0,1

axes<-matrix(c(-1,0,1,0,0,-1,0,1), nrow = 4,ncol = 2)
theta = 20
rot_mat = axes*c(cos(theta),sin(theta))
rot_axes<-rot_mat*axes
ggplot()+
  geom_point(data = as.data.frame(axes), 
               aes(x = V1, y=V2), color = "black")+
  geom_point(data = as.data.frame(rot_axes),
             aes(x=V1, y=V2), color = "red")
  #geom_segment(data = rot_axes, 
  #             aes(x = x_start, y=y_start, xend = x_end, yend = y_end), color = "red")
  
rot_mat<-matrix(c(0.99969465, 0.02471028, -0.02471028, 0.99969465),
                ncol = 2, nrow = 2)
rot_axes<-rot_mat %*% t(axes)
colnames(rot_axes)<-c(x1,x2,y1,y2)
ggplot()+
  geom_segment(data = as.data.frame(t(axes)),
               aes(x=V1, xend=V3,y=V2,yend=V4))+
  geom_segment(data = as.data.frame(rot_axes),
               aes(x=V1, xend=V3,y=V2,yend=V4),color="red")
rot_mat<-matrix(c(-0.9755890, 0.2196042, 0.2196042, 0.9755890),
                ncol = 2, nrow = 2)
rot_axes<-rot_mat %*% t(axes)
axes1.slope = (rot_axes[2,3] - rot_axes[2,1])/(rot_axes[1,3] - rot_axes[1,1])

ggplot()+
  geom_segment(data = as.data.frame(t(axes))%>% add_column(axis_name = c("PC1.orig","PC2.orig")),
               aes(x=V1, xend=V3,y=V2,yend=V4, color = axis_name))+
  geom_segment(data = (as.data.frame(rot_axes) %>% add_column(axis_name = c("PC1","PC2"))),
               aes(x=V1, xend=V3,y=V2,yend=V4, color = axis_name))+
  geom_abline(slope = axes1.slope, intercept = 0)


x_point.1 = c(-1,0)
x_point.2 = c(1,0)
y_point.1 = c(0,-1)
y_point.2 = c(0,1)

x_prime.1 = x_point.1 %*% rot_mat
x_prime.2 = x_point.2 %*% rot_mat
y_prime.1 = y_point.1 %*% rot_mat
y_prime.2 = y_point.2 %*% rot_mat

x_prime.slope = (x_prime.2[,2] - x_prime.1[,2])/(x_prime.2[,1]-x_prime.1[,1])
y_prime.slope = (y_prime.2[,2] - y_prime.1[,2])/(y_prime.2[,1]-y_prime.1[,1])

ggplot()+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_abline(slope = x_prime.slope, color = "red")+
  geom_abline(slope = y_prime.slope, color = "blue")
