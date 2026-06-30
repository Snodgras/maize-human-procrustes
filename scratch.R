lm_results<-tibble(snp_ID = NA, r.squared = NA, adj.r.squared = NA, model_formula =NA)

for(m in models){
  lm_model<-paste0("pull(t_maize_genotypes[,snp]) ~ ", m)
  print(paste0("Using this formula: ", lm_model))
        for(snp in 1:ncol(t_maize_genotypes)){
          if(str_detect(string = colnames(t_maize_genotypes)[snp], pattern = "S[0-9]")){ 
            placeholder.summary<-summary(lm(data = t_maize_genotypes, lm_model))
            lm_results<-add_row(lm_results, 
                                snp_ID = colnames(t_maize_genotypes)[snp],
                                r.squared = placeholder.summary$r.squared,
                                adj.r.squared = placeholder.summary$adj.r.squared,
                                model_formula = m)
            if(snp %% 100 == 0){print(paste("On SNP",snp,Sys.time()))}
          }
        }
}

lm_results<-na.omit(lm_results)
write_tsv(lm_results, "lm_results/21km_raw_LDprune_MAF0.01.genotypes.1000.allmodels.lm_results.tsv")

maize21km_human_lm_Rsq<-add_row(maize21km_human_lm_Rsq, snp_ID = lm_results$snp_ID, r.squared = lm_results$r.squared, adj.r.squared = lm_results$adj.r.squared, model_formula = lm_results$model_formula)

## make the linear model line of best fit plot
test_model<-lm(data = t_maize_genotypes, S1_6120151 ~ PC_1)
#https://stackoverflow.com/questions/65009244/creating-a-line-of-best-fit-in-r
test_model$coefficients
left_join(t_maize_genotypes, human_21km_pairs, by=c("sample_id"="Sample_ID_of_DNA_from_single_plants_used_in_GWAS.maize")) %>%
ggplot(aes(x=PC_1, y=S1_6120151))+
  geom_abline(slope = test_model$coefficients[2], intercept = test_model$coefficients[1])+
  geom_point(aes(color = Region.human))+
  xlab("Human PC 1")+ylab("Maize Genotype at 1 SNP")+
  scale_y_continuous(breaks = c(0,1,2))+
  scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  theme_bw()

#Just to double check plotting multiple ways
#plot(t_maize_genotypes$PC_1,t_maize_genotypes$S1_6120151)  # plots the points
#abline(a=coef(test_model)[1], b=coef(test_model)[2])

protest(X=select(centroidPCs_Mexhuman21km, longitude.human,latitude.human), 
        Y=select(centroidPCs_Mexhuman21km, PC1.maize, PC2.maize))
protest(X=select(centroidPCs_Mexhuman21km, PC1.human,PC2.human), 
        Y=select(centroidPCs_Mexhuman21km, PC1.maize, PC2.maize))

centroid_21km_procrustes<-procrustes(X=select(centroidPCs_Mexhuman21km, PC1.human,PC2.human), 
                                     Y=select(centroidPCs_Mexhuman21km, PC1.maize, PC2.maize))

tibble(Region = centroidPCs_Mexhuman21km$Region.human,
       PC1.human = centroidPCs_Mexhuman21km$PC1.human,
       PC2.human = centroidPCs_Mexhuman21km$PC2.human,
       rotated_PC1.maize = centroid_21km_procrustes$Yrot[,1],
       rotated_PC2.maize = centroid_21km_procrustes$Yrot[,2],
       PC1.maize = centroidPCs_Mexhuman21km$PC1.maize,
       PC2.maize = centroidPCs_Mexhuman21km$PC2.maize) %>%
  ggplot()+
  #geom_point(aes(x=PC1.human, y=PC2.human, color = Region), shape = 3)+
  #xlab("Human PC 1")+ylab("Human PC 2")+
  geom_point(aes(x=rotated_PC1.maize, y=rotated_PC2.maize,color = Region))+
  #geom_point(aes(x=PC1.maize, y=PC2.maize,color = Region))+
  xlab(paste0("Maize PC 1 (",signif(maize_21km_transformed.eigens[1,2]*100,3),"%)")) + 
  ylab(paste0("Maize PC 2 (",signif(maize_21km_transformed.eigens[2,2]*100,3),"%)"))+
  theme_bw()+scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))

test_rotation<-tibble(Region = centroidPCs_Mexhuman21km$Region.human,
       PC1.human = centroidPCs_Mexhuman21km$PC1.human,
       PC2.human = centroidPCs_Mexhuman21km$PC2.human,
       rotated_PC1.maize = centroid_21km_procrustes$Yrot[,1],
       rotated_PC2.maize = centroid_21km_procrustes$Yrot[,2],
       PC1.maize = centroidPCs_Mexhuman21km$PC1.maize,
       PC2.maize = centroidPCs_Mexhuman21km$PC2.maize) 

test_rotation %>% 
  mutate(PC1.maize.meancenter = PC1.maize - mean(PC1.maize,na.rm=TRUE),
         PC2.maize.meancenter = PC2.maize - mean(PC2.maize,na.rm=TRUE)) %>%
  ggplot()+
  geom_point(aes(x=rotated_PC1.maize, y=rotated_PC2.maize,color = Region))+
  geom_point(aes(x=PC1.maize.meancenter*centroid_21km_procrustes$scale, 
                 y=PC2.maize.meancenter*centroid_21km_procrustes$scale,color = Region), shape = 3)+
  xlab(paste0("Maize PC 1 (",signif(maize_21km_transformed.eigens[1,2]*100,3),"%)")) + 
  ylab(paste0("Maize PC 2 (",signif(maize_21km_transformed.eigens[2,2]*100,3),"%)"))+
  theme_bw()+scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))

axes<-tibble(x=c(-100,0),
             xend=c(100,0),
             y=c(0,-100),
             yend=c(0,100))

p_meancentered<-test_rotation %>% 
  mutate(PC1.maize.meancenter = PC1.maize - mean(PC1.maize,na.rm=TRUE),
         PC2.maize.meancenter = PC2.maize - mean(PC2.maize,na.rm=TRUE)) %>%
  ggplot()+
  geom_point(aes(x=PC1.maize.meancenter*centroid_21km_procrustes$scale, 
                 y=PC2.maize.meancenter*centroid_21km_procrustes$scale,color = Region), shape = 3)+
  xlab(paste0("Maize PC 1 (",signif(maize_21km_transformed.eigens[1,2]*100,3),"%)")) + 
  ylab(paste0("Maize PC 2 (",signif(maize_21km_transformed.eigens[2,2]*100,3),"%)"))+
  theme_bw()+scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  ggtitle("Mean Centered maize PC1 and PC2")

p_rotated<-test_rotation %>% 
  ggplot()+
  geom_point(aes(x=rotated_PC1.maize, y=rotated_PC2.maize,color = Region))+
  xlab(paste0("Maize PC 1 (",signif(maize_21km_transformed.eigens[1,2]*100,3),"%)")) + 
  ylab(paste0("Maize PC 2 (",signif(maize_21km_transformed.eigens[2,2]*100,3),"%)"))+
  theme_bw()+scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  ggtitle("Rotated maize PC1 and PC2")

p_meancentered+geom_segment(data=axes,
                            aes(x=x,xend = xend,y=y,yend=yend))

as.matrix(centroid_21km_procrustes$rotation)%*%as.matrix(axes)
p_rotated+geom_segment(data=tibble(x=c(-75.98275,65.01248),
                                   xend=c(75.98275,-65.01248),
                                   y=c(-65.01248,-75.98275),
                                   yend=c(65.01248,75.98275)),
                       aes(x=x,xend = xend,y=y,yend=yend))

rotated_axes<-as.matrix(centroid_21km_procrustes$rotation)%*%as.matrix(axes) %>%
  as_tibble() %>% 
  mutate(maizePC = c("PC1","PC2"),
         slope = (y-yend)/(x-xend))

p_rotated+geom_abline(data = rotated_axes, 
                    aes(slope = slope, intercept = 0,linetype=maizePC))

axes<-mutate(axes, maizePC = c("PC1","PC2"),
             slope = (y-yend)/(x-xend))
#won't work because PC2 slope = infinite
p_meancentered+#geom_abline(data = axes, aes(slope=slope,intercept = 0))
    geom_vline(aes(xintercept =0),linetype=2)+geom_hline(aes(yintercept =0),linetype = 1)

#now let's try adding on humans
ggplot(test_rotation)+
  geom_point(aes(x=PC1.human, y=PC2.human, color = Region), shape=3, alpha =0.8)+
  xlab("Human PC 1")+ylab("Human PC 2")+
  geom_point(aes(x=rotated_PC1.maize, y=rotated_PC2.maize,color = Region), alpha =0.6)+
  geom_abline(data = rotated_axes, 
               aes(slope = slope, intercept = 0, linetype = maizePC))+
  xlab("Human PC 1 ") + 
  ylab("Human PC 2")+
  theme_bw()+scale_color_manual(
    values = c("#fca207","#cb4d8e","#268189","#2d1a77",
               "#b100ea","#386651","#a7c957"),
    name = "Region", labels = c("Mexico, Center","Mexico, Gulf", 
                                "Mexico, Mayan","Mexico, \nNorth of Mesoamerica","Mexico, North",
                                "Mexico, Oaxaca","Mexico, West"))+
  #scale_y_continuous(limits = c(-175,705))+
  #scale_x_continuous(limits = c(-175,705))+
  coord_fixed()

min(test_rotation$PC2.human)

## TO INCLUDE IN THE R CODE OF FINAL SCRIPT ##
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
