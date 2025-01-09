#!/usr/bin/env Rscript

library(tidyverse)
library(vegan)
library(geosphere, lib.loc = "/home/snodgras/R_Packages/R4.2.3/")

args = commandArgs(trailingOnly = TRUE)

#Load in data
human_meta<-read_csv(args[1]) #human metadata with lat/long
human_pcs<-read_csv(args[2]) #human PCs

maize_meta<-read_csv(args[3]) #maize metadata with lat/long
maize_pca<-read_tsv(args[4]) #maize PCs
colnames(maize_pca)[1]<-"ind" #change first column name

tag<-args[5] #path to directory where output should be written
iteration<-args[6] #run number

print("Arguments loaded...")

####Function to find the nearest maize sample####

calculate_distance<-function(h.long, h.lat, m.long, m.lat){
  return(distGeo(c(h.long,h.lat),c(m.long,m.lat)) / 1000)
}

anchors<-tibble(human_id=NA, maize_id=NA, distance=NA) #final df
human_index<-sample(1:nrow(human_meta),nrow(human_meta),replace = F)
for(j in human_index){ #for each human sample
  distance<-c() #make distance vector
  for(i in 1:nrow(maize_meta)){ #for each maize sample
    #calculate distance from human sample
    distance<-c(distance, calculate_distance(human_meta$longitude[j], 
                                             human_meta$latitude[j],
                                             maize_meta$locations_longitude[i],
                                             maize_meta$locations_latitude[i]))
  }
  index<-order(distance) #provides all indices (to be able to get maize ID) ordered from min to max distance
  min_samples<-maize_meta$Sample_ID_of_DNA_from_single_plants_used_in_GWAS[index[1]]
  k<-1
  if(min_samples %in% anchors$maize_id){
    #only run the while loop if there's an overlap in sample IDs
    while(min_samples %in% anchors$maize_id){
      #while the maize sample id is in the anchors list
      k<-k+1 #iterate
      min_samples<-maize_meta$Sample_ID_of_DNA_from_single_plants_used_in_GWAS[index[k]]
    }
  }
  #add row to the final df
  anchors<-add_row(anchors, human_id = human_meta$sample_id[j],
                   maize_id = min_samples,
                   distance = distance[k]) #k comes from the index for the while loop
}
anchors<-na.omit(anchors)

### what if we did our best to minimize distance throughout?
find_min_distance<-function(h.long,h.lat,m.long,m.lat,m.df){
  distance<-c()
  for(m in 1:nrow(m.df)){
    distance<-c(distance,
                calculate_distance(h.long,h.lat,pull(m.df[m,m.long]),pull(m.df[m,m.lat])))
  }
  return(min(distance))
}

human_meta$distance<-NA
for(i in 1:nrow(human_meta)){
  distance<-c()
  for(m in 1:nrow(maize_meta)){
    distance<-c(distance,calculate_distance(human_meta$longitude[i],
                                            human_meta$latitude[i], 
                                            maize_meta$locations_longitude[m],
                                            maize_meta$locations_latitude[m]))
  }
  human_meta$distance[i]<-min(distance)
}

human_meta$distance %>% summary()

filter(human_meta, distance < 21) %>% nrow() #348
filter(human_meta, distance < 31) %>% nrow() #421
filter(human_meta, distance < 21 & sample_id %in% human_noOutliers$sample_name) %>% nrow() #332
filter(human_meta, distance < 31 & sample_id %in% human_noOutliers$sample_name) %>% nrow() #402

filter(human_meta, distance < 31 & sample_id %in% human_noOutliers$sample_name) %>%
  ggplot(aes(x=longitude, y=latitude, color = Region))+
  geom_point()+
  theme_bw()

trial_human_meta<-filter(human_meta, distance < 21 & sample_id %in% human_noOutliers$sample_name)

trial_anchors<-tibble(human_id=NA, maize_id=NA, distance=NA) #final df
human_index<-sample(1:nrow(trial_human_meta),nrow(trial_human_meta),replace = F)
for(j in human_index){ #for each human sample
  distance<-c() #make distance vector
  for(i in 1:nrow(maize_meta)){ #for each maize sample
    #calculate distance from human sample
    distance<-c(distance, calculate_distance(trial_human_meta$longitude[j], 
                                             trial_human_meta$latitude[j],
                                             maize_meta$locations_longitude[i],
                                             maize_meta$locations_latitude[i]))
  }
  index<-order(distance) #provides all indices (to be able to get maize ID) ordered from min to max distance
  min_samples<-maize_meta$Sample_ID_of_DNA_from_single_plants_used_in_GWAS[index[1]]
  k<-1
  if(min_samples %in% trial_anchors$maize_id){
    #only run the while loop if there's an overlap in sample IDs
    while(min_samples %in% trial_anchors$maize_id){
      #while the maize sample id is in the anchors list
      k<-k+1 #iterate
      min_samples<-maize_meta$Sample_ID_of_DNA_from_single_plants_used_in_GWAS[index[k]]
    }
  }
  #add row to the final df
  trial_anchors<-add_row(trial_anchors, human_id = trial_human_meta$sample_id[j],
                   maize_id = min_samples,
                   distance = distance[k]) #k comes from the index for the while loop
}
trial_anchors<-na.omit(trial_anchors)
summary(trial_anchors$distance)
#still high (Q1 jumps to 300 km after the min = 20.96 (BUT I KNOW THERE"S SOME THAT ARE RIGHT ON TOP OF EACH OTHER!))

inner_join(select(trial_human_meta, sample_id, latitude, longitude), trial_anchors, by=c("sample_id" = "human_id"))%>%
  ggplot(aes(x=longitude, y=latitude, color = distance))+
  geom_point()+
  scale_color_viridis_c()+
  theme_bw()

#check how many maize samples are there to choose from in that vicinity around each of the human samples
maize_meta$humans_within_20km<-NA
maize_meta$humans_within_30km<-NA
for(m in 1:nrow(maize_meta)){
  distance<-c()
  for(h in 1:nrow(trial_human_meta)){
    distance<-c(distance, calculate_distance(trial_human_meta$longitude[h], 
                                             trial_human_meta$latitude[h],
                                             maize_meta$locations_longitude[m],
                                             maize_meta$locations_latitude[m]))
  }
  maize_meta$humans_within_20km[m]<-sum(distance < 21, na.rm = T)
  maize_meta$humans_within_30km[m]<-sum(distance < 31, na.rm = T)
}

filter(maize_meta, humans_within_30km > 0) %>% nrow()

#so there's enough maize samples for each human sample? I suppose a multiple maize samples could be that close to the same human sample? 
#seems like a problem with the actual matching function then...

print("Anchors generated...")

#By randomizing the order of human samples, 
#you will change the distribution of distances at the lower end

####Pull out the PC values for the anchors####
colnames(human_pcs)[which( str_detect(colnames(human_pcs), "PC"))]<-colnames(human_pcs)[which( str_detect(colnames(human_pcs), "PC"))] %>% str_c(".human")
colnames(maize_pca)[which( str_detect(colnames(maize_pca), "PC"))]<-colnames(maize_pca)[which( str_detect(colnames(maize_pca), "PC"))] %>% str_c(".maize")

anchors<-inner_join(x=anchors, y=human_pcs, by = join_by("human_id" == "sample_name"))
anchors<-mutate(maize_pca, ID = str_split(ind, ":", simplify = T)[,1]) %>%
  inner_join(x=anchors, y=., by = join_by("maize_id" == "ID"))

####Run procrustes on the matched anchors####

anchor.protest<-
  protest(Y = select(anchors, c("PC1.maize","PC2.maize")),
           X = select(anchors, c("PC_1.human","PC_2.human"))
          )
print("Procrustes run...")
####Pull out fit metrics of protest####
#anchor.protest$ss #sum of squares
#anchor.protest$t0 #correlation
#anchor.protest$signif #significance of correlation with permutation

#anchor_output<-
tibble(SumOfSquares = anchor.protest$ss,
       Correlation = anchor.protest$t0,
       Significance = anchor.protest$signif,
       Dist.min = min(anchors$distance),
       Dist.max = max(anchors$distance),
       Dist.mean = mean(anchors$distance),
       Dist.median = median(anchors$distance),
       Dist.Q1 = quantile(anchors$distance)[2],
       Dist.Q3 = quantile(anchors$distance)[4]) %>%
  write_csv(path = paste0(tag,"anchor.",iteration,".out.csv"))
write_tsv(anchors, path = paste0(tag,iteration,".anchors.tsv"))

print("Output written...All done :)")

####Pull out transformation values####
#Not sure if we'll want this or not#