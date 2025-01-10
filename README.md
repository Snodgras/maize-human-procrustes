# maize-human-procrustes

This document outlines the general steps and scripts involved in doing the various procrustes analyses.

# 1. Identifying and removing outliers

PC values that are greater than 3 standard deviations away from the mean PC value are flagged as outliers. 
Outlier samples are removed from input files and PCs are recalculated without those outlier samples. 
Individuals who are the product of recent migrations (< 5 generations) are also removed. 

# 2. Calculating PCA

## Human data
Humans are categorized into 12 regions based on ethno-linguistic groups that were self-identified during sample collection. 
PCA was calculated using `fastPCA` after LD pruning.

The different subsets tried:
- all data
- outliers and recent migrants removed
- only data from Mexico with all outliers and recent migrants removed

## Maize data
Maize groups are defined either by geopolitical boundaries (i.e., country of origin) or their geographical distance to human ethnolinguistic groups. 

The different subsets tried: 
- all data
- outliers removed
- only data from Mexico
- mexicana ancestry removed and outliers removed, restricted to Mexican samples

# 3. Removal of mexicana ancestry from maize samples

# 4. Pairing of maize and human samples

# 5. Procrustes analyses

Procrustes calculated using the R package `vegan`. 

## geography vs. population structure
X dataframe is longitude and latitude (the non-transformed positions) and Y dataframe is PC 1 and PC 2 of either human or maize PCA data (the transformed positions). Thus PC space is converted to lat/long coordinates. 
