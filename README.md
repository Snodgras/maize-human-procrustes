# maize-human-procrustes

This outlines the general steps and scripts involved in doing the various procrustes analyses.
Most of this is completed in the script `plottingPCbyGeography.R` which follows the same layout. 
Exceptions noted. 

# 1. Identifying and removing outliers

PC values that are greater than 3 standard deviations away from the mean PC value are flagged as outliers. 
Outlier samples are removed from input files and PCs are recalculated without those outlier samples. 
Individuals who are the product of recent migrations (< 5 generations) are also removed. 

# 2. Calculating PCA

For more detailed notes of scripts and code usage for each step, see `Unix_Notes.md`

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

See `21km_admixture_removal.R`for how admixture was removed while results are plotted in `plottingPCbyGeography.R`

Essentially missing genotypes are mean imputed and the effects of estimated Mexicana admixture is regressed out of the genotypes. 

# 4. Pairing of maize and human samples

This was done in QGIS:

## QGIS creating intersections

1. make sure that all data is in equal area projection 
	Vector > Data Management Tools > Reproject layer
2. Draw buffer around human data points
	Vector > Geoprocessing Tools > Buffer
3. Intersect buffer layer with maize points
4. Merge attribute tables

Adding features about geographic location groupings to sample metadata came from the `MapMaking.R` script (as well as general geography figures).

# 5. Procrustes analyses

Procrustes calculated using the R package `vegan`. 

## geography vs. population structure
X dataframe is longitude and latitude (the non-transformed positions) and Y dataframe is PC 1 and PC 2 of human or PC 2 and 3 of maize PCA data (the transformed positions). Thus PC space is converted to lat/long coordinates. 

## maize vs. human PCs
X dataframe is the human centroid PCs and the Y dataframe is the maize samples' PCs
Pairings based off the QGIS step (4)

# 6. Correlation of PCs and linear regressions
Correlations are completed in `plottingPCbyGeography.R` and `PCcorrelations.R`
The linear regressions are run separately in parallelized scripts.
Where genotypes are used as the response variable, see: `MaizePVEbyHumansPCs.parallelable.R`. 
Where maize PCs are used as the response variable, see: `MaizePVEbyHumansPCs.ResponseToReview2026.R`
 
