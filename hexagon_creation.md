# Finding Galls-Peter equal area cylindrical projections of continents


# Making paired points between maize and humans

`slurm_anchored_procrustes.array.sh`

```
#!/bin/bash -l
#SBATCH -D /group/jrigrp11/snodgras_maizePopulations/anchored_procrustes
#SBATCH -J paired-procrustes
#SBATCH -p low2
#SBATCH -o /group/jrigrp11/snodgras_maizePopulations/pairedprocrustes_out-%A.%a.txt
#SBATCH -e /group/jrigrp11/snodgras_maizePopulations/pairedprocrustes_error-%A.%a.txt
#SBATCH -t 1:00:00
#SBATCH -c 1
#SBATCH --array=1-99

module load R
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/PairPointsByDistance.R \
	/group/jrigrp11/snodgras_maizePopulations/human_PCs/indigenous-american-genomes-popinfo.csv \
	/group/jrigrp11/snodgras_maizePopulations/human_PCs/all-indigenous.pca.csv \
	/group/jrigrp11/snodgras_maizePopulations/clean_seedspassport.csv \
	/group/jrigrp11/snodgras_maizePopulations/filtered_centeredMaizeGBS.pca \
	/group/jrigrp11/snodgras_maizePopulations/anchored_procrustes/ \
	$SLURM_ARRAY_TASK_ID
```

# Making hexagon map

QGIS
1. Open world-administrative-boundaries shapefile
2. Change projection of project to ESRI:54034 (Equal Area projection)
3. Vector/Research Tools/Create Grid
    - Grid Type: hexagon
    - Draw the extent
    - Horizontal and vertical spacing: 250 km
4. Vector/Research Tools/Select by Location to get fit grid to landmass
5. Vector/Geometry Tools/Centroids on grid 
	- Make sure that you click the box for only selected features
6. Vector/Geometry Tools/Add Geometry Attributes to centroids
	- Make sure to choose CRS WGS84
7. Vector/Data Management Tools/Join Attributes By Location to add x/y coords to grid cells
****8. Calculate neighbors for each cell (https://gis.stackexchange.com/questions/318187/generate-list-of-neighboring-polygons-qgis) and output as TSV
	- select a."Grid_500km", group_concat(b."Grid_500km") as neighbours
	- from Joined_layer a , Joined_layer b
	- where touches(a.geometry, b.geometry)
	- group by a."Grid_500km"
9. Load sample locations
10. Vector/Data Management Tools/Join Attributes By Location to add hex to sample locations

_Note to self_:
I'm restructuring the order of things above:
4. add in samples
5. select hexes that intersect samples (Vector/Geometry Tools/Select by location)
6. create a new layer with those selected hexes 
	- Right click layer name > export > save selected features
	- save as a shape file (NOT CSV)
7. Vector/Geometry Tools/ Centroids on the selected hexes
8. Vector/Geometry Tools/Add Geometry Attributes to centroids
	- make sure to choose CRS WGS84 
9. Vector/Data Management Tools/Join Attributes By Location to add x/y coords to grid cells
	- top: 500km hexes
	- bottom: added geom info-centroids
	- this makes it a hexagon rather than a point
10. Vector/Data Management Tools/Join Attributes By Location to add hex to sample locations
	- top: sample location data
	- bottom: shapes from step 9
	- `id_2`=hex id and the xcoord and ycoord are the centroid long and lat!
11. Save data with the hex sample locations as a csv
	- change the column names manually

With a 250km hexagon grid-- there are 201 hexes with maize samples and 53 with human samples
With a 500km hexagon grid-- there are 96 hexes with maize samples and 31 with human samples
With a 1000km hexagon grid-- there's 40 hexes with maize samples and 23 with human samples
_Make sure you're using the clean SeeDs dataset_

## Separate the VCF by samples within each of the unique hexes

In R
```
library(tidyverse)
maize_with_hexes<-read_csv("Desktop/RI_Lab/maize_samples/maize-seeds-passport-with-500km-hexes.csv")
#check number of unique hex_ids
maize_with_hexes$hex_id %>% unique() %>% length()
#should be 96
for(i in unique(maize_with_hexes$hex_id)){
	filter(maize_with_hexes, hex_id == i) %>% 
	select(c(Sample_ID_of_DNA_from_single_plants_used_in_GWAS,locations_latitude, locations_longitude, locations_elevation,countries_country_name, hex_id, centroid_longitude, centroid_latitude)) %>% 
	write_tsv(file = paste0("Desktop/RI_Lab/maize_samples/maizeSampleID_hex_",i,".tsv"))
	}
```
Then send files to farm
Join with the sample key to convert maize ID to VCF sample IDs
```
grep -v "^##" ../fullMaizeGBS.vcf | head -n 1 | tr '\t' '\n' > ../vcfIDs_fullMaizeGBS.txt
for f in maizeSampleID_hex*tsv ; do 
	cut -f 1 ${f} > temp.txt
	grep -w -f temp.txt ../vcfIDs_fullMaizeGBS.txt >> vcfID_${f}
	rm temp.txt
done

ml bcftools
for s in vcf*tsv ; do 
	bcftools view -S ${s} --max-alleles 2 --exclude-types indels ../fullMaizeGBS.vcf.gz > ${s#vcfID_maizeSampleID_}.vcf
done
```
The VCF has AC (allele count in genotypes (if allele is 1)) and AN (total number of alleles called)
It's unclear if there's missing data
Let's check to see if all the VCF's have the same line count (which should mean they have the same number of SNPs)
(Which means that missing data shouldn't have been filtered out)

```
bcftools view --max-alleles 2 --exclude-types indels ../fullMaizeGBS.vcf.gz | wc -l 
#should be the number of SNPs we see for each set of hex samples if there's no issues with missing data
717892
wc -l hex*vcf
```

_Assuming that it all checks out ok_

```
#print the vcf fields for SNP ID ($3) and INFO ($8)
#remove the letters and "=" from INFO field
#change ; to tab (makes allele count field 2 and total alleles field 3)
#print the SNP ID field ($1), allele counts ($2), total alleles field ($3), and allele frequency ($2/$3)
grep -v "^#" ${hexvcf} |awk -v OFS='\t' '{print $3, $8}' - | sed 's/A[C-N]=//g' | tr "\;" "\t" | awk -v OFS='\t' '{print $1, $2, $3, $2/$3}' - > ${hexfile%.tsv.vcf}.alleleFreqAndCnts.tsv
```
Then concatenate all into a single file that can be loaded into R and pivoted
```
echo SNP_ID ALLELE_COUNT ALLELE_TOTAL ALLELE_FREQUENCY HEX_ID > hex_all.alleleFreqAndCnts.tsv
sed -i 's/ /\t/g' hex_all.alleleFreqAndCnts.tsv

for h in hex_[0-9]*tsv ; do n=$(echo $h | cut -f 1 -d "." | cut -f 2 -d "_") ; sed -i "s/$/\t$n/" $h ; done

for h in hex_[0-9]*tsv ; do cat $h >> hex_all.alleleFreqAndCnts.tsv ; done

```
In R

```pivot_data_wide.R
#!/usr/bin/env Rscript
library(tidyverse)
data<-read_tsv("/group/jrigrp11/snodgras_maizePopulations/SpaceMix_files/hex_all.alleleFreqAndCnts.tsv")
pivot_wider(data, names_from = SNP_ID, values_from = c(ALLELE_COUNT, ALLELE_TOTAL, ALLELE_FREQUENCY)) %>%
	write_tsv("/group/jrigrp11/snodgras_maizePopulations/SpaceMix_files/INPUT_hex_all.alleleFreqAndCnts.tsv")
```
# make hex key file
```
cut -f 42,56,57 -d "," maize-seeds-passport-with-500km-hexes.csv > temp.csv
head -n 1 temp.csv > hex-500km-key.csv 
tail -n +2 temp.csv | sort | uniq| grep -v Collection >> hex-500km-key.csv 
```

# New pairing method between maize and human samples
Just wanting to visualize the 20km radius around human samples
And how those do or don't intersect with maize sample positions

1. Be in equal area projection
2. project human samples to equal area projection
3. Vector > Geoprocessing Tools > Buffer
4. select maize samples which intersect the buffer layer
5. save as a separate layer
6. select buffered regions that intersect with the maize layer from step 5

# To collapse multiple human samples and multiple maize samples
1. for each set of (overlapping) buffer features:
	- select
	- vector > geoprocessing tools > dissolve (make sure you check the box for selected features only)
2. merge all the dissolved layers together [vector > data management tools > merge layers]
3. edit the fully merged layer's attribute table to have an area id

# join attributes tables so each maize is assigned to a human area
1. vector > data management tools > join attributes by location
	- merge layer contains maize intersects with buffer layer
	
Do the same thing but with the buffer layer so each human is treated as an individual

# Find the samples that are unique to the 17K maize set that aren't in the GBS clean maize passport
Vector > geoprocessing tools > difference
- input layer = 17K layer, overlay = clean seeds passport

