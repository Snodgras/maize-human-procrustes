# Dataset information

v54.1.p1_HO_public.anno: 1240K+HO comes from AADR: `https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data`
It contains ancient and present-day individuals with various coverages plus present day individuals types on the Human Origins array with 597,573 sites

There are 38 Argentinan, 29 Bahamas, 192 Barbados, 24 Belize, 28 Bolivia, 72 Brazil,
(22 Canada), (5 Canary Islands), (Channel Islands 13), 25 Chile, 202 Colombia,
58 Cuba, 123 Dominican Republic,2 Haiti, 134 Mexico, 16 Panama, 265 Peru,
234 Puerto Rico, (St. Lucia? 12), (820 USA), 9 Venezuela, 

## Mexicana admixture removal from samples
use the admixture proportions already calculated
` /group/jrigrp10/mbmenon/aGWAS/Maize/MexicanaAnc_STR/admixProp.txt`
derived from `https://figshare.com/articles/dataset/mexicana_admixture_proportions_in_maize_traditional_varieties/21718076/1`

mexicana genotypes (to get allele frequency) comes from `/home/mbmenon/aGWAS/Maize/GenotypesDan`

0. Calculate mexicana allele frequencies

Filter `fullMaizeGBS.vcf` for the SNPs in the mexicana+other accessions vcf from mbmenon

```
#slurm_01.filterVCFbySNPs.sh
#uses snpID  from /home/mbmenon/aGWAS/Maize/GenotypesDan

ml bcftools

bcftools view --include ID==@snpID /group/jrigrp11/snodgras_maizePopulations/fullMaizeGBS.headerless.vcf > snpFiltered_MaizeGBS.vcf
```

Filter `AllChr_ZeaGBS.vcf` for only mexicana samples
Before this step, remove the "PL" info line from the vcf header or else it'll cause issues

```
sed -i '7d' AllChr_ZeaGBS.vcf #this removes the 7th line, only do this one time
```

Additionally, contig information needs to be added to the vcf header

```
#is this v3 or v4? v4 has longer pseudomolecules than v3
wget ftp://ftp.ensemblgenomes.org/pub/release-25/plants/fasta/zea_mays/dna/Zea_mays.AGPv3.25.dna.genome.fa.gz
ml samtools
zcat Zea_mays.AGPv3.25.dna.genome.fa.gz | samtools faidx - > Zea_mays.AGPv3.25.dna.genome.fa.gz.fai
tail Zea_mays.AGPv3.25.dna.genome.fa.gz.fai #this gives you chromosome lengths for the v3 assembly that was indexed in previous command
tail AllChr_ZeaGBS.vcf | cut -f 1-6 #check against at least chr10 coordinates from vcf

## looks like the vcf has coordinates beyond the length of the v3 chromosomes, so it must be v4!

#grab the v4 fasta
wget ftp://ftp.ensemblgenomes.org/pub/release-36/plants/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz
#then index to get the lengths needed for vcf header
zcat Zea_mays.AGPv4.dna.toplevel.fa.gz | samtools faidx - 
mv -- -.fai Zea_mays.AGPv4.dna.toplevel.fa.gz.fai

#add in the following lines to the vcf

##contig=<ID=1, length=307041717>
##contig=<ID=2, length=244442276>
##contig=<ID=3, length=235667834>
##contig=<ID=4, length=246994605>
##contig=<ID=5, length=223902240>
##contig=<ID=6, length=174033170>
##contig=<ID=7, length=182381542>
##contig=<ID=8, length=181122637>
##contig=<ID=9, length=159769782>
##contig=<ID=10, length=150982314>

```

Then move on to the filling of tags and subsetting

```

#slurm_02.filterVCFbySample.sh
#uses AllSamples_grouping.txt and AllChr_ZeaGBS.vcf from /home/mbmenon/aGWAS/Maize/GenotypesDan

ml bcftools

grep "Mexicana" AllSamples_grouping.txt | cut -f 1 > mexicana_sample_ids.txt
bcftools view -S mexicana_sample_ids.txt AllChr_ZeaGBS.vcf | bcftools +fill-tags - > AllChr_mexicanaGBS.vcf

```

The allele frequencies for each snp using all the samples in the vcf will be calculated from the `+fill-tags` option
But I'll need to be able to pull out those allele frequencies and match them to the existing data structure in the maize vcf

1. Match SNP with allele frequency in Mexicana

- add a column to the vcf with the mexicana frequency or add it as an info field 

```
ml tabix 

bgzip AllChr_mexicanaGBS.vcf 
tabix -p vcf AllChr_mexicanaGBS.vcf.gz 

bgzip AllChr_ZeaGBS.vcf
tabix -p vcf AllChr_ZeaGBS.vcf.gz

ml bcftools
bcftools annotate -a AllChr_mexicanaGBS.vcf.gz -c INFO/AF AllChr_ZeaGBS.vcf.gz > ZeaGBS_withMexAF.vcf

```

2. Match samples from vcf with samples from admixture estimate file

```
#subset the admixture data to just samples of interest
grep "SEEDGWAS" admixProp.txt > admixProp.GBSsamples.txt


ml bcftools

#subset to samples with admixture proportions for SEEDGWAS samples
cut -f 1 admixProp.GBSsamples.txt > GBSsample.id.txt
bcftools view -S GBSsample.id.txt ZeaGBS_withMexAF.vcf > SubsetZeaGBS_withMexAF.vcf

#minor allele frequency filter (removes 80021 snps [440095-360074])
bcftools view -q 0.01:minor SubsetZeaGBS_withMexAF.vcf > MAF0.01_SubsetZeaGBS_withMexAF.vcf

#to create subset with samples that have lat/long data
bcftools view -S /group/jrigrp11/snodgras_maizePopulations/clean_seedspassport.sampleIDsOnly.txt MAF0.01_SubsetZeaGBS_withMexAF.vcf > subset_with_LatLong/LatLong_MAF0.01_SubsetZeaGBS_withMexAF.vcf

#to subset only samples from Mexico
mkdir Mex_subset
#head -n 1 clean_seedspassport.csv > Mex_subset/Mex_seedspassport.csv
#grep "MEXICO" clean_seedspassport.csv >> Mex_subset/Mex_seedspassport.csv 
cd Mex_subset
cut -f 27 -d "," Mex_seedspassport.csv > Mex_SampleIDs.txt
grep -w -f Mex_SampleIDs.txt /group/jrigrp10/maize-linguistics/selected_genotypeIDs.csv > MexSampleKey.csv
cut -f 1 -d "," MexSampleKey.csv | sed 's/\"//g' |  awk -v OFS='\t' '{if($1 ~ /SEEDGWAS528[7-9]|SEEDGWAS529[0-9]|SEEDGWAS53[0-7][0-9]/) gsub(/:/, ".D17PEACXX.3.",$1); else gsub(/:/,".MRG.4.",$1 ); print $1}' -  | bcftools view -S - ../MAF0.01_SubsetZeaGBS_withMexAF.vcf > Mex_MAF0.01_SubsetZeaGBS_withMexAF.vcf

#create version without header so it can read into R
grep -v ^## MAF0.01_SubsetZeaGBS_withMexAF.vcf > headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf

grep -v ^## subset_with_LatLong/LatLong_MAF0.01_SubsetZeaGBS_withMexAF.vcf > subset_with_LatLong/headerless_LatLong_MAF0.01_SubsetZeaGBS_withMexAF.vcf
grep -v ^## Mex_MAF0.01_SubsetZeaGBS_withMexAF.vcf > headerless_Mex_MAF0.01_SubsetZeaGBS_withMexAF.vcf

```

4. Plink LD prune
```slurm_04.plinkpca.sh 
VCF=/group/jrigrp11/snodgras_maizePopulations/MAF0.01_SubsetZeaGBS_withMexAF.vcf
VCF=/group/jrigrp11/snodgras_maizePopulations/admix_removal/Mex_subset/Mex_MAF0.01_SubsetZeaGBS_withMexAF.vcf

ml conda/plink/1.90b6.21
#identify prunning sites due to LD threshold 
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out maize

#prune (I don't think we need this step)
plink --vcf $VCF --double-id --allow-extra-chr --set-missing-var-ids @:# --extract maize.prune.in --make-bed --pca --out maize
```

3. Run transformation on genotypes
G' = (G - 2 x mexicanaAlleleFrequency x alpha) / (1 - alpha)

For a snp, find the mexicana allele frequency for that snp
For each individual with a genotype at that SNP, 
- find the admix proportion (alpha) for that individual
- apply the transformation to the genotype
- write the new genotype

Pilot test:
create random allele frequencies for mexicana and random genotypes for 10 individuals and random alphas for those individuals

Subset test: 
head -n 5000 ` MAF0.01_SubsetZeaGBS_withMexAF.vcf` | grep -v ^## as test.vcf

Final draft: will output transformed genotypes to `transformed_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv`


Then do the LD filtering on transformed genotypes (360046 --> 134813 (225233 removed))
```
head -n 1 transformed_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv > LDprune_transformed_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
grep -w -f maize.prune.in transformed_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv >> LDprune_transformed_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/0\/0/0/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/0\/1/1/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/1\/0/1/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/1\/1/2/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/\.\/\./NA/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
```

This didn't work but this was done in the directory `subset_with_LatLong` along with the full set untransformed

5. Run PCA
Make sure that the centering is flipped to T
Use the script `maizePCA.R`
Invoked in slurm as: 
`Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/maizePCA.R /group/jrigrp11/snodgras_maizePopulations/admix_removal/LDprune_transformed_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv /group/jrigrp11/snodgras_maizePopulations/admix_removal/LDprune_transformed_MAF0.01_SubsetZeaGBS.pca /group/jrigrp11/snodgras_maizePopulations/admix_removal/LDprune_transformed_MAF0.01_SubsetZeaGBS.eigenvec`

### Make test set for weird PCA results

20000 Random SNPs
```
library(tidyverse)
LD_prune<-read_tsv("maize.prune.in", col_names="ID")
full<-read_tsv("headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf") %>% filter(ID %in% LD_prune$ID)
full[sample(nrow(full), 20000),] %>% write_tsv("test/20000RandomSNPs_withHeader.tsv")
```

Untransformed test
```
head -n 1 headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv > non_transformed_pca/LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
grep -w -f maize.prune.in headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv >> non_transformed_pca/LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/0\/0/0/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/0\/1/1/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/1\/0/1/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/1\/1/2/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv
sed -i 's/\.\/\./NA/g' LDprune_headerless_MAF0.01_SubsetZeaGBS_withMexAF.tsv

```

Now run on the maize samples within 21 km of a human sample
This makes the VCF of just those maize samples
```
cd /group/jrigrp11/snodgras_maizePopulations/21km_maize/

#This should subset the samples to just the ones you want and only write biallelic SNPs from the VCF
ml bcftools 
cut -f 1  2024-attempt/21km_maize_samples_ids.tsv| tail -n +2 | sed 's/\"//g'| awk -v OFS='\t' '{if($1 ~ /SEEDGWAS528[7-9]|SEEDGWAS529[0-9]|SEEDGWAS53[0-7][0-9]/) gsub(/:/, ".D17PEACXX.3.",$1); else gsub(/:/,".MRG.4.",$1 ); print $1}' -   | bcftools view -S - ../admix_removal/MAF0.01_SubsetZeaGBS_withMexAF.vcf > MAF0.01_21km_MaizeGBS.vcf
```

Then LD prune:

```
ml conda/plink/1.90b6.21
#identify prunning sites due to LD threshold 
plink --vcf MAF0.01_21km_MaizeGBS.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out maize

grep "^#CHROM" MAF0.01_21km_MaizeGBS.vcf > LDprune_MAF0.01_21km_MaizeGBS.tsv
grep -w -f maize.prune.in MAF0.01_21km_MaizeGBS.vcf  >> LDprune_MAF0.01_21km_MaizeGBS.tsv

#then change genotype scores to be 0, 1, 2 (numeric)
sed -i 's/0\/0/0/g' LDprune_MAF0.01_21km_MaizeGBS.tsv
sed -i 's/0\/1/1/g' LDprune_MAF0.01_21km_MaizeGBS.tsv
sed -i 's/1\/0/1/g' LDprune_MAF0.01_21km_MaizeGBS.tsv
sed -i 's/1\/1/2/g' LDprune_MAF0.01_21km_MaizeGBS.tsv
sed -i 's/\.\/\./NA/g' LDprune_MAF0.01_21km_MaizeGBS.tsv
```

then run `21km_admixture_removal.R`

To add in sample ID column for genotypes:

```
head -n 1 LDprune_MAF0.01_21km_MaizeGBS.tsv | cut -f 10- | tr '\t' '\n' > 21km_transformed_LDprune_MAF0.01.genotypes.sampleIDs.txt
```

Edit the PVE R script to read in the above file and add as a column to the genotype file

### PVE R Script
To make it run faster, divide the input genotype file into multiple files to run on low

```
genotype_in=$1

i=0
k=$(wc -l $genotype_in | cut -f 1 -d " " )

while [ $i -lt $k  ] ; do
        head -n 1 $genotype_in > divided_raw_genotypes/${genotype_in%.tsv}.${i}.txt
        j=`expr $i + 1`
        tail -n +${j} ${genotype_in} | head -n 1000 >> divided_raw_genotypes/${genotype_in%.tsv}.${i}.txt
        i=`expr $i + 1000`
done
```

```
#SBATCH -J LM_maizePVEbyHumanPC_1
#SBATCH -p high2
#SBATCH -o /group/jrigrp11/snodgras_maizePopulations/21km_maize/runMaizePVEbyHumanPCs.1.raw.out
#SBATCH -e /group/jrigrp11/snodgras_maizePopulations/21km_maize/runMaizePVEbyHumanPCs.1.raw.err
#SBATCH -t 12:00:00
#SBATCH --ntasks=16
#SBATCH --mem 128G

ml R

Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R raw /group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv PC_1 /group/jrigrp11/snodgras_maizePopulations/21km_maize/PC1.raw_LDprune_MAF0.01.lm_results.tsv

```
```
w.write("#!/bin/bash\n")
        w.write("#SBATCH -p low\n")
        w.write("#SBATCH -n 16\n")
        w.write("#SBATCH -t 01:00:00\n") #TIME REQUEST
        w.write("#SBATCH --mem 128G\n") #MEMORY REQUEST
        w.write("#SBATCH -J "+jobname+"_"+str(filecount)+"\n")
        w.write("#SBATCH -o "+jobname+"_"+str(filecount)+".o%j\n")
        w.write("#SBATCH -e "+jobname+"_"+str(filecount)+".e%j\n")
        
for i in divided_raw_genotypes/* ; do 
	echo ml R \; Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R raw /group/jrigrp11/snodgras_maizePopulations/21km_maize/${i} PC_1 /group/jrigrp11/snodgras_maizePopulations/21km_maize/divided_raw_lm_output/${i#divided_raw_genotypes/}.lm_results.tsv >> raw_lm_PC1.cmds ; 
done
for i in 0 ; do 
	while read -r line; do 
		bash makeSLURM.sh slurm_${i}_raw_PC1.cmd raw_PC1_${i} ; 
		echo ${line} >> slurm_${i}_raw_PC1.cmd ;  
		i=`expr $i + 1` ;  
	done < raw_lm_PC1.cmds ; 
done

bash create_slurms.sh PC_1+PC_2 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_2.cmds
bash create_slurms.sh PC_1+PC_2+PC_3 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_3.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_4.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_5.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_6.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_7.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_8.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_9.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_PC1_10.cmds

for i in {3..10} ; do for j in {0..138} ; do sbatch slurm_${j}_raw_lm_PC1_${i}.cmds ; done ; done

mkdir divided_transformed_genotypes
mkdir divided_transformed_lm_output

for i in 1 ; do 
	n=$(awk '{print NF}' 21km_transformed_LDprune_MAF0.01.genotypes.tsv | head -n 1) ; 
	while [ $i -lt $n ] ; do 
		j=`expr $i + 999` ; 
		cut -f ${i}-${j} 21km_transformed_LDprune_MAF0.01.genotypes.tsv > divided_transformed_genotypes/21km_transformed_LDprune_MAF0.01.genotypes.${i}.txt ; 
		i=`expr $i + 1000` ; 
		echo $i
	done ; 
done

bash create_slurms.sh PC_1 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1.cmds
bash create_slurms.sh PC_1+PC_2 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_2.cmds
bash create_slurms.sh PC_1+PC_2+PC_3 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_3.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_4.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_5.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_6.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_7.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_8.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_9.cmds
bash create_slurms.sh PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_PC1_10.cmds

for i in {2..10} ; do for j in {0..138} ; do sbatch slurm_${j}_transformed_lm_PC1_${i}.cmds ; done ; done


```

*FOR SOME REASON SNP FILE 1000 of RAW IS NOT RUNNING, SAYS THAT THERE'S AN NA IN THE DATAFRAME NAMES*
Last line was empty (full of NAs)

```
human_21km_pairs<-read_csv("Individual-human-20km-maize-samples-joined.csv")
```

*WHY IS THE PC1-10 SO DIFFERENT FROM THE FIRST LM RUN THROUGH?* They aren't actually that different...?

Adding in latitude, longitude

_edit MaizePVEbyHumansPCs.parallelable.R so that transposed genotype object also includes lat/long/elevation columns_
_edit create slurms script to avoid weird output file names_ 

slurm scripts named: `slurm_[SNPdivisionNumber0-138]_[raw/transformed]_lm_latlong.cmds`

low jobs currently hung on: `transformed_locations_latitude_longitude_2`

```
#transformed
#_TESTING ON THIS COMMAND, division 0#
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlong.cmds

bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_2.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_3.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_4.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_5.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_6.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_7.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_8.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_9.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10 transformed divided_transformed_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/transformed_lm_latlongPC1_10.cmds

#raw
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlong.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_2.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_3.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_4.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_5.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_6.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_7.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_8.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_9.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongPC1_10.cmds

```
Including elevation

```
#raw
bash create_slurms.sh locations_elevation.maize raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_elev.cmds

bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelev.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_2.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2+PC_3 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_3.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2+PC_3+PC_4 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_4.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2+PC_3+PC_4+PC_5 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_5.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_6.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_7.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_8.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_9.cmds
bash create_slurms.sh locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10 raw divided_raw_genotypes/ /group/jrigrp11/snodgras_maizePopulations/21km_maize/slurms_and_logs/raw_lm_latlongelevPC1_10.cmds

```
To make concatenated result files:
```
head -n 1 divided_transformed_lm_output/21km_transformed_LDprune_MAF0.01.genotypes.1.PC_1.lm_results.tsv > 21km_transformed_LDprune_MAF0.01.genotypes.all_models.lm_results.tsv
for i in divided_transformed_lm_output/*.tsv ; do tail -n +2 $i >> 21km_transformed_LDprune_MAF0.01.genotypes.all_models.lm_results.tsv ; done

head -n 1 divided_raw_lm_output/LDprune_MAF0.01_21km_MaizeGBS.1.PC_1.lm_results.tsv > 21km_raw_LDprune_MAF0.01.genotypes.all_models.lm_results.tsv
for i in divided_raw_lm_output/*.tsv ; do tail -n +2 $i >> 21km_raw_LDprune_MAF0.01.genotypes.all_models.lm_results.tsv ; done

```
Using Forrest's code to parallelize (`foreach()`)
```
#!/bin/bash
#SBATCH -p high2
#SBATCH -n 16
#SBATCH -t 03:00:00
#SBATCH --mem 128G
#SBATCH -J LLEA_raw_lm
#SBATCH -o slurms_and_logs/LLEA_raw_lm.out
#SBATCH -e slurms_and_logs/LLEA_raw_lm.err
ml R

Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmix_lm_results.tsv
```
```
#!/bin/bash
#SBATCH -p high2
#SBATCH -n 16
#SBATCH -t 03:00:00
#SBATCH --mem 128G
#SBATCH -J LLEA_raw_lm
#SBATCH -o slurms_and_logs/LLEA_PCs_raw_lm.out
#SBATCH -e slurms_and_logs/LLEA_PCs_raw_lm.err
ml R

Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9+PC_10 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_10_lm_results.tsv
	
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8+PC_9 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_9_lm_results.tsv
	
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7+PC_8 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_8_lm_results.tsv

Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6+PC_7 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_7_lm_results.tsv

Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2+PC_3+PC_4+PC_5+PC_6 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_6_lm_results.tsv

Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2+PC_3+PC_4+PC_5 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_5_lm_results.tsv
	
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2+PC_3+PC_4 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_4_lm_results.tsv

Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2+PC_3 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_3_lm_results.tsv
	
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1+PC_2 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_2_lm_results.tsv
	
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/21km_maize/MaizePVEbyHumansPCs.parallelable.R \
	raw \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/LDprune_MAF0.01_21km_MaizeGBS.tsv \
	locations_latitude.maize+locations_longitude.maize+locations_elevation.maize+Admix+PC_1 \
	/group/jrigrp11/snodgras_maizePopulations/21km_maize/21km_raw_LDprune_MAF0.01_LatLongElevAdmixPC1_lm_results.tsv
```

Top Snps to potentially look at:
```
Chromosome        bp    snp Difference LLEA10
   <fct>          <dbl>  <dbl>      <dbl>  <dbl>
 1 7          127679528  99430      0.727  0.863
 2 1           35374876   4724      0.655  0.806
 3 1          306080024  21387      0.612  0.685
 4 3          187345731  48634      0.612  0.685
 5 5          206102482  79283      0.612  0.685
 6 10         135307922 134812      0.612  0.685
 7 1            6526734   1156      0.478  0.545
 8 8           62427842 109884      0.462  0.533
 9 2            4233687  22525      0.401 NA    
10 1          194362111  11470      0.381 NA    
11 1          297413402  20056      0.359 NA    
12 8          164858047 114822      0.347 NA    
13 2          238199640  37610      0.337 NA    
14 4          115551780  58826      0.336 NA    
15 9          148606533 126514     NA      0.649
16 8          126885347 111916     NA      0.539
17 1          172712468  10026     NA      0.525
18 5          188595628  77338     NA      0.516
19 3           84443353  43690     NA      0.514
20 5          200859705  78648     NA      0.513

so SNP names should be: (in snps_of_interest.txt)
S7_127679528 
S1_35374876  
S1_306080024 
S3_187345731 
S5_206102482 
S10_135307922
S1_6526734   
S8_62427842  
S2_4233687   
S1_194362111 
S1_297413402 
S8_164858047 
S2_238199640 
S4_115551780 
S9_148606533 
S8_126885347 
S1_172712468 
S5_188595628 
S3_84443353  
S5_200859705

head -n 1 LDprune_MAF0.01_21km_MaizeGBS.tsv > snps_of_interest.genotypes
grep -f snps_of_interest.txt LDprune_MAF0.01_21km_MaizeGBS.tsv >> snps_of_interest.genotypes

head -n 1 /group/jrigrp11/snodgras_maizePopulations/admix_removal/headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf > snps_of_interest.allSamples.genotypes
grep -f snps_of_interest.txt /group/jrigrp11/snodgras_maizePopulations/admix_removal/headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf >> snps_of_interest.allSamples.genotypes
```

## 100km human-maize pairs
```
# 1. Generate human-maize pairs file

## create 100km_maize_human_pairs.csv using qgis
## replace the header names since qgis cuts them down

head -n 1 clean_seedspassport.csv 
head -n 1 indigenous-american-genomes-popinfo.csv

#use vi to replace header

# 2. Create MAF 0.01 filtered genotype file of those specific maize samples

tail -n +2 100km_maize_human_pairs.csv | cut -f 27 -d , | sort | uniq | grep -v "NA" > 100km_maize_ids.txt

#because the sample ids don't exactly match between files, need to change my ids to that of the vcf
grep "SEEDGWAS" /group/jrigrp11/snodgras_maizePopulations/admix_removal/MAF0.01_SubsetZeaGBS_withMexAF.vcf | cut -f 10- | tr "\t" "\n" > vcf_maize_sample_ids.txt
awk -v OFS='\t' -F "." '{print $0,$1}' vcf_maize_sample_ids.txt > maize_sample_id_key.txt

join -1 2 -2 1 maize_sample_id_key.txt 100km_maize_ids.txt > 100km_vcf_maize_ids.txt

ml bcftools

cut -f 2 -d " " /group/jrigrp11/snodgras_maizePopulations/100km_maize/100km_vcf_maize_ids.txt | bcftools view -S -\
	/group/jrigrp11/snodgras_maizePopulations/admix_removal/MAF0.01_SubsetZeaGBS_withMexAF.vcf \
	> /group/jrigrp11/snodgras_maizePopulations/100km_maize/100km_MAF0.01_SubsetZeaGBS_withMexAF.vcf

# 3. LD prune maize genotypes

ml plink
# plink/1.9-beta7.1

#identify prunning sites due to LD threshold 
plink --vcf /group/jrigrp11/snodgras_maizePopulations/100km_maize/100km_MAF0.01_SubsetZeaGBS_withMexAF.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out maize

grep "^#CHROM" 100km_MAF0.01_SubsetZeaGBS_withMexAF.vcf  > LDprune_MAF0.01_100km_MaizeGBS.tsv
grep -w -f maize.prune.in 100km_MAF0.01_SubsetZeaGBS_withMexAF.vcf >> LDprune_MAF0.01_100km_MaizeGBS.tsv

#then change genotype scores to be 0, 1, 2 (numeric)
sed -i 's/0\/0/0/g' LDprune_MAF0.01_100km_MaizeGBS.tsv
sed -i 's/0\/1/1/g' LDprune_MAF0.01_100km_MaizeGBS.tsv
sed -i 's/1\/0/1/g' LDprune_MAF0.01_100km_MaizeGBS.tsv
sed -i 's/1\/1/2/g' LDprune_MAF0.01_100km_MaizeGBS.tsv
sed -i 's/\.\/\./NA/g' LDprune_MAF0.01_100km_MaizeGBS.tsv


# 4. Run lm script

```
Top SNPs of interest:
```
S6_53598250
S9_148365668
S1_172713896
S9_13141167
S5_188595628
S1_23618793
S1_172961727
S5_104847918
S5_110897416
S3_48601830
S2_77476338
S5_200859705
S5_107545947
S6_108202924

head -n 1 /group/jrigrp11/snodgras_maizePopulations/admix_removal/headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf > 100kmsnps_of_interest.allSamples.genotypes
grep -f 100km_snps_of_interest.txt /group/jrigrp11/snodgras_maizePopulations/admix_removal/headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf >> 100kmsnps_of_interest.allSamples.genotypes

```
Top SNPs that are top for difference between LLEA and LLEA+10PCs
```
S9_141814044
S7_147837873
S2_235925301
S3_18734573
S1_100651817
S8_172039281
S8_2566123
S2_145562528
S1_235942075
S8_164811862
S9_68301736
S2_1622935
S6_160173021
S2_54984271

head -n 1 /group/jrigrp11/snodgras_maizePopulations/admix_removal/headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf > 100km_diff_snps_of_interest.allSamples.genotypes
grep -f 100km_diff_snps_of_interest.txt /group/jrigrp11/snodgras_maizePopulations/admix_removal/headerless_MAF0.01_SubsetZeaGBS_withMexAF.vcf >> 100km_diff_snps_of_interest.allSamples.genotypes

```


## QGIS creating intersections

1. make sure that all data is in equal area projection 
	Vector > Data Management Tools > Reproject layer
2. Draw buffer around human data points
	Vector > Geoprocessing Tools > Buffer
3. Intersect buffer layer with maize points
4. Merge attribute tables?