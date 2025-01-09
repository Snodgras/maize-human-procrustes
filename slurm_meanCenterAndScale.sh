#!/bin/bash 
#SBATCH -J meanCenterAndScale
#SBATCH -p bigmemh
#SBATCH -o /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.out
#SBATCH -e /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.err
#SBATCH -t 4:00:00
#SBATCH --ntasks=16

module load R

#replace the genotype strings with numbers
#sed -i 's/0|0/0/g' fullMaizeGBS.headerless.vcf
#sed -i 's/0|1/1/g' fullMaizeGBS.headerless.vcf
#sed -i 's/1|0/1/g' fullMaizeGBS.headerless.vcf
#sed -i 's/1|1/2/g' fullMaizeGBS.headerless.vcf

#remove any lines that are not biallelic (have alt alleles 2, 3)
#grep -v -e "0|2" -e "0|3" -e "1|2" -e "1|3" -e "2|0" -e "2|1" -e "2|2" -e "2|3" -e "3|0" -e "3|1" -e "3|2" -e "3|3" fullMaizeGBS.headerless.vcf > biallelic.headerless.vcf


#Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R biallelic.headerless.vcf /group/jrigrp11/snodgras_maizePopulations/

#Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.1-100000.
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.2.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.100001-200000.
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.3.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.200001-300000.
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.4.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.300001-400000.
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.5.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.400001-500000.
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.6.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.500001-600000.
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.7.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.600001-700000.
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.8.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.700001-800000.
Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/meanCenterAndScale.R test.9.vcf /group/jrigrp11/snodgras_maizePopulations/ temp.800001-851597.

cat temp.1-100000.centeredMaizeGBS.vcf > centeredMaizeGBS.vcf
for i in temp.[1-9]0* ; do tail -n +2 $i >> centeredMaizeGBS.vcf ; done
rm temp* test*.vcf
