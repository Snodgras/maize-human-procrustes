#!/bin/bash 
#SBATCH -J filterVCF
#SBATCH -p med
#SBATCH -o /group/jrigrp11/snodgras_maizePopulations/filterVCF_out
#SBATCH -e /group/jrigrp11/snodgras_maizePopulations/filterVCF_error
#SBATCH -t 2:00:00

ml bcftools 
#echo Starting zipping job...
#bgzip -c fullMaizeGBS.vcf > fullMaizeGBS.vcf.gz
#echo Done zipping and starting indexing...
#tabix -p vcf fullMaizeGBS.vcf.gz
#echo Done indexing and starting filtering...
#cut -f 1 -d "," maizeSampleKey.csv | bcftools view --samples-file - fullMaizeGBS.vcf.gz > maizeGBS.subset.vcf

cut -f 1 -d "," maizeSampleKey.method2.csv | bcftools view --samples-file - fullMaizeGBS.vcf.gz > maizeGBS.subset.method2.vcf
