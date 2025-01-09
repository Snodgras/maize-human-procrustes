#!/bin/bash 
#SBATCH -J tasselPCA
#SBATCH -p bigmemh
#SBATCH -o /group/jrigrp11/snodgras_maizePopulations/tasselPCA.out
#SBATCH -e /group/jrigrp11/snodgras_maizePopulations/tasselPCA.err
#SBATCH -t 12:00:00
#SBATCH --ntasks=16
ml tassel
ml jdk
#run_pipeline.pl -Xmx126g -fork1 -importGuess maizeGBS.subset.vcf -PrincipalComponentsPlugin -covariance true -endPlugin -export maizeGBS.subset.pca -runfork1
#run_pipeline.pl -Xmx64g -fork1 -importGuess maizeGBS.subset.method2.vcf -PrincipalComponentsPlugin -covariance true -endPlugin -export maizeGBS.subset.pca -runfork1
run_pipeline.pl -Xmx126g -fork1 -importGuess fullMaizeGBS.vcf -PrincipalComponentsPlugin -covariance true -endPlugin -export maizeGBS.full.pca -runfork1

