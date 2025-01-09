#!/bin/bash 
#SBATCH -J runPCA
#SBATCH -p bigmemh
#SBATCH -o /group/jrigrp11/snodgras_maizePopulations/runPCA.out
#SBATCH -e /group/jrigrp11/snodgras_maizePopulations/runPCA.err
#SBATCH -t 4:00:00
#SBATCH --ntasks=16

module load R

Rscript --vanilla /group/jrigrp11/snodgras_maizePopulations/maizePCA.R /group/jrigrp11/snodgras_maizePopulations/filtered_centeredMaizeGBS.vcf /group/jrigrp11/snodgras_maizePopulations/filtered_centeredMaizeGBS.pca /group/jrigrp11/snodgras_maizePopulations/filtered_centeredMaizeGBS.eigenval

