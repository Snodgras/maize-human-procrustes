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
