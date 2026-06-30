#!/bin/bash -l

#SBATCH -D /group/jrigrp10/maize-linguistics/scripts
#SBATCH -o /home/fli21/slurm-log/regression-linguistics-%j.txt
#SBATCH -e /home/fli21/slurm-log/regression-linguistics-%j.txt
#SBATCH -J regression-linguistics
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu 4GB
#SBATCH -n 16
#SBATCH --account=jrigrp
#SBATCH --partition=high2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=frrli@ucdavis.edu

module load R

Rscript regressionSNPs.R mayan NL