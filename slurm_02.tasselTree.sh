#!/bin/bash 
#SBATCH -J tasselTree
#SBATCH -p bigmemm
#SBATCH -o /group/jrigrp11/snodgras_maizePopulations/tasselTree.out
#SBATCH -e /group/jrigrp11/snodgras_maizePopulations/tasselTree.err
#SBATCH -t 1:00:00
#SBATCH --ntasks=8
ml tassel
ml jdk
#run_pipeline.pl -Xmx64g -fork1 -importGuess maizeGBS.subset.vcf -tree Neighbor -export maizeGBS.subset.tree -runfork1 
run_pipeline.pl -Xmx64g -fork1 -importGuess maizeGBS.subset.method2.vcf -tree Neighbor -export maizeGBS.subset.tree -runfork1 

