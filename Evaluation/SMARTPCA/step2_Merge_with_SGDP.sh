#!/bin/bash
#SBATCH --job-name=SGDP_Merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00 
#SBATCH --output=SGDP_Merge_%a.out
#SBATCH --error=SGDP_Merge_%a.err

#List containing all Imputed and Validation files
LIST=merging_files

#Merging the SGDP reference panel with Imputed and Validation samples
plink --bfile SGDP --merge-list $LIST --allow-no-sex --make-bed --out SGDP_MERGED

