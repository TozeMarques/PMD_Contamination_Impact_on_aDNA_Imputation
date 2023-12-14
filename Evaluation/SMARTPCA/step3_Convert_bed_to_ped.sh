#!/bin/bash
#SBATCH --job-name=Bed2Ped
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00 
#SBATCH --output=Bed2Ped_%a.out
#SBATCH --error=Bed2Ped_%a.err

#Bed file containing all SGDP, Imputed and Validation samples
BED=SGDP_MERGED

#Output in ped format
OUT=SGDP_MERGED

#Converting Merged bed file into ped
$PLINK --bfile $BED --recode --tab --out $OUT

