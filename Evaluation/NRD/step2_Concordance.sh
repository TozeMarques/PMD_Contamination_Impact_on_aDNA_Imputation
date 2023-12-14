#!/bin/bash
#SBATCH --job-name=Concordance
#SBATCH --array 1-N
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00 
#SBATCH --output=Concordance_%a.out
#SBATCH --error=Concordance_%a.err

#Path to GLIMPSE_concordance_static
BIN=/glimpse2_alpha0.7b3_stable/GLIMPSE_concordance_static

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Allele frequency label
AF=AF

#Input files, lists created in step1
LST=lists/${ID}/${ID}_concordance.lst

#Output prefix
OUT=${ID}_

#Imputation performance evaluation using GLIMPSE_concordance
$BIN --input $LST --min-val-dp 8 --output $OUT --min-val-gl 0.9999 --bins 0.000 0.001 0.010 0.020 0.050 0.100 0.200 0.300 0.400 0.500 --af-tag $AF

