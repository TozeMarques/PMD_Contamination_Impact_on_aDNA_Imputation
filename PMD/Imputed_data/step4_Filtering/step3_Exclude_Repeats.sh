#!/bin/bash
#SBATCH --job-name=Repeats
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00 
#SBATCH --output=Repeats_%a.out
#SBATCH --error=Repeats_%a.err


#Location of repeats mask 
RPT=/Repeats/Allchr.bed.gz
HEA2=/Repeats/header.txt

INP2=Imputed_MAF_Filtered_Masks.vcf.gz
OUT2=Imputed_MAF_Filtered_Masks_Repeats.vcf.gz

#Annotating variants in repeat regions and removing them
bcftools annotate -a $RPT -m REPEATS=repeats -h $HEA2 -c CHROM,POS,FROM,TO $INP2 | bcftools view --exclude 'REPEATS="repeats"' -Oz -o $OUT2
bcftools index -f $OUT2

