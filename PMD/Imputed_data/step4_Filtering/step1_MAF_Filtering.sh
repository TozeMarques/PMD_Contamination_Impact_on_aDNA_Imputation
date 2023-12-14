#!/bin/bash
#SBATCH --job-name=MAF_Filter
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00 
#SBATCH --output=MAF_Filter_%a.out
#SBATCH --error=MAF_Filter_%a.err


#MAF-filtering with vcftools
vcftools --vcf Imputed.vcf --maf 0.05 --recode --out Imputed_MAF_Filtered

#Compressing the output vcf file
bgzip -c Imputed_MAF_Filtered.recode.vcf > Imputed_MAF_Filtered.vcf.gz

bcftools index Imputed_MAF_Filtered.vcf.gz
