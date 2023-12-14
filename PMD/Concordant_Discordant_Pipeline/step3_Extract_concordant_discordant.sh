#!/bin/bash
#SBATCH --job-name=Extract_SNPs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00 
#SBATCH --output=Extract_SNPs_%a.out
#SBATCH --error=Extract_SNPs_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Bcftools GLs vcf file
BCF=${SPL}_Bcf.vcf.gz

#Atlas GLs vcf file
ATLAS=${SPL}_Atlas.vcf.gz

#Using the list of concordantly imputed SNPs, extracts them into a single vcf file and orders them
bcftools view -T conc_list.txt $BCF -Oz -o ${SPL}_Concordant.vcf.gz
bcftools sort ${SPL}_Concordant.vcf.gz -Oz -o ${SPL}_Concordant_sorted.vcf.gz
bcftools index ${SPL}_Concordant_sorted.vcf.gz

#Using the list of discordantly imputed SNPs, unique to Bcftools, extracts them into a single vcf file and orders them
bcftools view -T bcfdisc_list.txt $BCF -Oz -o ${SPL}_Bcf_Discordant.vcf.gz
bcftools sort ${SPL}_Bcf_Discordant.vcf.gz -Oz -o ${SPL}_Bcf_Discordant_sorted.vcf.gz
bcftools index ${SPL}_Bcf_Discordant_sorted.vcf.gz

#Using the list of discordantly imputed SNPs, unique to Atlas, extracts them into a single vcf file and orders them
bcftools view -T atlasdisc_list.txt $ATLAS -Oz -o ${SPL}_Atlas_Discordant.vcf.gz
bcftools sort ${SPL}_Atlas_Discordant.vcf.gz -Oz -o ${SPL}_Atlas_Discordant_sorted.vcf.gz
bcftools index ${SPL}_Atlas_Discordant_sorted.vcf.gz

#Launch next step
sbatch step4_Change_header.sh $1

