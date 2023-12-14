#!/bin/bash
#SBATCH --job-name=Header
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00 
#SBATCH --output=Header_%a.out
#SBATCH --error=Header_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Concordant and discordant files
CON=${SPL}_Concordant_sorted.vcf.gz
BCF=${SPL}_Bcf_Discordant_sorted.vcf.gz
ATLAS=${SPL}_Atlas_Discordant_sorted.vcf.gz

#Extracts the the header of the files and appends _Imp_... to the sample ID
bcftools query -l $CON | awk '{print $1 "_Imp_Con"}' > Con_samples.txt
bcftools query -l $BCF | awk '{print $1 "_Imp_Bcf_Disc"}' > Bcf_samples.txt
bcftools query -l $ATLAS | awk '{print $1 "_Imp_Atlas_Disc"}' > Atlas_samples.txt

#Changes the header of the files
bcftools reheader -s Con_samples.txt -o ${SPL}_Concordant_sorted_reheader.vcf.gz $CON
bcftools reheader -s Bcf_samples.txt -o ${SPL}_Bcf_Discordant_sorted_reheader.vcf.gz $BCF
bcftools reheader -s Atlas_samples.txt -o ${SPL}_Atlas_Discordant_sorted_reheader.vcf.gz $ATLAS

#Indexing
bcftools index ${SPL}_Concordant_sorted_reheader.vcf.gz
bcftools index ${SPL}_Bcf_Discordant_sorted_reheader.vcf.gz 
bcftools index ${SPL}_Atlas_Discordant_sorted_reheader.vcf.gz 

#Merges the concordant and the two discordant files into a single vcf file
bcftools merge ${SPL}_Concordant_sorted_reheader.vcf.gz ${SPL}_Bcf_Discordant_sorted_reheader.vcf.gz  ${SPL}_Atlas_Discordant_sorted_reheader.vcf.gz -Oz -o ${SPL}_Merged.vcf.gz
bcftools index ${SPL}_Merged.vcf.gz

