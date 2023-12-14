#!/bin/bash
#SBATCH --job-name=Extract_Cols
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00 
#SBATCH --output=Extract_Cols_%a.out
#SBATCH --error=Extract_Cols_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Bcftools GLs vcf file
BCF=${SPL}_Bcf.vcf.gz

#Atlas GLs vcf file
ATLAS=${SPL}_Atlas.vcf.gz

#Extracting the %CHROM\t%POS\t[%GT]\n columns from both Bcftools and Atlas GLs files
bcftools query -f '%CHROM\t%POS\t[%GT]\n' $BCF > ${SPL}_Bcf_Extracted_GTs.txt
bcftools query -f '%CHROM\t%POS\t[%GT]\n' $ATLAS > ${SPL}_Atlas_Extracted_GTs.txt

#Launch next step
sbatch step2_python_CompareGTs.sh $1
