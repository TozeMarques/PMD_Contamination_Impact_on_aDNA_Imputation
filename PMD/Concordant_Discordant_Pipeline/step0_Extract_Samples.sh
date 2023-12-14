#!/bin/bash
#SBATCH --job-name=Extract_Sample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00 
#SBATCH --output=Extract_Sample_%a.out
#SBATCH --error=Extract_Sample_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

mkdir -p $PDIR/${SPL}

#Extracting samples from vcf file containing all samples
#Extracting from bcftools calls
bcftools view -s $SPL Bcftools_AllSamples.vcf.gz -Oz -o /${SPL}/${SPL}_Bcf.vcf.gz 
bcftools index -f /${SPL}/${SPL}_Bcf.vcf.gz

#Extracting from Atlas calls
bcftools view -s $SPL Atlas_AllSamples.vcf.gz -Oz -o /${SPL}/${SPL}_Atlas.vcf.gz 
bcftools index -f ${SPL}/${SPL}_Atlas.vcf.gz

#Launch next step
sbatch step1_Extract_columns.sh $1

