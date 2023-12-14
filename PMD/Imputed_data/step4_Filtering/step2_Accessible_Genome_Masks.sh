#!/bin/bash
#SBATCH --job-name=Masks
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00 
#SBATCH --output=Masks_%a.out
#SBATCH --error=Masks_%a.err


#Location of 1000 Genomes masks in curnagl
MSK=/Accessible_Genome_Masks/Allchr.bed.gz
HEA=/Accessible_Genome_Masks/header.txt

INP1=Imputed_MAF_Filtered.vcf.gz
OUT1=Imputed_MAF_Filtered_Masks.vcf.gz

#Annotating and keeping only the 1000G Accessible Genome strict Masks
bcftools annotate -a $MSK -m MASK=strict -h $HEA -c CHROM,POS,FROM,TO $INP1 | bcftools view -i 'INFO/MASK="strict"' -Oz -o $OUT1
bcftools index -f $OUT1

