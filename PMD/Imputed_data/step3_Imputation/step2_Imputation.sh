#!/bin/bash
#SBATCH --array=1-N
#SBATCH --job-name Imputation
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --export=NONE
#SBATCH --mem 4G
#SBATCH --time=2:00:00
#SBATCH --output=Imputation_%a.out
#SBATCH --error=Imputation_%a.err

#Path to GLIMPSE_phase
IMP=/GLIMPSE_v1.1.1/phase/bin/GLIMPSE_phase

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Chunks file with the following columns and no header: Chunk number, Chromosome, Chunk with buffer coordinates, Chunk without buffer coordinates
LST2=referencel.panel.chr${CHR}.txt

DATA=$(cat $LST2 | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

#Chunk number
IDG=$(echo $DATA | cut -d" " -f1)

#Chromosome
CHR=$(echo $DATA | cut -d" " -f2)

#Chunk with buffer coordinates
IREG=$(echo $DATA | cut -d" " -f3)

#Chunk without buffer coordinates
OREG=$(echo $DATA | cut -d" " -f4)

#Reference panel
REF=chr${CHR}.reference.panel.bcf

#Genotype likelihoods
GLS=${ID}_chr${CHR}_GLs.vcf.gz

#Genetic map (if build 37 of the reference genome)
MAP=chr${CHR}.b37.gmap.gz

#Output file
OUT=/${ID}/chr$CHR.reg$IDG.vcf.gz

mkdir -p /${ID}/

#Imputation using GLIMPSE_phase
$IMP --input $GLS --reference $REF --input-region $IREG --output-region $OREG --map $MAP --output $OUT
bcftools index -f $OUT

