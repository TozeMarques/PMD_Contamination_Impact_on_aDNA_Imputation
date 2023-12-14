#!/bin/bash
#SBATCH --job-name=Chunking
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00 
#SBATCH --output=Chunking_%a.out
#SBATCH --error=Chunking_%a.err


#Chromosome
CHR=${SLURM_ARRAY_TASK_ID}

#Path to GLIMPSE_chunk 
BIN=/GLIMPSE-v1.1.1/chunk/bin/GLIMPSE_chunk

#Reference panel
REF=chr${CHR}.reference.panel.bcf

#Output file, chunk coordinates
OUT=/coordinates/referencel.panel.chr${CHR}.txt

#Generate chunks, with window size of 1Mb with a buffer of 200kb size
$BIN --input $REF --region chr$CHR --window-size 1000000 --buffer-size 200000  --output $OUT

