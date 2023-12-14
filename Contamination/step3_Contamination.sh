#!/bin/bash
#SBATCH --job-name=Contamination
#SBATCH --array 1-N
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00 
#SBATCH --output=Contamination_%a.out
#SBATCH --error=Contamination_%a.err

#File with 2 columns: the fraction of reads to extract, the contamination rate
LST=Fractions_Modern.txt

#Fraction of reads to extract matching the X number of reads removed from the Ancient genome for the desired contamination rate
FRAC=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $1}')

#Contamination rate
ContRate=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $2}')

#Merging/contaminating the ancient genome with reads from a modern genome
samtools merge -f Ancient_Contaminated_${ContRate}%.bam Ancient_Removed.bam Contaminating_Reads.bam
samtools index Ancient_Contaminated_${ContRate}%.bam

