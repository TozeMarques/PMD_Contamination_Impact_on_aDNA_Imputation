#!/bin/bash
#SBATCH --job-name=Removing
#SBATCH --array 1-N
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00 
#SBATCH --output=Removing_%a.out
#SBATCH --error=Removing_%a.err

#File with 2 columns: the fraction of reads to remove, the contamination rate
LST=Fractions_Ancient.txt

#Fraction of reads to remove
FRAC=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $1}')

#Contamination rate
ContRate=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $2}')

#Reference genome
FASTA=reference_genome.fasta

#Downsampling/removing an X amount of reads
samtools view -T $FASTA -s $FRAC -bo Ancient_Removed.bam Ancient.bam
samtools index Ancient_Removed.bam

