#!/bin/bash
#SBATCH --job-name=Extracting
#SBATCH --array 1-N
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=72:00:00 
#SBATCH --output=Extracting_%a.out
#SBATCH --error=Extracting_%a.err

#File with 2 columns: the fraction of reads to extract, the contamination rate
LST=Fractions_Modern.txt

#Fraction of reads to extract matching the X number of reads removed from the Ancient genome for the desired contamination rate
FRAC=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $1}')

#Contamination rate
ContRate=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $2}')

#Reference genome
FASTA=reference_genome.fasta

#Downsampling/extracting the same amount of X reads extracted from the ancient genome in step1
samtools view -T $FASTA -s $FRAC -bo Contaminating_Reads.bam Modern.bam
samtools index Contaminating_Reads.bam

