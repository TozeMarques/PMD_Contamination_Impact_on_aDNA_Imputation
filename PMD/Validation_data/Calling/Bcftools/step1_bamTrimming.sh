#!/bin/bash
#SBATCH --job-name=Trimming
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=16:00:00 
#SBATCH --output=Trimming_%a.out
#SBATCH --error=Trimming_%a.err

#Path to software
BIN=/work/FAC/FBM/DBC/amalaspi/popgen/bmota/binaries/bamUtil/bam

#File with 2 columns: bam filename and sample ID
LST=bamlist.txt

#Bam file
BAM=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Sample ID
SPL=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $2}')

INP=$BAM

OUT=bams/${SPL}.trim10.bam

#Trimming 10 bases
$BIN trimBam $INP $OUT 10
samtools index $OUT
