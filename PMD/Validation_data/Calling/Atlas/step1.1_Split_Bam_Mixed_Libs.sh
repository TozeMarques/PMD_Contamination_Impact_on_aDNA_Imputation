#!/bin/bash
#SBATCH --job-name=Split_Mixed
#SBATCH --array=1-N
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00 
#SBATCH --output=Split_Mixed_%a.out
#SBATCH --error=Split_Mixed_%a.err


#File with 2 columns: bam filename and sample ID
LST=samples.txt

#Directory containing the bam files
BAMDIR=/bams

#Bam file name
BAM=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Sample ID
SPL=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $2}')

#File containing all the libraries names of the mixed library samples
$LST2=${SPL}_Libs.txt

#Library name
LIB=$(cat $LST2 | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $1}')

#Splitting the paired reads of the mixed library into a bam file
samtools view -b -r $LIB -f 1 -F 0 $BAMDIR/$BAM > ${SPL}_${LIB}_paired.bam

#Splitting the unpaired reads of the mixed library into a bam file
samtools view -b -r $LIB -f 0 -F 1 $BAMDIR/$BAM > ${SPL}_${LIB}_unpaired.bam

