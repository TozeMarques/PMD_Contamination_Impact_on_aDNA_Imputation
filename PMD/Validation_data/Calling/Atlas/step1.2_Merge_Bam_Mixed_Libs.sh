#!/bin/bash
#SBATCH --job-name=Merge_Libs
#SBATCH --array=1-N
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00 
#SBATCH --output=Merge_Libs_%a.out
#SBATCH --error=Merge_Libs_%a.err
 
#File with 2 columns: bam filename and sample ID
LST=samples.txt

#Sample ID
SPL=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $2}')

#Directory containing the bam files
BAMDIR=/bams

#Merging the paired and unpaired bam files resulting in no longer having mixed libraries
samtools merge $PDIR/${SPL}.bam $PDIR/${SPL}_*
samtools index $PDIR/${SPL}.bam

#Removing temporary files
rm $PDIR/${SPL}_*

