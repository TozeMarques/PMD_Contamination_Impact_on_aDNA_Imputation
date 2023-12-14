#!/bin/bash
#SBATCH --job-name=Repeats
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00 
#SBATCH --output=Repeats_%a.out
#SBATCH --error=Repeats_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
SID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Chromosome
CHR=${SLURM_ARRAY_TASK_ID}

#Location of repeats mask 
RPT=/Repeats/chr${CHR}.bed.gz
HEA2=/Repeats/header.txt

INP=${SID}.chr${CHR}.raw.Q20.q30.1000G.mask.calls.vcf.gz
OUT=${SID}.chr${CHR}.raw.Q20.q30.1000G.mask.noRepeats.calls.vcf.gz

#Annotating variants in repeat regions and removing them
bcftools annotate -a $RPT -m REPEATS=repeats -h $HEA -c CHROM,POS,FROM,TO $INP | bcftools view --exclude 'REPEATS="repeats"' -Oz -o $OUT
bcftools index -f $OUT

