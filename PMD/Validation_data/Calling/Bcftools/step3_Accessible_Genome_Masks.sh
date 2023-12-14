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

#File containing the sample IDs
LST=samples.txt

#Sample ID
SID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Chromosome
CHR=${SLURM_ARRAY_TASK_ID}

#Location of 1000 Genomes masks in curnagl
MSK=Accessible_Genome_Masks/chr${CHR}.bed.gz
HEA=Accessible_Genome_Masks/header.txt

INP=${SID}.chr${CHR}.raw.Q20.q30.calls.vcf.gz
OUT=${SID}.chr${CHR}.raw.Q20.q30.1000G.mask.calls.vcf.gz

#Annotating and keeping only the 1000G Accessible Genome strict Masks
bcftools annotate -a $MSK -m MASK=strict -h $HEA -c CHROM,POS,FROM,TO $INP | bcftools view -i 'INFO/MASK="strict"' -Oz -o $OUT
bcftools index -f $OUT
