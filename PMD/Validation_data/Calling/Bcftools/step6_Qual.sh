#!/bin/bash
#SBATCH --job-name=QUAL
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00 
#SBATCH --output=QUAL_%a.out
#SBATCH --error=QUAL_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
SID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Chromosome
CHR=${SLURM_ARRAY_TASK_ID}


INP=${SID}.chr${CHR}.raw.Q20.q30.1000G.mask.noRepeats.DP.calls.vcf.gz
OUT=${SID}.chr${CHR}.raw.Q20.q30.1000G.mask.noRepeats.DP.qual.calls.vcf.gz

#Excludes positions having a value of QUAL less than 30
bcftools filter --exclude "QUAL<30" $INP -Oz -o $OUT
bcftools index -f $OUT
