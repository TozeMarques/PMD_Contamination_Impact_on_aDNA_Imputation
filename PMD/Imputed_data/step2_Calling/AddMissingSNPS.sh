#!/bin/bash
#SBATCH --job-name=FixMissSNPs
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00 
#SBATCH --output=FixMissSNPs_%a.out
#SBATCH --error=FixMissSNPs_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
SPL=$(cat samples | head -n $1 | tail -n 1 | awk '{print $1}')

#Chromosome
CHR=${SLURM_ARRAY_TASK_ID}

#Reference file
REF=reference.vcf.gz

#Atlas file
ATLAS=/calls/${SPL}/cov1.0/${SPL}_chr${CHR}_atlas_call.vcf.gz

#Output files
OUT=${SPL}_chr${CHR}_atlas_call_temp.vcf.gz
OUT2=${SPL}_chr${CHR}_atlas_call.vcf.gz

#Adding missing positions into Atlas vcf file as "placeholder" SNPs
python AddMissingSNPS.py $REF $ATLAS $OUT
bcftools view $OUT -Oz -o $OUT2
bcftools index -f $OUT2

#Remove temporary file
rm $OUT

