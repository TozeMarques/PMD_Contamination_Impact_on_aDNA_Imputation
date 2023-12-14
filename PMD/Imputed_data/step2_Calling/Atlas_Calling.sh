#!/bin/bash
#SBATCH --job-name=Atlas_Call
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00 
#SBATCH --output=Atlas_Call_%a.out
#SBATCH --error=Atlas_Call_%a.err

#File containing the sample IDs
LST=samples.txt 

#Sample ID
SPL=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Chromosome
CHR=${SLURM_ARRAY_TASK_ID}

#Coverage
COV=1.0

#Reference genome
FASTA=reference_genome.fasta

#List of variants in reference panel with the following columns and no header: Chromosome, Position, REF allele, ALT allele
alleles=chr${CHR}positions.sites

#Location of downsampled bam files
BAMDIR0=/bams/cov${COV}
#Downsampled (1x) bam file (input of step 1)
BAM0=$BAMDIR0/${SPL}.bam

#Splitmerged bam file (output of step 1)
BAM=/calls/${SPL}/cov${COV}/${SPL}.chr${CHR}.temp_mergedReads.bam

#Output file of step 3
OUT=/calls/${SPL}/cov${COV}/${SPL}_chr${CHR}_call

mkdir -p /calls/${SPL}/cov${COV}/pmd/

## Step 1: Splitmerge ##
atlas task=splitMerge bam=${BAM0} fasta=${FASTA} chr=$CHR readGroupSettings=${SPL}ReadGroups.txt out=/calls/${SPL}/cov${COV}/${SPL}.chr$CHR.temp

## Step 2: PMD estimation ##
atlas task=PMD bam=${BAM} fasta=${FASTA} length=50 chr=$CHR out=/calls/${SPL}/cov${COV}/pmd/${SPL}.chr${CHR}

## Step 3: Calling ##
atlas task=call method=MLE bam=${BAM} fasta=${FASTA} infoFields=DP chr=$CHR formatFields=GT,DP,GL,PL alleles=$alleles out=${OUT} indName=${SPL} pmdFile=/calls/${SPL}/cov${COV}/pmd/${SPL}.chr${CHR}_PMD_input_Empiric.txt

#Rename
bcftools view ${OUT}_MaximumLikelihood.vcf.gz -Oz -o ${OUT}.vcf.gz
bcftools index -f ${OUT}.vcf.gz

#Remove temporary files
rm ${OUT}_MaximumLikelihood.vcf.gz*
rm /calls/${SPL}/cov${COV}/${SPL}.chr$CHR.temp*

