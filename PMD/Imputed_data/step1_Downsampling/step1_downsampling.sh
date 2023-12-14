#!/bin/bash
#SBATCH --job-name=Downsampling1
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00 
#SBATCH --output=Downsampling1_%a.out
#SBATCH --error=Downsampling1_%a.er

#File with 2 columns: bam filename and sample ID
LST=bamlist.txt

#Bam file
BAM=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Sample ID
SPL=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $2}')

#Directory containing the bam files
BAMDIR=/bams

#Reference genome
FASTA=reference_genome.fasta

#Chromosome
CHR=${SLURM_ARRAY_TASK_ID}


#Position files
#List of variants in reference panel with the following header: #CHROM\tPOS\tID\tREF\tALT
#Can be generated with bcftools view -G reference.panel.genotypes.bcf
VPOS=chr${CHR}.referencePanel.vcf.gz

#List of variants in reference panel with the following columns and no header: CHROM\tPOS\tREF,ALT
TSV=chr${CHR}.referencePanel.tsv.gz


#Output VCF files
SPL1=calls/chr${CHR}/${SPL}.raw.calls.spl
OUT=calls/chr${CHR}/${SPL}.raw.calls.vcf.gz
TMP=calls/chr${CHR}/${SPL}.raw.calls.tmp.vcf.gz
LOG=calls/chr${CHR}/${SPL}.raw.calls.log

#Call genotypes using bcftools
echo ${BAMDIR}/${BAM} $SPL > $SPL1
bcftools mpileup -f $FASTA -I -E -a 'FORMAT/DP' --ignore-RG -T $VPOS -r $CHR ${BAMDIR}/${BAM} | bcftools call -Aim -C alleles -T $TSV -Oz -o $OUT
bcftools reheader -s $SPL1 -o $TMP $OUT
mv $TMP $OUT
bcftools index -f $OUT

