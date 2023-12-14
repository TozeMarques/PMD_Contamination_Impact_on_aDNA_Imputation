#!/bin/bash
#SBATCH --job-name=Validation
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=16:00:00 
#SBATCH --output=Validation_%a.out
#SBATCH --error=Validation_%a.err

#File with 2 columns: bam filename and sample ID
LST=bamlist.txt

#Sample ID
SID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $2}')

#Directory containing the bam files
BAMDIR=/bams

#Bam file
BAM=${SPL}.trim10.bam

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
SID1=chr${CHR}/${SID}.raw.Q20.q30.calls.spl
OUT=chr${CHR}/${SID}.raw.Q20.q30.calls.vcf.gz
TMP=chr${CHR}/${SID}.raw.Q20.q30.calls.tmp.vcf.gz

#Call genotypes using bcftools with 3 filters (-Q 20 -q 30 -C 50)
echo ${BAMDIR}/${BAM} $SID > $SID1
bcftools mpileup -f $FASTA -I -E -a 'FORMAT/DP' --ignore-RG -T $VPOS -Q 20 -q 30 -C 50 -r $CHR ${BAMDIR}/${BAM} | bcftools call -Aim -C alleles -T $TSV -Oz -o $OUT
bcftoolsreheader -s $SID1 -o $TMP $OUT
mv $TMP $OUT
bcftools index -f $OUT


