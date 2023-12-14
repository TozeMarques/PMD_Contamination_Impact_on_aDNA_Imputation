#!/bin/bash
#SBATCH --job-name=Downsampling2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00 
#SBATCH --output=Downsampling2_%a.out
#SBATCH --error=Downsampling2_%a.err


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


#Extract the depth per position across all chromosomes for the original genome
:> ${SPL}.cov.gz
for CHR in {1..22}; do
	bcftools query -f '%INFO/DP\n' /chr${CHR}/${SPL}.raw.calls.vcf.gz | gzip -c >> ${SPL}.cov.gz
done

#Downsample & calling
COV=1
mkdir -p /calls/cov${COV}/bams
	
#Downsampled bam file name
DS_BAM=/calls/cov${COV}/bams/${SPL}.bam

#Compute fraction of reads from INFO/DP
FRAC=$(zcat ${SPL}.cov.gz | awk -v c=$COV 'BEGIN { s = 0; l=0; } { s+=$1; l++; } END { print 1+c*l/s; }')	
echo $FRAC > /calls/cov${COV}/bams/${SPL}.frac
echo "Fraction=" $FRAC

#Downsampling sequencing reads
echo "Downsampling Reads Coverage = " $COV "x ********************************"
echo "DS_BAM = " $DS_BAM  
	
#Downsampling
samtools view -T $FASTA -s $FRAC -bo $DS_BAM ${BAMDIR}/${BAM}
samtools index $DS_BAM


