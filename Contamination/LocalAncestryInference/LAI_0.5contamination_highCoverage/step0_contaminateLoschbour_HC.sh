#!/bin/bash
#SBATCH --job-name=sample
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=1:30:00
#SBATCH --output=sample.out
#SBATCH --error=sample.er

# as specified for a slurm system
module load gcc; module load samtools/1.15.1;

# Ancient genome bam file (original high-coverage Loschbour 18x)
BAM=PATH0/Loschbour.hg19_1000g.bam
# Contamination source: Dinka-3 genome (SGDP)
BAMD=PATHSGDP/SS6004480.srt.aln.bam
# Working directory
PDIR1=PATH1/step1_local_ancestry/contaminatedHighCoverage/
TMP=${PDIR1}/Loschbour_subsampling0.5.bam
TMP2=${PDIR1}/Dinka-3.chr1.bam
TMP3=${PDIR1}/Dinka-3.chr1.subsampled.bam
OUTBAM=${PDIR1}/chr1.Loschbour_Dinka-3_0.5_hc.bam

# Subsample high-coverage ancient genome to keep half of the reads
echo "Downsampling Loschbour genome with subsampling s=1.5 where s=seed.fraction (chr1)" 
samtools view -s 1.5 $BAM 1 -bo $TMP
samtools index $TMP

echo "Downsampling is done"

#Extract chromosome 1 from Dinka-3 genome

samtools view $BAMD 1 -bo $TMP2
samtools index $TMP2

#Number of reads mapping to chr1 in Dinka-3 genome
ND=$(samtools view -c $TMP2)
#Number of reads mapping to chr1 in subsampled Loschbour genome
NL=$(samtools view -c $TMP)

# Calculate fraction of reads to subsample from Dinka-3 genome to achieve 50% contamination
FRAC=$(python -c "print( 1+$NL / float($ND) )")

# Subsample Dinka-3 genome to extract the number of reads required to achieve 50% contamination
echo "Downsampling genome with subsampling s="${FRAC}" where s=seed.fraction" 
samtools view -s $FRAC $TMP2 1 -bo $TMP3
samtools index $TMP3
echo "Downsampling is done"

# Merge reads from Loschbour with reads from Dinka-3

echo "START MERGING BAMS"
samtools merge -o $OUTBAM $TMP $TMP3
samtools index -f $OUTBAM

echo "FINAL BAM IS READY"

##----- Call genotypes for contaminated bam file --------
#Chromosome
CHR=1
#Sample ID
SPL=Loschbour_Dinka

#reference genome
FASTA=PATH2/reference_human/hs.build37.1/hs.build37.1.fa

#Just SNPs 
VTYPE=snps_only
REF=1000GP_nygc_umich
#Position file	
VPOS=PATH3/1000Genomes_umich/noSingletons/1000GP_nygc_umich_chr${CHR}_snps_only_sites.vcf.gz
TSV=PATH3/1000Genomes_umich/noSingletons/1000GP_nygc_umich_chr${CHR}_snps_only_sites.tsv.gz
	
#Output VCF and LOG files
SPL1=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.calls.spl
OUT=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.calls.vcf.gz
TMP=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.calls.tmp.vcf.gz
LOG=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.calls.log

echo "START GENOTYPE CALLING"

#Call genotypes using bcftools
echo ${OUTBAM} $SPL > $SPL1
bcftools mpileup -f $FASTA -I -E -a 'FORMAT/DP' --ignore-RG -T $VPOS  -Q 20 -q 30  -r $CHR ${OUTBAM} | bcftools call -Aim -C alleles -T $TSV -Oz -o $OUT
bcftools reheader -s $SPL1 -o $TMP $OUT
mv $TMP $OUT
bcftools index -f $OUT
echo $CHR $OUTBAM $? >> $LOG

echo "GENOTYPE CALLING COMPLETED"

# Apply 1000G accessible mask

MSK=PATH4/accessible_genome_masks/chr${CHR}.bed.gz
HEA=PATH4/accessible_genome_masks/header.txt

INP1=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.calls.vcf.gz
OUT1=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.calls.vcf.gz

echo "START APPLYING 1000G MASK"
bcftools annotate -a $MSK -m MASK=strict -h $HEA -c CHROM,POS,FROM,TO $INP1 | bcftools view -i 'INFO/MASK="strict"' -Oz -o $OUT1
bcftools index -f $OUT1

echo "1000G mask completed"

# remove repeat regions

RPT=PATH5/Repeats/by_chr_sorted/chr${CHR}.bed.gz
HEA2=PATH5/Repeats/by_chr_sorted/header.txt

INP2=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.calls.vcf.gz
OUT2=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.noRepeats.calls.vcf.gz

echo "START EXCLUDING REPEAT REGIONS"
bcftools annotate -a $RPT -m REPEATS=repeats -h $HEA2 -c CHROM,POS,FROM,TO $INP2 | bcftools view --exclude 'REPEATS="repeats"' -Oz -o $OUT2
bcftools index -f $OUT2

echo "REPEAT REGIONS REMOVED"

# Filter for depth
RAW=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.calls.vcf.gz
INP3=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.noRepeats.calls.vcf.gz
OUT3=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.noRepeats.DP.qual.calls.bcf

echo "CALCULATE DOC"
DOC=$(bcftools query -f '%INFO/DP\n' $RAW |  awk 'BEGIN { s = 0; l=0; } { s+=$1; l++; } END { print s/l;}')

echo $DOC; 

LOW=$(python3 limitsDoC.py $DOC | awk '{print $1}')
UPP=$(python3 limitsDoC.py $DOC | awk '{print $2}')
echo $LOW $UPP

echo "REMOVE SITES WITH EXTREME DEPTH AND QUAL<30"
bcftools filter --exclude "FORMAT/DP<$LOW |  FMT/DP>${UPP}" $INP3 | bcftools filter --exclude "QUAL<30" | bcftools view -e 'F_MISSING<0' -Ob -o $OUT3
bcftools index -f $OUT3
echo "DEPTH AND QUAL FILTERING DONE"

#remove intermediate files
rm ${PDIR1}/*vcf.gz*
