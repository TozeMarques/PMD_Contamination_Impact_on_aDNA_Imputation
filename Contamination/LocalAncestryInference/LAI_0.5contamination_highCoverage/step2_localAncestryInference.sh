#!/bin/bash
#SBATCH --job-name rfmix
#SBATCH --array 1
#SBATCH --time=5:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --export=NONE
#SBATCH --mem 60G
#SBATCH --output=rfmix.out
#SBATCH --error=rfmix.er

#This script allows inferring local ancestry with RFmix for the contaminated genome of Loschbour with 50% of Dinka reads (high-coverage version)

#chromosome
CHR=1
# Sample ID
SPL=Loschbour_Dinka

#binary file for RFmix
BIN=PATH1/rfmix/rfmix

#project directory
PDIR1=PATH0
#output directory
ODIR=${PDIR1}/local_ancestry
mkdir -p $ODIR

#input phased genotype calls 
INP=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.noRepeats.DP.qual.phased.bcf
#reference panel for local ancestry inference
REF=PATH2/localAncestry_panel/1000GP_CEU_YRI_LWK_chr${CHR}_snps_only_genotypes.bcf
#populations file
POP=PATH2/inds_pops_YRI_CEU_LWK_localAncestry.txt
#genetic map
MAP=PATH3/GMAP_shapeit4/allChrs.b37.rfmix.gmap
#output file prefix
OUT=${ODIR}/chr${CHR}.${SPL}.highCoverage.CEU.YRI.LWK

# Infer local ancestry with RFMix
$BIN -f $INP -r $REF -m $POP -g $MAP -o ${OUT} --chromosome=${CHR} --n-threads=10
