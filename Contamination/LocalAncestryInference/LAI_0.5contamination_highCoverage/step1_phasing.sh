#!/bin/bash
#SBATCH --job-name phase
#SBATCH --time=1:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --export=NONE
#SBATCH --mem 8G
#SBATCH --output=phase.out
#SBATCH --error=phase.er

#working directory
PDIR1=PATH0

#CHR
CHR=1
#Sample ID
SPL=Loschbour_Dinka

#Genetic map
MAP=PATH1/GMAP_shapeit4/chr${CHR}.b37.gmap
REF=PATH2/1000Genomes_umich/noSingletons/1000GP_nygc_umich_chr${CHR}_snps_only_genotypes.bcf
INP=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.noRepeats.DP.qual.calls.bcf

OUT1=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.noRepeats.DP.qual.noMissingness.bcf
OUT=${PDIR1}/chr${CHR}.${SPL}.Q20.q30.1000G.mask.noRepeats.DP.qual.phased.bcf

BIN=PATH3/SHAPEIT5/phase_common_static

# Remove SNPs with missing data
bcftools view $INP -e 'F_MISSING>0' | bcftools +fill-tags - -Ob -o $OUT1 -- -t AN,AC
bcftools index -f $OUT1

#Phase with SHAPEIT5
$BIN --input $OUT1 -H $REF --filter-maf 0 --region ${CHR} --map $MAP --output $OUT --thread 4
