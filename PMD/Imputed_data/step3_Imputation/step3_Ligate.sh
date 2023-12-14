#!/bin/bash
#SBATCH --job-name=Ligate
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00 
#SBATCH --output=Ligate_%a.out
#SBATCH --error=Ligate_%a.err

#Path to GLIMPSE_ligate
BIN=/GLIMPSE_v1.1.1/ligate/bin/GLIMPSE_ligate

#Chromosome
CHR=$SLURM_ARRAY_TASK_ID

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

mkdir -p /step3_ligate/lists/${ID}/
mkdir -p /step3_ligate/ligated/${ID}/

#List containing all imputed chunks of a chromosome
LST2=/step3_ligate/lists/${ID}/chr$CHR.list.txt
ls -1v /${ID}/chr$CHR.reg*.vcf.gz > $LST2

#Outputs
OUT=/step3_ligate/ligated/${ID}/chr$CHR.ligated.vcf.gz
OUL=/step3_ligate/ligated/${ID}/chr$CHR.ligated.log

#Ligating all the chunks of a chromosome with GLIMPSE_ligate
$BIN --input $LST --output $OUT --log $OUL
bcftools index -f $OUT

