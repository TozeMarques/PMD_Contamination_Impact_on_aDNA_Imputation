#!/bin/bash
#SBATCH --job-name=Filter_DP
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00 
#SBATCH --output=Filter_DP_%a.out
#SBATCH --error=Filter_DP_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
SID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Chromosome
CHR=${SLURM_ARRAY_TASK_ID}


RAW=${SID}.chr${CHR}.raw.Q20.q30.calls.vcf.gz
INP=${SID}.chr${CHR}.raw.Q20.q30.1000G.mask.noRepeats.calls.vcf.gz
OUT=${SID}.chr${CHR}.raw.Q20.q30.1000G.mask.noRepeats.DP.calls.vcf.gz

#Depth of coverage measurement
DOC=$($BCFT query -f '%INFO/DP\n' $RAW |  awk 'BEGIN { s = 0; l=0; } { s+=$1; l++; } END { print s/l;}')

#Take as a lower threshold the highest value between on third of the DOC, or 8 DP
LOW=$(python3 limitsDoC.py $DOC | awk '{print $1}')

#Takes as a upper threshold twice the DOC
UPP=$(python3 limitsDoC.py $DOC | awk '{print $2}')

#Excludes positions having less than the LOW value, or more than the UPP value
bcftools filter --exclude "FORMAT/DP<$LOW |  FMT/DP>${UPP}" $INP -Oz -o $OUT
bcftools index -f $OUT
