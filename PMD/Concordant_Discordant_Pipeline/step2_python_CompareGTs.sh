#!/bin/bash
#SBATCH --job-name=python_CompareGTs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00 
#SBATCH --output=python_CompareGTs_%a.out
#SBATCH --error=python_CompareGTs_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Python script to compare the GT of each SNP between the two files. If the same -> writes to a file, if different -> write to 2 different files
python CompareGTs.py ${SPL}_Bcf_Extracted_GTs.txt ${SPL}_Atlas_Extracted_GTs.txt

#Launch next step
sbatch step3_Extract_concordant_discordant.sh $1
