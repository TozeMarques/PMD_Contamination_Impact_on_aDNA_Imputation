#!/bin/bash
#SBATCH --job-name=Vcf2Bed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00 
#SBATCH --output=Vcf2Bed_%a.out
#SBATCH --error=Vcf2Bed_%a.err

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Vcf file (Imputed or Validation vcf file)
FILE=${SPL}_Merged.vcf.gz

#Converting vcf file (Imputed or Validation file) into bed
$PLINK --vcf $FILE --make-bed --double-id --out ${SPL}_Merged

#Changing the format of the 2nd column of the .bim file
awk '{OFS="\t"; $2=$1"_"$4; print}' ${SPL}_Merged.bim > ${SPL}_Merged_modified.bim
mv ${SPL}_Merged_modified.bim ${SPL}_Merged.bim

#Changing the 6th column of the fam file to add a 2 (2=Case sample)
awk '{ $6=2; print}' ${SPL}_Merged.fam > ${SPL}_Merged_modified.fam 
mv ${SPL}_Merged_modified.fam ${SPL}_Merged.fam

