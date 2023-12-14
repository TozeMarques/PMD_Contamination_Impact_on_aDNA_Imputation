#!/bin/bash

#File containing the sample IDs
LST=samples.txt

#Sample ID
ID=$(cat $LST | head -n $1 | tail -n 1 | awk '{print $1}')

#Directory where lists will be saved
DIR=lists/${ID}

mkdir -p ${DIR}

#File name
LST=${DIR}/${ID}_concordance.lst

:>$LST1
for CHR in {1..22}; do

	#Reference panel
	REF=chr${CHR}.reference.panel.bcf

	#Validation data set
	VAL=Validation.vcf.gz

	#Imputed data set
	IMPBCF=Imputed.vcf.gz
	
	#Making a file containing 4 columns: Chromosome, Reference panel file path, Validation date file path, Imputed data file path
	echo "$CHR $REF $VAL $IMPBCF" >> $LST

done

