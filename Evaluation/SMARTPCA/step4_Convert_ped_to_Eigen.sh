#!/bin/bash
#SBATCH --job-name=Ped2Eigen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=12:00:00 
#SBATCH --output=Bed2Ped_%a.out
#SBATCH --error=Bed2Ped_%a.err

#Path to software
eigen=/eigensoft/7.2.1/bin/convertf

OUT=SGDP_MERGED

#Build par file
PAR=$OUT\.par
rm $PAR
echo "genotypename:" $OUT\.ped >> $PAR
echo "snpname:" $OUT\.map >> $PAR
echo "indivname:" $OUT\.fam >> $PAR
echo "outputformat: EIGENSTRAT" >> $PAR
echo "genotypeoutname:" $OUT\.geno >> $PAR
echo "snpoutname:" $OUT\.snp >> $PAR
echo "indivoutname:" $OUT\.ind >> $PAR
echo "familynames: NO" >> $PAR

#Converting ped to eigen
$eigen -p $PAR

