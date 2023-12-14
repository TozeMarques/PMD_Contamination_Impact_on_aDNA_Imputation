#!/bin/bash
#SBATCH --job-name=Ped2Eigen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=16:00:00 
#SBATCH --output=Bed2Ped_%a.out
#SBATCH --error=Bed2Ped_%a.err

#Path to software
eigen=/eigensoft/7.2.1/bin/smartpca

INP=SGDP_MERGED
OUT=/Results/SGDP_MERGED

#Build PAR file
PAR=$OUT\.par
rm $PAR
echo "genotypename:" $INP\.geno >> $PAR
echo "snpname:" $INP\.snp >> $PAR
echo "indivname:" $INP\.ind >> $PAR
echo "evecoutname:" $OUT\.evec >> $PAR
echo "evaloutname:" $OUT\.eval >>$PAR 
echo "poplistname: pops.txt" >> $PAR #pops.txt contains the single word "Control"
echo "lsqproject: YES" >> $PAR
echo "numoutevec: 10" >> $PAR
echo "numthreads: 1" >> $PAR
echo "outliermode: 2" >> $PAR

#Running SmartPCA
$eigen -p $PAR

