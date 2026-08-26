#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J extract

#SBATCH -o logs/my-job-%j.out
#SBATCH -e logs/my-job-%j.err
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 50G
#SBATCH -c 16

snpfile="cojo/CA_cojo_snps.txt"

for chr in 1 3 4 5 8 10 11 14 17 19 22; do
  plink \
    --bfile /well/visscher-wray/shared/UKBiobank/v3PPP_impQC/ukbPPP_imp_chr${chr}_v3_impQC \
    --extract "${snpfile}" \
    --recode A \
    --out subset_cojo_CA_chr${chr}
done