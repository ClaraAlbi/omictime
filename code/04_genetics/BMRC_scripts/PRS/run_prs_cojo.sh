#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J prs

#SBATCH --array=1-22
#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 50G
#SBATCH -c 16

CHROM=${SLURM_ARRAY_TASK_ID}

plink2 \
  --bfile /well/visscher-wray/shared/UKBiobank/v3All_impQC/ukbAll_imp_chr${CHROM}_v3_impQC \
  --score weights_COJO_prots.txt 1 2 header-read \
  --score-col-nums 3-35 \
  --out results/prs/prs_prots_chr${CHROM}