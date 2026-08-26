#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J prs

#SBATCH --array=1-21
#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 50G
#SBATCH -c 16

CHROM=${SLURM_ARRAY_TASK_ID}

plink --bfile "/well/visscher-wray/shared/UKBiobank/v3All_impQC/ukbAll_imp_chr${CHROM}_v3_impQC" \
      --score results/sbayesrc/gctb.snpRes 2 5 8 header \
      --out results/prs/prs_CA_norep_chr${CHROM}
