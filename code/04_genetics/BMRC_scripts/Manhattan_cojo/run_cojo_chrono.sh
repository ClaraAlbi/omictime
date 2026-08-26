#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J cojo

#SBATCH --array=1-21
#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 50G
#SBATCH -c 16

chr=${SLURM_ARRAY_TASK_ID}

GCTA_BIN=/users/visscher-wray/uus177/bin/gcta-1.95.0-linux-kernel-3-x86_64/gcta64

cojo_f="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/morning_person-5.1.cojo"

head $cojo_f

# Run COJO for this chromosome
${GCTA_BIN} \
  --bfile /well/visscher-wray/shared/UKBiobank/v3PPP_impQC/ukbPPP_imp_chr${chr}_v3_impQC \
  --chr ${chr} \
  --maf 0.01 \
  --cojo-file ${cojo_f} \
  --cojo-slct \
  --out results/cojo/morning_person-5.1_chr${chr}