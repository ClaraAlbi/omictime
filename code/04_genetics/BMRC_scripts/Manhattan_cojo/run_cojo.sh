#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J cojo

#SBATCH --array=1-22
#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 50G
#SBATCH -c 16

# Map the array task ID to the chromosome number
chr=${SLURM_ARRAY_TASK_ID}

# Path to GCTA binary
GCTA_BIN=/users/visscher-wray/uus177/bin/gcta-1.95.0-linux-kernel-3-x86_64/gcta64


cojo_f="data/GWAS_res_v5/res_eur_chr${chr}.fastGWA"

tmp=$(mktemp --tmpdir="$(dirname "$f")" "$(basename "$f").tmp.XXXXXXXX") || exit 1

awk -v OFS='\t' '{
    print $2, $4, $5, $7, $8, $9, $10, $6
  }' "$cojo_f" > "$tmp"

head $tmp

# Run COJO for this chromosome
${GCTA_BIN} \
  --bfile /well/visscher-wray/shared/UKBiobank/v3PPP_impQC/ukbPPP_imp_chr${chr}_v3_impQC \
  --chr ${chr} \
  --maf 0.01 \
  --cojo-file $tmp \
  --cojo-slct \
  --out results/cojo/res_v5_chr${chr}

rm -f "$tmp"