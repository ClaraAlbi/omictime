#!/bin/bash
#selected traits

#SBATCH -A visscher-wray.prj
#SBATCH -J ldsc

#SBATCH -o logs/my-job-%j.out
#SBATCH -e logs/my-job-%j.err
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 50G
#SBATCH -c 16


RESDIR="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/v3"
LDSC="/users/visscher-wray/uus177/bin/ldsc/ldsc.py"
LDREF="/well/visscher-wray/shared/eur_w_ld_chr"

traits=(
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_M10_TIME.munged.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_L5_TIME.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_MIDP.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_DIURNAL_INACT.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_N_SLEEP_EPISODES.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_DUR.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_DUR_SD.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_EFF.munged.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v4/res_GWAS.munged.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_v2_GWAS.munged.sumstats.gz"  
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/RAgroup.linear.munge.sumstats.gz" 
)


chrono="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/m"

for trait in "${traits[@]}"; do
  name=$(basename "$trait" .munged.sumstats.gz)

  python "$LDSC" \
    --rg ${trait},${chrono} \
    --ref-ld-chr "$LDREF/" \
    --w-ld-chr "$LDREF/" \
    --samp-prev -1,0.58 \
    --pop-prev -1,0.50 \
    --out "$RESDIR/"rg_${name}_vs_morning
done