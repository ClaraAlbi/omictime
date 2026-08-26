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

mkdir -p "$RESDIR"

traits=(
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/res_GWAS.munged.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_M10_TIME.munged.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_L5_TIME.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_MIDP.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_DIURNAL_INACT.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_N_SLEEP_EPISODES.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_DUR.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_DUR_SD.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_EFF.munged.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_v2_GWAS.munged.sumstats.gz"  
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_2018_morning_person.munged.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_b.munged.sumstats.gz"
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/RAgroup.linear.munge.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/RAgroup.logistic.munge.sumstats.gz" 
  "${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/morning_person-5.1.munged.sumstats.gz"
  )

n=${#traits[@]}

for ((i=0; i<n; i++)); do
  for ((j=i+1; j<n; j++)); do

    t1="${traits[i]}"
    t2="${traits[j]}"

    b1=$(basename "$t1" .munged.sumstats.gz)
    b2=$(basename "$t2" .munged.sumstats.gz)

    out="${b1}_vs_${b2}"

    echo "=================================="
    echo "Trait 1: $b1"
    echo "Trait 2: $b2"
    echo "Output : $out"
    echo "=================================="

    python "$LDSC" \
      --rg "${t1},${t2}" \
      --ref-ld-chr "$LDREF/" \
      --w-ld-chr "$LDREF/" \
      --out "$RESDIR/$out"

  done
done
