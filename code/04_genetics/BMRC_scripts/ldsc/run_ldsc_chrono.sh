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


python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
  --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v4/res_GWAS.munged.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_2018_morning_person.munged.sumstats.gz  \
  --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
  --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
  --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/res_GWAS_Jones_2018_morning_person
  
python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
  --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v4/res_GWAS.munged.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_v2_GWAS.munged.sumstats.gz  \
  --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
  --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
  --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/res_GWAS_chrono_v2
  