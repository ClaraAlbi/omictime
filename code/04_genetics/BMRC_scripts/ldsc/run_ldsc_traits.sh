#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J cojo

#SBATCH -o logs/my-job-%j.out
#SBATCH -e logs/my-job-%j.err
#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH --mem 50G
#SBATCH -c 16

chrono=(${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_UKB_GWAS_merged_munged.fastGWA.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res/res_GWAS.munged.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_cond/res_GWAS_cond_munged.fastGWA.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_pca/res_GWAS_pPCA.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_2018_morning_person.munged.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_UKB.sumstats.gz
)
# #chrono=(${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_cond/res_GWAS_cond_munged.fastGWA.sumstats.gz)
# #chrono=(${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_pca/res_GWAS_pPCA.sumstats.gz)
# chrono=(${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_2018_morning_person.munged.sumstats.gz ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_UKB.sumstats.gz)

#psych=(BIP.sumstats.gz pgc-mdd2024_Clin_eur.sumstats.gz pgc-mdd2024_EHR_eur.sumstats.gz pgc-mdd2024_Quest_eur.sumstats.gz SCZ.sumstats.gz)
psych=(CAD.sumstats.gz)

for c in "${chrono[@]}"; do
  for p in "${psych[@]}"; do
    cbase=$(basename "$c"); cbase="${cbase%%.*}"
    pbase=$(basename "$p"); pbase="${pbase%%.*}"
    out="${cbase}_vs_${pbase}"
    
    echo $p
    echo $c
    echo $out
    python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
      --rg ${c},${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_other/${p} \
      --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
      --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
      --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/${out}
  done
done