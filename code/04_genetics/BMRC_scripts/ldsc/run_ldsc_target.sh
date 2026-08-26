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

pref="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/"
RESDIR="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/ldsc/"

mkdir -p "$RESDIR"

chrono=(${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/res_GWAS.munged.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_v2_GWAS.munged.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_L5_TIME.munged.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_M10_TIME.munged.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/ACC_SLEEP_MIDP.munged.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_b.munged.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/RAgroup_lin.linear.munge.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/RAgroup_log.logistic.munge.sumstats.gz 
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/morning_person-5.1.munged.sumstats.gz 
 )
 
chrono=(
${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/res_GWAS.munged.sumstats.gz 
 )

psych=(
GWAS_other/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv_tidy.ma.full.sumstats.gz
GWAS_MDD/bip2024.munged.sumstats.gz
GWAS_MDD/pgc-mdd2024_Clin_eur.sumstats.gz
GWAS_MDD/pgc-mdd2024_EHR_eur.sumstats.gz
GWAS_MDD/pgc-mdd2024_Quest_eur.sumstats.gz
23andme/Jansen_2018_insomnia_broad.munged.sumstats.gz
GWAS_other/AD_sumstats_Jansenetal_2019sept_tidy.ma.full.sumstats.gz
arth/GCST90132223_buildGRCh37.munged.sumstats.gz 
selected/BP-ICE_EUR_DBP_transformed_15-04-2020.munge.sumstats.gz
selected/BP-ICE_EUR_HTN_15-04-2020.munge.sumstats.gz
selected/BP-ICE_EUR_PP_transformed_15-04-2020.munge.sumstats.gz
selected/BP-ICE_EUR_SBP_transformed_15-04-2020.munge.sumstats.gz
GWAS_other/GCST90132314_buildGRCh37.tsv.sumstats.gz
GWAS_other/EUR_Metal_LDSC-CORR_Neff.v2.txt_tidy.ma.full.sumstats.gz
GWAS_other/mtDNA_CN.ALLm2.bgen.stats_imp.ma.sumstats.gz
GWAS_other/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt_imp.ma.sumstats.gz
GWAS_other/MAGIC_FastingGlucose.txt_imp.ma.sumstats.gz
GWAS_other/Insomnia_01/1200_imp.ma.sumstats.gz
GWAS_other/GCST012229.tsv_imp.ma.sumstats.gz
)

psych=(
  
)

# Loop through combinations
for c in "${chrono[@]}"; do
  for p in "${psych[@]}"; do

    # Skip identical pairs (same file)
    if [[ "$c" == "$p" ]]; then
        continue
    fi

    # Strip extensions for output naming
    cbase=$(basename "$c"); cbase="${cbase%%.*}"
    pbase=$(basename "$p"); pbase="${pbase%%.*}"

    out="${cbase}_vs_${pbase}"

    echo "-----------------------------"
    echo "Chrono file: $c"
    echo "Psych file : $p"
    echo "Output tag : $out"
    echo "-----------------------------"

    python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
      --rg "${c},${pref}${p}" \
      --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
      --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
      --out "${RESDIR}/${out}"

  done
done
