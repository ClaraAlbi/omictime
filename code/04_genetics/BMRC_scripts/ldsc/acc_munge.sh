#!/bin/bash

OUTDIR=${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc

for TRAIT in \
  ACC_DIURNAL_INACT \
  ACC_L5_TIME \
  ACC_M10_TIME \
  ACC_N_SLEEP_EPISODES \
  ACC_SLEEP_DUR \
  ACC_SLEEP_DUR_SD \
  ACC_SLEEP_EFF \
  ACC_SLEEP_MIDP
do

  python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
    --sumstats ${OUTDIR}/${TRAIT}.sumstats \
    --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --signed-sumstats BETA,0 \
    --p P \
    --N-col N \
    --info INFO \
    --frq AF1 \
    --out ${OUTDIR}/${TRAIT}.munged

done