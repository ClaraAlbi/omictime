#!/bin/bash

#SBATCH -A visscher-wray.prj
#SBATCH -J cojo

#SBATCH -o logs/my-job-%j.out
#SBATCH -e logs/my-job-%j.err
#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH --mem 50G
#SBATCH -c 16


# Paths
MUNGE=/users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py
ALLELES=/well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist
OUTDIR=${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_other/

FILES=(
  "/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Quantitative_Traits/mtCNT_01/mtDNA_CN.ALLm2.bgen.stats_imp.ma" 
  "/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Quantitative_Traits/BMI_02/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt_imp.ma" 
  "/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Quantitative_Traits/FastingGlucose_02/MAGIC_FastingGlucose.txt_imp.ma" 
  "/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Quantitative_Traits/Insomnia_01/1200_imp.ma" 
  "/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Quantitative_Traits/HipCir_01/GCST012229.tsv_imp.ma"
)


for FILE in "${FILES[@]}"
do
    echo "Processing: $FILE"

    BASENAME=$(basename "$FILE")
    OUTFILE="${OUTDIR}/${BASENAME}"

    python "$MUNGE" \
        --sumstats "$FILE" \
        --merge-alleles "$ALLELES" \
        --snp SNP \
        --a1 A1 \
        --a2 A2 \
        --signed-sumstats b,0 \
        --p p \
        --frq freq \
        --N-col N \
        --out "$OUTFILE"
done
