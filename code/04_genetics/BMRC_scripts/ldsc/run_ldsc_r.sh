#!/bin/bash

#SBATCH -A visscher-wray.prj
#SBATCH -J cojo

#SBATCH -o logs/my-job-%j.out
#SBATCH -e logs/my-job-%j.err
#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH --mem 50G
#SBATCH -c 16

TARGET="/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Binary_Traits"

# Paths
MUNGE=/users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py
ALLELES=/well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist
OUTDIR=${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_other/

# Find all tidy.ma files recursively
find "$TARGET" -type f -name "*tidy.ma*" | while read FILE
do
    echo "Processing: $FILE"

    # Extract base filename without extension(s)
    BASENAME=$(basename "$FILE")
    BASENAME=${BASENAME%_tidy.ma}

    # Create output file prefix
    OUTFILE="${OUTDIR}/${BASENAME}"

    # Run munge_sumstats
    python $MUNGE \
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
