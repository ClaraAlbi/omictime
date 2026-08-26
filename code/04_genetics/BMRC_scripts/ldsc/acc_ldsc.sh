#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J cojo

#SBATCH -o logs/my-job-%j.out
#SBATCH -e logs/my-job-%j.err
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 50G
#SBATCH -c 16

# Directories
OUTDIR="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/acc/"
RESDIR="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/v2/"
mkdir -p "$RESDIR"

# Collect all LDSC-munged sumstats files in OUTDIR

#chrono=(${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_v2_GWAS.munged.sumstats.gz ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v4/res_GWAS.munged.sumstats.gz )
chrono=(${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_2018_morning_person.munged.sumstats.gz )
psych=( $(find "$OUTDIR" -type f -name "*.sumstats.gz") )

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
      --rg "${c},${p}" \
      --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
      --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
      --out "${RESDIR}/${out}"
  
  done
done
