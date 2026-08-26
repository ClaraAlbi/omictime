#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J sbayesrc_res 

#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p long
#SBATCH -t 34:00:00
#SBATCH --mem 200G
#SBATCH -c 32


ldm="/well/visscher-wray/shared/LD_ref/ukbEUR_Imputed"
gwas="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_NOREP/res_v2_eur_norep.fastGWA"
anot="/well/visscher-wray/shared/LD_ref/annot_baseline2.2.txt"
# 
# tmp=$(mktemp --tmpdir="$(dirname "$f")" "$(basename "$f").tmp.XXXXXXXX") || exit 1
# 
# awk -v OFS='\t' '{
#     print $2, $4, $5, $7, $8, $9, $10, $6
#   }' "$gwas" > "$tmp"
# 
# 
# /users/visscher-wray/uus177/bin/gctb_2.5.4_Linux/gctb \
#    --impute-summary \
#    --ldm-eigen "${ldm}" \
#    --gwas-summary "${tmp}" \
#    --out "data/res_v2_eur_norep" \
#    --thread 16
# 
# rm -f "$tmp"

# Run GCTB
/users/visscher-wray/uus177/bin/gctb_2.5.4_Linux/gctb \
  --ldm-eigen "${ldm}" \
  --gwas-summary "data/res_v2_eur_norep.imputed.ma" \
  --sbayes RC \
  --annot "${anot}" \
#  --num-chains 32 \
  --out "results/sbayesrc/res_v2_eur_norep_v2" \
  --thread 32
