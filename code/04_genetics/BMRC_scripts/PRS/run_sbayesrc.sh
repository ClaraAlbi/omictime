#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J sbayesrc_res 

#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p long
#SBATCH -t 24:00:00
#SBATCH --mem 200G
#SBATCH -c 48


ldm="/well/visscher-wray/shared/LD_ref/ukbEUR_Imputed"
gwas="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v4/res_v4_GWAS.fastGWA"
anot="/well/visscher-wray/shared/LD_ref/annot_baseline2.2.txt"
genemap="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/gene_map_hg38_hg19.txt"

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
#    --out "data/res_GWAS" \
#    --thread 16
# 
# rm -f "$tmp"

# Run GCTB
# /users/visscher-wray/uus177/bin/gctb_2.5.4_Linux/gctb \
#   --ldm-eigen "${ldm}" \
#   --gwas-summary "data/res_GWAS.imputed.ma" \
# #  --sbayes RC \
#   --gwfm RC \
#   --gene-map "${genemap}" \
#   --annot "${anot}" \
#   --num-chains 16 \
#   --out "results/res_GWAS_v4_RC" \
#   --thread 16

/users/visscher-wray/uus177/bin/gctb_2.5.5_Linux/gctb \
  --ldm-eigen "${ldm}" \
  --gwas-summary "data/res_GWAS.imputed.ma" \
  --gwfm RC \
  --gene-map "${genemap}" \
  --annot "${anot}" \
  --out "results/res_GWAS_v4_RC" \
  --thread 48
