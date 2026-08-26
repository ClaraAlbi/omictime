#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J sbayesrc_CHRONO 

#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p long
#SBATCH -t 24:00:00
#SBATCH --mem 200G
#SBATCH -c 16


ldm="/well/visscher-wray/shared/LD_ref/ukbEUR_Imputed"
gwas="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/morning_person.af.tsv"
anot="/well/visscher-wray/shared/LD_ref/annot_baseline2.2.txt"
genemap="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/gene_map_hg38_hg19.txt"

# tmp=$(mktemp --tmpdir="$(dirname "$f")" "$(basename "$f").tmp.XXXXXXXX") || exit 1
# 
# #SNP A1 A2 freq b se p N 
# 
# awk -v OFS='\t' '{
#     print $17, $18, $19, $21, $4, $5, $3, $20
#   }' "$gwas" > "$tmp"
# 
# 
# /users/visscher-wray/uus177/bin/gctb_2.5.4_Linux/gctb \
#    --impute-summary \
#    --ldm-eigen "${ldm}" \
#    --gwas-summary "${tmp}" \
#    --out "data/chrono_GWAS" \
#    --thread 16
# 
# head "$tmp"
# 
# rm -f "$tmp"

/users/visscher-wray/uus177/bin/gctb_2.5.5_Linux/gctb \
  --ldm-eigen "${ldm}" \
  --gwas-summary "data/chrono_GWAS.imputed.ma" \
  --sbayes RC \
  --annot "${anot}" \
  --out "results/chrono_GWAS_23andme" \
  --thread 16

# #Run GCTB
# /users/visscher-wray/uus177/bin/gctb_2.5.4_Linux/gctb \
#   --ldm-eigen "${ldm}" \
#   --gwas-summary "data/chrono_GWAS.imputed.ma" \
#   --sbayes RC \
#   --annot "${anot}" \
#   --num-chains 16 \
#   --out "results/chrono_GWAS_23andme" \
#   --thread 16
