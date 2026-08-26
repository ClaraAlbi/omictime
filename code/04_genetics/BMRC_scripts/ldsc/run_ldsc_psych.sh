#!/bin/bash 

# #SBATCH -A visscher-wray.prj
# #SBATCH -J cojo
# 
# #SBATCH -o logs/my-job-%j.out
# #SBATCH -e logs/my-job-%j.err
# #SBATCH -p short
# #SBATCH -t 12:00:00
# #SBATCH --mem 50G
# #SBATCH -c 16

# Psychiatric traits schizophrenia (rG = −0.11, P = 1 × 10−7), 
# depressive symptoms (rG = −0.16; P = 2 × 10−6), 
# major depressive disorder (rG = −0.19; P = 3 × 10−5) and 
# intelligence (rG = −0.11; P = 8 × 10−6) were all negatively correlated with the morning chronotype.

#name="pgc-mdd2024_Clin_eur"
#name="pgc-mdd2024_EHR_eur"
# name="pgc-mdd2024_Quest_eur"
# 
# gwas="/well/visscher-wray/shared/PGC_MDD/${name}_v3-49-24-11.tsv.gz"
# 
# tmpfile=$(mktemp)
# gzip -dc "$gwas" | awk '/^#CHROM/ {start=1} start'   > "$tmpfile"
# 
# #head -n 1 "$tmpfile" | tr '\t' '\n' | nl -ba
# #head "$tmpfile"
# 
# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats "$tmpfile" \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp ID \
#     --a1 EA \
#     --a2 NEA \
#     --signed-sumstats BETA,0 \
#     --p PVAL \
#     --N-col NEFF \
#     --info IMPINFO \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_MDD/${name}
# 
# rm "$tmpfile"

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats /exafs1/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Binary_Traits/SCZ_02/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.ma \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats b,0 \
#     --p p \
#     --frq freq \
#     --N-col N \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_MDD/SCZ

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats /exafs1/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Binary_Traits/BIP_02/pgc-bip2021-all.vcf.tsv_tidy.ma \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats b,0 \
#     --p p \
#     --frq freq \
#     --N-col N \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_MDD/BIP
