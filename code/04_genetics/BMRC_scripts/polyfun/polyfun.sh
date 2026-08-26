#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J polyfun

#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem 16G
#SBATCH -c 16
##SBATCH --array=1-48

# python ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/polyfun/polyfun/munge_polyfun_sumstats.py \
#   --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/CA_v5.fastGWA \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/res_v5_GWAS.parquet
# 
# python ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/polyfun/polyfun/polyfun.py \
#     --compute-h2-L2 \
#     --no-partitions \
#     --output-prefix results/polyfun/res_v5_GWAS \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/res_v5_GWAS.parquet \
#     --ref-ld-chr /well/visscher-wray/shared/annotations/baselineLF2.2.UKB/baselineLF2.2.UKB. \
#     --w-ld-chr /well/visscher-wray/shared/annotations/baselineLF2.2.UKB/weights.UKB.


# python ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/polyfun/polyfun/create_finemapper_jobs.py \
#     --sumstats results/polyfun/res_v5_GWAS.${chr}.snpvar_ridge_constrained.gz \
#     --n 44547 \
#     --method finemap \
#     --finemap-exe /exafs1${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/polyfun/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 \
#     --pvalue-cutoff 5e-8 \
#     --max-num-causal 5 \
#     --out-prefix results/polyfun/polyfun_chr${chr} \
#     --jobs-file results/polyfun/polyfun_all_chr${chr}.txt

# for chr in {1..22}; do
#   python ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/polyfun/polyfun/create_finemapper_jobs.py \
#     --sumstats results/polyfun/res_v5_GWAS.${chr}.snpvar_ridge_constrained.gz \
#     --n 44547 \
#     --method finemap \
#     --finemap-exe /exafs1${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/polyfun/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64 \
#     --pvalue-cutoff 5e-8 \
#     --max-num-causal 5 \
#     --out-prefix results/polyfun/polyfun_v2 \
#     --jobs-file results/polyfun/polyfun_all_chr${chr}.txt
# done

#cat results/polyfun/polyfun_all_chr*.txt > results/polyfun/polyfun_all_jobs_v2.txt

#JOBS_FILE="results/polyfun/polyfun_all_jobs_v2.txt"

#CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOBS_FILE}")
#bash -lc "${CMD}"


# 
 python ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/polyfun/polyfun/aggregate_finemapper_results.py \
     --out-prefix results/polyfun/polyfun_v2 \
     --pvalue-cutoff 5e-8 \
     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/res_v5_GWAS.parquet \
     --out results/polyfun/polyfun_CA.txt.gz


