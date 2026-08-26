#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J cojo

#SBATCH -o logs/my-job-%j.out
#SBATCH -e logs/my-job-%j.err
#SBATCH -p short
#SBATCH -t 12:00:00
#SBATCH --mem 50G
#SBATCH -c 16

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/CA_v5.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v5/res_GWAS.munged

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v4/res_v4_GWAS.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v4/res_GWAS.munged


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_b.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_b.munged

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/logistic/RAgroup_noFA_noFU.assoc.logistic.gz \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --frq MAF \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N 71500 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/RAgroup.logistic.munge

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/linear/ra_full_noFA.assoc.linear.gz \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --frq MAF \
#     --signed-sumstats BETA,0 \
#     --N 77440 \
#     --p P \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/RA/RAgroup.linear.munge



# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_MDD/bip2024_eur_no23andMe \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --N-col Neff_half \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_MDD/bip2024.munged


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_v2_GWAS.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_v2_GWAS.munged
# 


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/arth/GCST90132223_buildGRCh37.tsv.gz \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp variant_id \
#     --a1 effect_allele \
#     --a2 other_allele \
#     --signed-sumstats beta,0 \
#     --p p_value \
#     --N 97173 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/arth/GCST90132223_buildGRCh37.munged

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats /well/visscher-wray/users/uus177/psych_sleep/data/GCST90399571.h.tsv  \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp rsid \
#     --a1 effect_allele \
#     --a2 other_allele \
#     --signed-sumstats beta,0 \
#     --p p_value \
#     --N 27200 \
#     --out /well/visscher-wray/users/uus177/psych_sleep/data/RL.munged

python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
    --sumstats /well/visscher-wray/users/uus177/psych_sleep/data/RL_23andme.txt  \
    --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
    --snp assay.name \
    --a1 A1 \
    --a2 A0 \
    --signed-sumstats effect,0 \
    --p pvalue \
    --N 65640 \
    --out /well/visscher-wray/users/uus177/psych_sleep/data/RL_23andme.munged


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res/res_GWAS_merged.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/res_GWAS.munged
# 


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_lasso/lasso_GWAS.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_lasso/lasso_GWAS.munged

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v3/res_v3_GWAS.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v3/res_v3_GWAS.munged


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_prots/res_prots_GWAS.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_prots/res_prots_GWAS.munged
# 

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_pPCA_v2/res_pPCA_v2_GWAS.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_pPCA_v2/res_pPCA_v2_GWAS.munged


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_pca/res_GWAS_pPCA.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_pca/res_GWAS_pPCA


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_cond/res_GWAS_cond.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_cond/res_GWAS_cond_munged.fastGWA


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_UKB_GWAS_merged.fastGWA \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats BETA,0 \
#     --p P \
#     --N-col N \
#     --info INFO \
#     --frq AF1 \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_UKB_GWAS_merged_munged.fastGWA

    
# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/morning_person-5.1.cojo \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats b,0 \
#     --p P \
#     --N-col N \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/morning_person-5.1.munged

# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/23andme/insomnia-5.2.with_rs_a_n.dat \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp assay.name \
#     --a1 A2 \
#     --a2 A1 \
#     --signed-sumstats effect,0 \
#     --p pvalue \
#     --N-col im.num.sum \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/23andme/Jansen_2018_insomnia_broad.munged


# python /users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py \
#     --sumstats /exafs1/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Binary_Traits/CNT_02/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt_tidy.ma.full \
#     --merge-alleles /well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist \
#     --snp SNP \
#     --a1 A1 \
#     --a2 A2 \
#     --signed-sumstats b,0 \
#     --p p \
#     --frq freq \
#     --N-col N \
#     --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_UKB

# 
# 
# 
# INPUT_DIR="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/selected"
# OUT_DIR="${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/selected"
# LDSC="/users/visscher-wray/uus177/bin/ldsc/munge_sumstats.py"
# HM3="/well/visscher-wray/shared/eur_w_ld_chr/w_hm3.snplist"
# 
# for f in "$INPUT_DIR/BP-ICE_EUR_HTN_15-04-2020.txt"
# do
#   echo "Processing $f"
# 
#   TMPFILE=$(mktemp)
# 
#   # Remove MarkerName and weight columns
#   awk -F'\t' 'BEGIN{OFS="\t"}
#   NR==1{
#     for(i=1;i<=NF;i++){
#       col=tolower($i)
#       if(col!="markername" && col!="weight")
#         keep[++k]=i
#     }
#     for(j=1;j<=k;j++) printf "%s%s", $(keep[j]), (j<k?OFS:ORS)
#     next
#   }
#   {
#     for(j=1;j<=k;j++) printf "%s%s", $(keep[j]), (j<k?OFS:ORS)
#   }' "$f" > "$TMPFILE"
# 
#   OUT=$(basename "$f" .txt)
# 
#   # Run LDSC munge
#   python "$LDSC" \
#     --sumstats "$TMPFILE" \
#     --merge-alleles "$HM3" \
#     --a1 Alelle1 \
#     --a2 Allele2 \
#     --signed-sumstats Zscore,0 \
#     --p "P-value" \
#     --N-col N \
#     --frq Freq1 \
#     --out "$OUT_DIR/${OUT}.munge"
# 
#   rm "$TMPFILE"
# done
# 
# 

# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_UKB.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res/res_GWAS.munged.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/chrono_vs_res_rg

# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_2018_morning_person.munged.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_cond/res_GWAS_cond_munged.fastGWA.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/morning_vs_res_cond_rg
# 
# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_UKB_GWAS_merged_munged.fastGWA.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res/res_GWAS.munged.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/chronoq_vs_res_rg

# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_UKB_GWAS_merged_munged.fastGWA.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_cond/res_GWAS_cond_munged.fastGWA.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/chronoq_vs_res_cond_rg
# 
# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_UKB.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res/res_GWAS.munged.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/jones_vs_res_rg
# 
# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_UKB.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_cond/res_GWAS_cond_munged.fastGWA.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/jones_vs_res_cond_rg
# 




# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_UKB.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v3/res_v3_GWAS.munged.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/jones_vs_res_v3_rg
# 
# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/chrono_UKB_GWAS_merged_munged.fastGWA.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v3/res_v3_GWAS.munged.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/chronoq_vs_res_v3_rg
# 
# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_2018_morning_person.munged.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/data/GWAS_res_v3/res_v3_GWAS.munged.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/morning_vs_res_v3_rg
# 
# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v2/res_v2_GWAS.munged.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res_v3/res_v3_GWAS.munged.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/res_v2_vs_res_v3_rg




# python /users/visscher-wray/uus177/bin/ldsc/ldsc.py \
#   --rg ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/chronotype/Jones_2018_morning_person.munged.sumstats.gz,${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/data/GWAS_res/res_GWAS.munged.sumstats.gz \
#   --ref-ld-chr /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --w-ld-chr  /well/visscher-wray/shared/eur_w_ld_chr/ \
#   --out ${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/results/morning_vs_res_rg

