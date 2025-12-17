#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"

for chr in {1..22}; do
  run_gcta_fast_GWAS="
    cp gcta64 \$HOME/gcta64 && chmod +x \$HOME/gcta64
    \$HOME/gcta64 --bgen \"/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.bgen\" \
                        --sample \"/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.sample\" \
                        --pheno phenotypes_prots.txt \
                        --mpheno 2 \
                        --qcovar qcovar.txt \
                        --covar covar.txt \
                        --fastGWA-mlm \
                        --extract ukbEURu_imp_all_v3_impQC_maf01.snpList \
                        --grm-sparse /mnt/project/grm/sp_grm_eur_OX \
                        --covar-maxlevel 110 \
                        --out res_prots_chr${chr} \
                        --thread-num 8
    "

  dx run swiss-army-knife \
      -iin="gcta64" \
      -iin="${project}:/phenotypes_prots.txt" \
      -iin="${project}:/covar.txt" \
      -iin="${project}:/qcovar.txt" \
      -iin="${project}:/grm/ukbEURu_imp_all_v3_impQC_maf01.snpList" \
      -icmd="${run_gcta_fast_GWAS}" \
      --priority high \
      --cost-limit 10 \
      --instance-type mem3_ssd1_v2_x8 \
      --destination="${project}:/circadian/GWAS/" \
      --tag res_chr${chr} \
      -y \
      --brief
done
