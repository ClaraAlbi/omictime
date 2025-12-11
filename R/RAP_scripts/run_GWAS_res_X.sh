#!/bin/bash


project="blood_biomarkers - Jul 01, 2024"

COVARS="age_recruitment,batch,PC1-PC20"

run_gcta_fast_GWAS="
    plink2 --bgen \"/mnt/project/Bulk/Imputation/Imputation from genotype (TOPmed)/ukb21007_cX_b0_v1.bgen\" ref-first \
                        --sample \"/mnt/project/Bulk/Imputation/Imputation from genotype (TOPmed)/ukb21007_cX_b0_v1.sample\" \
                        --pheno-name res \
                        --pheno phenotypes_1.txt \
                        --covar-name ${COVARS} \
                        --covar covariates.txt  \
                        --no-input-missing-phenotype \
                        --glm hide-covar \
                        --covar-variance-standardize \
                        --out res_chr_X_1 \
                        --thread-num 8
    "

  dx run swiss-army-knife \
      -iin="${project}:/phenotypes_1.txt" \
      -iin="${project}:/covariates.txt" \
      -icmd="${run_gcta_fast_GWAS}" \
      --priority high \
      --cost-limit 10 \
      --instance-type mem3_ssd1_v2_x8 \
      --destination="${project}:/circadian/GWAS/" \
      --tag res_chr${chr} \
      -y \
      --brief

