#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"
phenotypes="phenotypes.txt"

COVARS="sex,age_recruitment,batch,PC1-PC20"

for pheno in res res_abs; do
  echo "==> Running GWAS for phenotype: $pheno"

  plink_command="plink2 \
    --pfile /mnt/project/data/ukbi_v3_qc \
    --pheno-name ${pheno} \
    --pheno $phenotypes  \
    --covar-name ${COVARS} \
    --covar $phenotypes  \
    --no-input-missing-phenotype \
    --glm hide-covar \
    --covar-variance-standardize \
    --threads 16   \
    --out gwas_${pheno}"

  dx run swiss-army-knife \
    -iin="${phenotypes}" \
    -icmd="${plink_command}" \
    --instance-type="mem1_ssd1_v2_x36" \
    --destination="${project}:/gwas/" \
    --priority high \
    --brief --yes

done
