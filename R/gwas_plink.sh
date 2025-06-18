#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"
phenotypes="phenotypes.txt"
covariates="covariates.txt"

for CHR in 22; do
  echo "Processing chromosome $CHR..."

  BGEN_FILE="ukb22828_c${CHR}_b0_v3.bgen"
  SAMPLE_FILE="ukb22828_c${CHR}_b0_v3.sample"

  plink_command="plink2 \
    --pfile /mnt/project/data/ukbi_chr${CHR}_v3_qc \
    --pheno-name gap \
    --pheno $phenotypes  \
    --covar-name sex,age_recruitment,batch,PC1-PC10 \
    --covar $covariates  \
    --no-input-missing-phenotype \
    --glm firth-fallback  \
    --threads 8   \
    --out gwas_gap_chr${CHR}"

  dx run swiss-army-knife \
    -iin="${phenotypes}" \
    -iin="${covariates}" \
    -icmd="${plink_command}" \
    --instance-type="mem1_ssd1_v2_x16" \
    --destination="${project}:/gwas/" \
    --priority high \
    --brief --yes

done
