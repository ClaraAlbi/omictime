#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"
imp_dir="/Bulk/Imputation/UKB imputation from genotype/"

for CHR in {1..22}; do
  echo "Processing chromosome $CHR..."

  BGEN_FILE="ukb22828_c${CHR}_b0_v3.bgen"
  SAMPLE_FILE="ukb22828_c${CHR}_b0_v3.sample"

  plink_command="plink2 \
    --bgen ${BGEN_FILE} ref-first \
    --sample ${SAMPLE_FILE} \
    --maf 0.01 \
    --hwe 1e-6 \
    --make-pgen \
    --out ukb_all_chr${CHR}_v3_qc"

  dx run swiss-army-knife \
    -iin="${imp_dir}${BGEN_FILE}" \
    -iin="${imp_dir}${SAMPLE_FILE}" \
    -icmd="${plink_command}" \
    --instance-type="mem1_ssd1_v2_x36" \
    --destination="${project}:/data/" \
    --priority high \
    --brief --yes

done
