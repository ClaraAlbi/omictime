#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"
snpfile="set_cojo.txt"


for chr in {19..19}; do
  plink_command="plink2 \
      --bgen \"/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.bgen\" \
      --ref-first \
      --sample \"/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.sample\" \
      --extract ${snpfile} \
      --recode A \
      --out subset_cojo_pQTL"

  dx run swiss-army-knife \
    -iin="${snpfile}" \
    -icmd="${plink_command}" \
    --instance-type="mem1_ssd1_v2_x36" \
    --destination="${project}:/snps/" \
    --priority high \
    --brief --yes
done
