#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"
snpfile="cojo_pqtls.txt"

plink_command="plink2 \
    --pfile /mnt/project/data/ukbi_v3_qc \
    --extract ${snpfile} \
    --recodeA \
    --out subset_cojo_pqtls"

dx run swiss-army-knife \
  -iin="${snpfile}" \
  -icmd="${plink_command}" \
  --instance-type="mem1_ssd1_v2_x36" \
  --destination="${project}:/snps/" \
  --priority high \
  --brief --yes
