#!/bin/bash

plink_step1="

plink2 \
  --pmerge-list /mnt/project/data/${MERGE_LIST} pfile \
  --set-all-var-ids @:#:\\\$r:\\\$a \
  --new-id-max-allele-len 662 \
  --make-bed \
  --out ukb_all_v3_qc
"

dx run swiss-army-knife \
  -icmd="${plink_step1}" \
  --instance-type="mem1_ssd1_v2_x72" \
  --destination="${project}:/data/" \
  --priority high \
  --brief --yes
