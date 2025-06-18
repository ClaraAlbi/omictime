#!/bin/bash

# Generate merge list (RAP-local)

plink_step1="

plink2 \
  --pmerge-list /mnt/project/test_merge.txt pfile \
  --set-all-var-ids @:#:\\\$r:\\\$a \
  --new-id-max-allele-len 662 \
  --make-bed \
  --out ukbi_test
"

dx run swiss-army-knife \
  -icmd="${plink_step1}" \
  --instance-type="mem1_ssd1_v2_x36" \
  --destination="${project}:/data/" \
  --priority high \
  --brief --yes
