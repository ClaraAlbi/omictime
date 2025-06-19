#!/bin/bash


# Generate merge list (RAP-local)
MERGE_LIST="merge_list.txt"
CMD="ls /mnt/project/data/ukbi_chr*_qc.pgen | sed 's/\.pgen//' > ${MERGE_LIST}"
#
# dx run swiss-army-knife \
#   -icmd="${CMD}" \
#   --instance-type="mem1_ssd1_v2_x16" \
#   --destination="${project}:/data/" \
#   --priority high \
#   --brief --yes


plink_step1="

plink2 \
  --pmerge-list /mnt/project/data/${MERGE_LIST} pfile \
  --set-all-var-ids @:#:\\\$r:\\\$a \
  --new-id-max-allele-len 662 \
  --make-bed \
  --out ukbi_v3_qc
"

dx run swiss-army-knife \
  -icmd="${plink_step1}" \
  --instance-type="mem1_ssd1_v2_x72" \
  --destination="${project}:/data/" \
  --priority high \
  --brief --yes

plink_step2="

plink2 \
  --bfile /mnt/project/data/ukbi_v3_qc \
  --geno 0.1 \
  --mind 0.1 \
  --make-bed \
  --out ukbi_v3_qc2
"
#
# dx run swiss-army-knife \
#   -icmd="${plink_step2}" \
#   --instance-type="mem1_ssd1_v2_x36" \
#   --destination="${project}:/data/" \
#   --priority high \
#   --brief --yes
#
