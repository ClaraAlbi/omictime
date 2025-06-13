#!/bin/bash


# Generate merge list (RAP-local)
MERGE_LIST="merge_list.txt"
CMD="ls /mnt/project/data/ukbi_chr*.pgen | sed 's/\.pgen//' > ${MERGE_LIST}"

dx run swiss-army-knife \
  -icmd="${CMD}" \
  --instance-type="mem1_ssd1_v2_x16" \
  --destination="${project}:/data/" \
  --priority high \
  --brief --yes


#
# plink2 \
#   --pfile $(head -n 1 "$MERGE_LIST") \
#   --merge-list <(tail -n +2 "$MERGE_LIST") \
#   --make-bed \
#   --out "${MERGED_OUT}"
#
# plink \
#   --bfile "/mnt/project/${MERGED_OUT}" \
#   --maf 0.01 \
#   --mac 20 \
#   --geno 0.1 \
#   --mind 0.1 \
#   --hwe 1e-6 \
#   --make-bed \
#   --out "/mnt/project/${QC_OUT}"
