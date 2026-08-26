#!/bin/bash

TAR_DIR=${OMICTIME_BMRC_DIR:-/well/visscher-wray/users/uus177/gwas_acceleration}/downloads

mkdir -p logs results
# 
# for tar in "${TAR_DIR}"/*.tar; do
#   echo "Submitting ${tar}"
#   sbatch code/cojo_tar.sh "${tar}"
# done

PROTEINS=("IGFBP4" "SERPINE1" "TNFRSF12A" "HSPB6" "SEMA3F" "PLAT" "HGF" "MB" "ADAMTS15" "TNFRSF1A" "FURIN" "IL6")

for prot in "${PROTEINS[@]}"; do
  for tar in "${TAR_DIR}/${prot}"*.tar; do
    [[ -f "$tar" ]] || continue
    echo "Submitting ${tar}"
    sbatch code/cojo_tar.sh "${tar}"
  done
done