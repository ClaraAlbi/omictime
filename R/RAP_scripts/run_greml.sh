#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"


run_gcta_fast_GWAS="
    cp gcta64 \$HOME/gcta64 && chmod +x \$HOME/gcta64
    \$HOME/gcta64 --reml  --grm-adj 0  --grm-cutoff 0.05 \
                        --pheno res_v2_eur_unrel.txt \
                        --qcovar qcovar.txt \
                        --covar covar.txt \
                        --extract ukbEURu_imp_all_v3_impQC_maf01.snpList \
                        --grm /mnt/project/grm/ukbPPP_imp_v3_impQC \
                        --out reml_res \
                        --thread-num 8
    "

dx run swiss-army-knife \
      -iin="gcta64" \
      -iin="${project}:/res_v2_eur_unrel.txt" \
      -iin="${project}:/covar.txt" \
      -iin="${project}:/qcovar.txt" \
      -iin="${project}:/grm/ukbEURu_imp_all_v3_impQC_maf01.snpList" \
      -icmd="${run_gcta_fast_GWAS}" \
      --priority high \
      --cost-limit 10 \
      --instance-type mem1_ssd1_v2_x72 \
      --destination="${project}:/circadian/GWAS/" \
      --tag res_chr${chr} \
      -y \
      --brief
