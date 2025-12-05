#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"

  run_gcta_fast_GWAS="
    cp gcta64 \$HOME/gcta64 && chmod +x \$HOME/gcta64
    \$HOME/gcta64 --mbfile geno_chrs.txt \
                        --pheno phenotypes_rel.txt \
                        --qcovar qcovar.txt \
                        --covar covar.txt \
                        --fastGWA-mlm \
                        --extract ukbEURu_imp_all_v3_impQC_maf01.snpList \
                        --grm-sparse /mnt/project/grm/sp_grm_eur_OX \
                        --covar-maxlevel 110 \
                        --out res_rel_chr${chr} \
                        --thread-num 8
    "

  dx run swiss-army-knife \
      -iin="gcta64" \
      -iin="${project}:/phenotypes_rel.txt" \
      -iin="${project}:/geno_chrs.txt" \
      -iin="${project}:/covar.txt" \
      -iin="${project}:/qcovar_prots.txt" \
      -iin="${project}:/grm/ukbEURu_imp_all_v3_impQC_maf01.snpList" \
      -icmd="${run_gcta_fast_GWAS}" \
      --priority high \
      --cost-limit 10 \
      --instance-type mem3_ssd1_v2_x8 \
      --destination="${project}:/circadian/GWAS/" \
      --tag res_chr${chr} \
      -y \
      --brief

