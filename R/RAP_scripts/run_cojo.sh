#!/bin/bash


project="blood_biomarkers - Jul 01, 2024"

#gcta64  --bfile test  --chr 1 --maf 0.01 --cojo-file test.ma --cojo-slct --out test_chr1


for chr in {19..19}; do
  run_cojo_GWAS="
    cp gcta64 \$HOME/gcta64 && chmod +x \$HOME/gcta64
    \$HOME/gcta64 --bgen \"/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.bgen\" \
                        --sample \"/mnt/project/Bulk/Imputation/UKB imputation from genotype/ukb22828_c${chr}_b0_v3.sample\" \
                        --chr ${chr} \
                        --cojo-file res_chr${chr}.fastGWA \
                        --cojo-slct \
                        --out cojo_res_chr${chr} \
                        --thread-num 8
    "

  dx run swiss-army-knife \
      -iin="gcta64" \
      -iin="${project}:/circadian/GWAS/res_chr${chr}.fastGWA" \
      -icmd="${run_cojo_GWAS}" \
      --priority high \
      --cost-limit 10 \
      --instance-type mem3_ssd1_v2_x8 \
      --destination="${project}:/circadian/GWAS/" \
      --tag res_chr${chr} \
      -y \
      --brief
done
