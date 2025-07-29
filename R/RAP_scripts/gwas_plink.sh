#!/bin/bash

project="blood_biomarkers - Jul 01, 2024"
phenotypes="phenotypes.txt"
covariates="covariates_cojo.txt"

COVARS="sex,age_recruitment,batch,PC1-PC20,1:151926201:G:A_G,1:155049402:T:C_T,1:155106054:A:G_A,1:172389198:A:G_A,1:178529898:G:A_G,1:178650869:G:A_G,1:178716163:G:A_G,1:178862779:A:G_A,1:178868454:C:A_C,1:178907461:C:T_C,4:1114909:A:G_A,4:1162694:T:G_T,4:77371434:G:C_G,4:103656377:A:G_A,5:147386552:C:T_C,5:147390110:T:C_T,5:147466342:A:C_A,5:147499250:G:A_G,5:147500883:C:T_C,5:147606425:G:A_G,5:180236824:C:CCA_C,11:32308936:G:T_G,11:61605499:T:A_T,11:73091913:G:A_G,11:73169662:C:CT_C,11:126250680:G:A_G,17:7074125:G:A_G,19:18482358:G:A_G,19:18498246:G:A_G,19:18499858:T:C_T,19:49206462:C:T_C,19:55999142:C:G_C"

for pheno in res; do
  echo "==> Running GWAS for phenotype: $pheno"

  plink_command="plink2 \
    --pfile /mnt/project/data/ukbi_v3_qc \
    --pheno-name ${pheno} \
    --pheno $phenotypes  \
    --covar-name ${COVARS} \
    --covar ${covariates}  \
    --no-input-missing-phenotype \
    --glm hide-covar \
    --covar-variance-standardize \
    --threads 16   \
    --out gwas_${pheno}_pQTL"

  dx run swiss-army-knife \
    -iin="${phenotypes}" \
    -iin="${covariates}"\
    -icmd="${plink_command}" \
    --instance-type="mem1_ssd1_v2_x36" \
    --destination="${project}:/gwas/" \
    --priority high \
    --brief --yes

done
