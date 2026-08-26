#!/bin/bash 

#SBATCH -A visscher-wray.prj
#SBATCH -J sbayesrc_res 

#SBATCH -o logs/my-job-%j.out 
#SBATCH -e logs/my-job-%j.err 
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH --mem 200G
#SBATCH -c 16


ldm="/well/visscher-wray/shared/LD_ref/ukbEUR_Imputed"

/users/visscher-wray/uus177/bin/gctb_2.5.5_Linux/gctb \
  --ldm-eigen "${ldm}" --get-pwld --rsq 0.5 --out ${ldm}/rsq0.5 --thread 12 > ${ldm}/rsq0.5.log 2>&1