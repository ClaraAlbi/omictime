#!/bin/bash 

ZIP=data/acc/accel_GWAS_all_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.zip

for pheno in ACC_DIURNAL_INACT ACC_L5_TIME ACC_M10_TIME ACC_N_SLEEP_EPISODES \
             ACC_SLEEP_DUR ACC_SLEEP_DUR_SD ACC_SLEEP_EFF ACC_SLEEP_MIDP; do

  case $pheno in
    ACC_DIURNAL_INACT) N=84757 ;;
    ACC_L5_TIME) N=85205 ;;
    ACC_M10_TIME) N=85670 ;;
    ACC_N_SLEEP_EPISODES) N=84810 ;;
    ACC_SLEEP_DUR) N=85449 ;;
    ACC_SLEEP_DUR_SD) N=84441 ;;
    ACC_SLEEP_EFF) N=84810 ;;
    ACC_SLEEP_MIDP) N=84810 ;;
  esac

  unzip -p $ZIP \
  | awk -v beta="${pheno}_RAW_SIN_BETA" \
        -v se="${pheno}_RAW_SIN_SE" \
        -v pval="${pheno}_RAW_SIN_P" \
        -v af="${pheno}_RAW_SIN_A1FREQ" \
        -v N=$N '
  BEGIN { OFS="\t" }
  NR==1 {
      for (i=1;i<=NF;i++) idx[$i]=i
      print "SNP","CHR","BP","A1","A2","BETA","SE","P","INFO","AF1","N"
      next
  }
  {
      b = $idx[beta]
      s = $idx[se]
      p = $idx[pval]
      f = $idx[af]
      info = $idx["INFO"]

      # Drop missing or invalid numeric values
      if (b=="NA" || b=="-" || p=="NA" || p=="-" ||
          f=="NA" || f=="-" || info=="NA" || info=="-") next

      print $idx["SNP"], \
            $idx["CHR"], \
            $idx["BP"], \
            $idx["ALLELE1"], \
            $idx["ALLELE0"], \
            b, \
            s, \
            p, \
            info, \
            f, \
            N
  }' \
  | gzip -c \
  > data/acc/${pheno}.sumstats.gz
done
