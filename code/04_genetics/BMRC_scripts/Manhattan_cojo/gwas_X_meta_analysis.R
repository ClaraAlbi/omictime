

##CHROM	POS	ID	REF	ALT	PROVISIONAL_REF?	A1	OMITTED	A1_FREQ	TEST	OBS_CT	BETA	SE	T_STAT	P	ERRCODE


x_Chr <- inner_join(data.table::fread("data/GWAS_res/res_chr_X_0.res_1MAF.glm.linear") %>%
  select(CHR = V1, POS = V2, BETA_F = V12, SE_F = V13, A1 = V7, A2 = V8),
  data.table::fread("data/GWAS_res/res_chr_X_1.res_1MAF.glm.linear") %>%
  select(CHR = V1, POS = V2, BETA_M = V12, SE_M = V13, A1 = V7, A2 = V8) %>% mutate(BETA_M = 2*BETA_M, SE_M = 2*SE_M),
  by = c("CHR", "POS", "A1", "A2")) 

weight_f = 1 / (x_Chr$SE_F^2)
weight_m = 1 / (x_Chr$SE_M^2)

#CHR	SNP	POS	A1	A2	N	AF1	BETA	SE	P	INFO

x_chr_meta <- x_Chr %>%
  mutate(BETA = (BETA_F * weight_f + BETA_M * weight_m) / (weight_f + weight_m),
         SE = sqrt(1 / (weight_f + weight_m)),
         P = 2 * pnorm( -abs(BETA / SE) ))

data.table::fwrite(x = x_chr_meta, file = "data/GWAS_res/res_chrX.fastGWA")
