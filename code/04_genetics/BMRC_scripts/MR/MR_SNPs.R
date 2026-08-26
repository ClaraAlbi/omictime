source("config/paths.R")

library(TwoSampleMR)
library(tidyr)
library(dplyr)
library(stringr)
library(ggrepel)
library(purrr)
### COJO

files <- list.files(bmrc_path("results/cojo/"), "jma", full.names = T)

cojo_CA <- bind_rows(lapply(files[str_detect(files, "cojo_res_v4")], data.table::fread)) %>%
  left_join(data.table::fread("annovar/CAannotate.variant_function") %>% select(Chr = V3, annotation = V1, genes = V2, bp = V4) %>% mutate(Chr = as.numeric(str_remove(Chr, "chr")))) %>%
  arrange(Chr)

gwas_CA %>%
  filter(P_CA < 5e-8) %>%
  group_by(CHR) %>% count()

#             
#             %>% group_by(V3, V1) %>% summarise(g = list(unique(V12))), by = c("bp" = "V3")) %>%
#   mutate(g = map_chr(g, ~paste0(.x, collapse = ",")))

cojo_CA %>% select(SNP, Chr, bp) %>% data.table::fwrite("cojo_snps.txt", row.names = F, sep = " ")

# hr1    1234567    1234567    A    G
# chr2    7654321    7654321    C    T
cojo_CA %>%
  left_join(gwas_CA %>% select(SNP, A1, A2)) %>%
  mutate(chr = paste0("chr", Chr), bp2 = bp) %>%
  select(chr, bp, bp2, A1, A2) %>%
  data.table::fwrite(., "cojo_CA.avinput", col.names = F, sep = "\t")

write(cojo_CA$SNP, "cojo_CA.txt")




cojo_chrono <- bind_rows(lapply(files[str_detect(files, "chrono_v2")], data.table::fread)) %>%
  left_join(data.table::fread("cojo/chronoannotate.variant_function") %>% select(Chr = V3, annotation = V1, genes = V2, bp = V4) %>% mutate(Chr = as.numeric(str_remove(Chr, "chr"))) %>%
               distinct(bp, .keep_all = T)) %>%
  arrange(Chr)

write(cojo_chrono$SNP, "cojo/cojo_chrono.txt")

cojo_chrono %>%
  left_join(gwas_chrono %>% select(SNP, A1, A2)) %>%
  mutate(chr = paste0("chr", Chr), bp2 = bp) %>%
  select(chr, bp, bp2, A1, A2) %>%
  data.table::fwrite(., "cojo/cojo_chrono.avinput", col.names = F, sep = "\t")


###

gwas_chrono <- data.table::fread("data/chronotype/chrono_v2_GWAS.fastGWA")

gwas_CA <- data.table::fread("data/GWAS_res_v5/CA_v5.fastGWA")

gwas_M10 <-  

#gwas_SCZ <- data.table::fread("/well/visscher-wray/shared/UKBiobank_info_from_UQ/SBayesRC/Uno_Traits/Binary_Traits/SCZ_02/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv_imp.ma")



exposure_dat <- TwoSampleMR::format_data(
  cojo_CA %>% 
    filter(!is.na(genes)) %>%
    as.data.frame() %>%
    mutate(
      Phenotype = "CA"),
  type = "exposure",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "refA",
  eaf_col = "freq",
  pval_col = "p",
  samplesize_col = "n"
)

outcome_dat <- TwoSampleMR::format_data(
  gwas_chrono %>% 
    filter(SNP %in% cojo_CA$SNP) %>%
    as_tibble() %>% 
    mutate(Phenotype = "Chronotype"),
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1",
  pval_col = "P",
  samplesize_col = "N"
)


harm <- TwoSampleMR::harmonise_data(
  exposure_dat,
  outcome_dat,
  action = 2
)

mr_res <- TwoSampleMR::mr(
  harm,
  method_list = c(
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)




#### CHRONO -> CA

exposure_dat_chrono <- TwoSampleMR::format_data(
  cojo_chrono %>% 
    filter(!is.na(genes)) %>%
    as.data.frame() %>%
    mutate(
      Phenotype = "CA"),
  type = "exposure",
  snp_col = "SNP",
  beta_col = "b",
  se_col = "se",
  effect_allele_col = "refA",
  eaf_col = "freq",
  pval_col = "p",
  samplesize_col = "n"
)

outcome_dat_CA <- TwoSampleMR::format_data(
  gwas_CA %>% 
    filter(SNP %in% cojo_chrono$SNP) %>%
    as_tibble() %>% 
    mutate(Phenotype = "Chronotype"),
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1",
  pval_col = "P",
  samplesize_col = "N"
)


harm_2 <- TwoSampleMR::harmonise_data(
  exposure_dat_chrono,
  outcome_dat_CA,
  action = 2
) %>%left_join(cojo_chrono, by = "SNP") %>%
  mutate(genes = str_remove_all(genes, "\\(.*?\\).*")) 

mr_res_2 <- TwoSampleMR::mr(
  harm_2,
  method_list = c(
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)

het    <- TwoSampleMR::mr_heterogeneity(harm_2)
pleio  <- TwoSampleMR::mr_pleiotropy_test(harm_2)

egger <- TwoSampleMR::mr_pleiotropy_test(harm_2) 

presso <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = harm_2,
  NbDistribution = 1000,
  SignifThreshold = 0.05
)





# 
# top_CA <- gwas_CA %>%
#   filter(P < 5e-8) %>%
#   rename_with(~ paste0(.x, "_CA"), .cols = -SNP) %>%
#   left_join(gwas_chrono %>% rename_with(~ paste0(.x, "_chrono"), .cols = -SNP)) %>%
#   left_join(cojo_CA) 


p1 <- harm %>%
  left_join(cojo_CA, by = "SNP") %>%
  mutate(genes = str_remove_all(genes, "\\(.*?\\).*")) %>%
  ggplot(aes(x = beta.exposure, y = beta.outcome)) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome, ymax = beta.outcome + 1.96*se.outcome)) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96*se.exposure, xmax = beta.exposure + 1.96*se.exposure)) +
  ggrepel::geom_label_repel(aes(label = genes), alpha = 0.4, max.overlaps = 100) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  #scale_x_continuous(limits = c(-0.2, 0.2)) +
  theme_classic() +
  labs(title = "CA → Chronotype", x = "Beta CA (± CI)", y = "Beta chronotype (± CI)") +
  theme(text = element_text(size = 16)) 

p1

# top_chrono <- gwas_chrono %>%
#   filter(P < 5e-8) %>%
#   rename_with(~ paste0(.x, "_chrono"), .cols = -SNP) %>%
#   left_join(gwas_CA %>% rename_with(~ paste0(.x, "_CA"), .cols = -SNP)) %>%
#   left_join(cojo_chrono) 

p2 <- harm_2 %>%
  ggplot(aes(x = beta.exposure, y = beta.outcome)) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome, ymax = beta.outcome + 1.96*se.outcome)) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96*se.exposure, xmax = beta.exposure + 1.96*se.exposure)) +
  geom_smooth(method = "lm") +
  ggrepel::geom_label_repel(data = harm_2 %>% filter(abs(beta.exposure) > 0.03), aes(label = genes), alpha = 0.4, max.overlaps = 100) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  scale_x_continuous(limits = c(-0.1, 0.1)) +
  theme_classic() +
  labs(title = "Chronotype → CA", x = "Beta chronotype (± CI)", y = "Beta CA (± CI)") +
  theme(text = element_text(size = 16)) 

p2

px <- cowplot::plot_grid(p1, p2)

ggsave("OUTPUTS/f_MR.png", px, width = 16, height = 8)

### ACCELEREM


gwas_M10 <- data.table::fread("data/acc/ACC_M10_TIME.sumstats.gz") %>%
  filter(SNP %in% cojo_CA$SNP)


outcome_M10 <- TwoSampleMR::format_data(
  gwas_M10 %>% 
    as_tibble() %>% 
    mutate(Phenotype = "M10"),
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1",
  pval_col = "P",
  samplesize_col = "N"
)


harm_M10 <- TwoSampleMR::harmonise_data(
  exposure_dat,
  outcome_M10,
  action = 2
)


mr_M10 <- TwoSampleMR::mr(
  harm_M10,
  method_list = c(
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)

p3 <- harm_M10 %>%
  left_join(cojo_CA, by = "SNP") %>%
  mutate(genes = str_remove_all(genes, "\\(.*?\\).*")) %>%
  ggplot(aes(x = beta.exposure, y = beta.outcome)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome, ymax = beta.outcome + 1.96*se.outcome)) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96*se.exposure, xmax = beta.exposure + 1.96*se.exposure)) +
  ggrepel::geom_label_repel(aes(label = genes), alpha = 0.4, max.overlaps = 100, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  theme_classic() +
  labs(title = "CA → M10", x = "Beta CA (± CI)", y = "Beta M10 (± CI)") +
  theme(text = element_text(size = 12)) 
p3





gwas_L5 <- data.table::fread("data/acc/ACC_L5_TIME.sumstats.gz") %>%
  filter(SNP %in% cojo_CA$SNP)


outcome_L5 <- TwoSampleMR::format_data(
  gwas_L5 %>% 
    as_tibble() %>% 
    mutate(Phenotype = "L5"),
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1",
  pval_col = "P",
  samplesize_col = "N"
)


harm_L5 <- TwoSampleMR::harmonise_data(
  exposure_dat,
  outcome_L5,
  action = 2
)


mr_L5 <- TwoSampleMR::mr(
  harm_L5,
  method_list = c(
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)

p4 <- harm_L5 %>%
  left_join(cojo_CA, by = "SNP") %>%
  mutate(genes = str_remove_all(genes, "\\(.*?\\).*")) %>%
  ggplot(aes(x = beta.exposure, y = beta.outcome)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome, ymax = beta.outcome + 1.96*se.outcome)) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96*se.exposure, xmax = beta.exposure + 1.96*se.exposure)) +
  ggrepel::geom_label_repel(aes(label = genes), alpha = 0.4, max.overlaps = 100, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  theme_classic() +
  labs(title = "CA → L5", x = "Beta CA (± CI)", y = "Beta L5 (± CI)") +
  theme(text = element_text(size = 12)) 
p4





gwas_MIDP <- data.table::fread("data/acc/ACC_SLEEP_MIDP.sumstats.gz") %>%
  filter(SNP %in% cojo_CA$SNP)


outcome_MIDP <- TwoSampleMR::format_data(
  gwas_MIDP %>% 
    as_tibble() %>% 
    mutate(Phenotype = "MIDP"),
  type = "outcome",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "AF1",
  pval_col = "P",
  samplesize_col = "N"
)


harm_MIDP <- TwoSampleMR::harmonise_data(
  exposure_dat,
  outcome_MIDP,
  action = 2
)


mr_MIDP <- TwoSampleMR::mr(
  harm_MIDP,
  method_list = c(
    "mr_ivw_mre",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)

p5 <- harm_MIDP %>%
  left_join(cojo_CA, by = "SNP") %>%
  mutate(genes = str_remove_all(genes, "\\(.*?\\).*")) %>%
  ggplot(aes(x = beta.exposure, y = beta.outcome)) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = beta.outcome - 1.96*se.outcome, ymax = beta.outcome + 1.96*se.outcome)) +
  geom_errorbar(aes(xmin = beta.exposure - 1.96*se.exposure, xmax = beta.exposure + 1.96*se.exposure)) +
  ggrepel::geom_label_repel(aes(label = genes), alpha = 0.4, max.overlaps = 100, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  scale_y_continuous(limits = c(-0.1, 0.1)) +
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  theme_classic() +
  labs(title = "CA → MIDP", x = "Beta CA (± CI)", y = "Beta MIDP (± CI)") +
  theme(text = element_text(size = 12)) 
p5


P_ACC <- cowplot::plot_grid(p3, p4, p5, ncol = 3)

ggsave("OUTPUTS/f_MR_acc.png", P_ACC, width = 10, height = 4)








