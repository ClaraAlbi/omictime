library(tidyr)
library(dplyr)
library(stringr)
library(pur)

rint <- function(x) {
  ranks <- rank(x, ties.method = "average")

  # Calculate the rank-inverse normal transformation
  n <- length(ranks)
  qnorm((ranks - 0.5) / n)
}

olink_cohort <- data.table::fread("/mnt/project/olink_instance_0.csv", select = "eid")

olink_cohort %>%
  mutate(IID = eid) %>%
  rename(FID = eid) %>%
  data.table::fwrite(., "olink_cohort.txt", sep = "\t", quote = FALSE, row.names = FALSE)

covs <- data.table::fread("/mnt/project/covariates.tsv") %>%
  mutate(across(c(2,3, 6, 7, 9), as.factor),
         bmi = `21002-0.0` / (`50-0.0`/100)^2,
         smoking = case_when(`20116-0.0` == "-3" ~ NA_character_,
                             TRUE ~ `20116-0.0`)) %>% select(-"50-0.0", -"21002-0.0", -"20116-0.0")

colnames(covs) <- c("eid", "sex", "birth_year",  "age_recruitment",  "assessment_centre", "month_attending", "bmi", "smoking")

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20", `22006-0.0`, `22021-0.0`, `22000-0.0`)
colnames(gen_covs) <- c("eid", paste0("PC", 1:20), "is_white", "rel", "batch")


df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  select(eid, time_day, pred_mean) %>%
  left_join(covs) %>%
  left_join(gen_covs) %>%
  filter(is_white == 1) %>%
  filter(rel == 0)

df$res <- residuals(lm(pred_mean ~ time_day, data = df))

df %>% mutate(FID = eid) %>%
  select(FID, eid, sex, age_recruitment, batch, contains("PC")) %>%
  mutate(across(c(sex, batch), as.factor)) %>%
  rename(IID = eid) %>%
  data.table::fwrite(., "covariates.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### COVARIATES COJO SNPS pQTLS

# snps <- data.table::fread("/mnt/project/snps/subset_cojo_pqtls.raw") %>%
#   mutate(across(contains(":"), ~round(.x,0)))
#
# data_b %>% mutate(FID = eid) %>%
#   select(FID, eid, sex, age_recruitment, batch, contains("PC")) %>%
#   mutate(across(c(sex, batch), as.factor)) %>%
#   left_join(snps %>% select(FID, contains(":")) %>% rename(eid = FID)) %>%
#   rename(IID = eid) %>%
#   data.table::fwrite(., "covariates_cojo.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#

df %>% mutate(FID = eid, res_abs = abs(res)) %>%
  select(FID, eid, res, res_abs) %>%
  rename(IID = eid) %>%
  data.table::fwrite(., "phenotypes.txt", sep = "\t", quote = FALSE, row.names = FALSE)



### gcta

covs %>%
  left_join(gen_covs) %>%
  mutate(FID = eid) %>%
  select(FID, eid,  age_recruitment,contains("PC")) %>%
  rename(IID = eid) %>%
  data.table::fwrite(., "qcovar.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")


covs %>%
  left_join(gen_covs) %>%
  mutate(FID = eid) %>%
  select(FID, eid,  sex, batch) %>%
  mutate(batch = case_when(is.na(batch) ~ NA_character_,
                           TRUE ~paste0("b",batch)),
         batch = as.factor(batch)) %>%
  rename(IID = eid) %>%
  data.table::fwrite(., "covar.txt", sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")






### snps pQTLS


### NMR


covs <- data.table::fread("/mnt/project/covariates.tsv") %>%
  mutate(across(c(2,3, 6, 7, 9), as.factor),
         bmi = `21002-0.0` / (`50-0.0`/100)^2,
         smoking = case_when(`20116-0.0` == "-3" ~ NA_character_,
                             TRUE ~ `20116-0.0`)) %>% select(-"50-0.0", -"21002-0.0", -"20116-0.0")

colnames(covs) <- c("eid", "sex", "birth_year",  "age_recruitment",  "assessment_centre", "month_attending", "bmi", "smoking")

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20", `22006-0.0`, `22021-0.0`, `22000-0.0`)
colnames(gen_covs) <- c("eid", paste0("PC", 1:20), "is_white", "rel", "batch")

l <- list.files("/mnt/project/biomarkers_3", full.names = T)

data_b <- tibble(f = l[str_detect(l, "predictions")]) %>%
  mutate(d = map(f, readRDS),
         type = stringr::str_extract(f, "(?<=predictions_)([^_]+)")) %>%
  unnest(d) %>%
  pivot_wider(id_cols = c(eid, time_day), names_from = type, values_from = contains("pred")) %>%
  filter(!is.na(pred_lgb_NMR)) %>%
  select(eid, time_day, pred_lgb_NMR) %>%
  left_join(covs) %>%
  left_join(gen_covs) %>%
  filter(is_white == 1) %>%
  filter(rel == 0)

data_b$res <- residuals(lm(pred_lgb_NMR ~ time_day, data = data_b))

data_b %>% mutate(FID = eid) %>%
  select(FID, eid, sex, age_recruitment, batch, contains("PC")) %>%
  mutate(across(c(sex, batch), as.factor)) %>%
  rename(IID = eid) %>%
  data.table::fwrite(., "covariates.txt", sep = "\t", quote = FALSE, row.names = FALSE)

data_b %>% mutate(FID = eid, res_abs = rint(abs(res))) %>%
  select(FID, eid, res, res_abs) %>%
  rename(IID = eid) %>%
  data.table::fwrite(., "phenotypes_NMR.txt", sep = "\t", quote = FALSE, row.names = FALSE)

