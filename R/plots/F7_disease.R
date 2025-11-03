#remotes::install_github("lcpilling/ukbrapR")
library(ukbrapR)
library(tidyr)
library(dplyr)
install.packages("broom")
library(purrr)

ckd <- ukbrapR:::codes_df_ckd

ukbrapR:::codes_df_test

bp <- tribble(~condition, ~vocab_id, ~code,
              "bp", "ICD10","F31.2",
              "bp", "ICD9", "296",
              "bp","Read2", "E11..11",
              "bp","Read2","E114.00",
              "bp","Read2","E115.00",
              "bp","Read2","E117.00",
              "bp","Read2","Eu31.00",
              "bp", "CTV3","E1176",
              "bp", "ukb_noncancer", "1291",
              "scz", "ukb_noncancer", "1289")

diagnosis_list <- get_diagnoses(disease_icd10)

diagnosis_df <- get_df(diagnosis_list, group_by="condition")


df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  left_join(diagnosis_df)

df$res <- residuals(lm(pred_mean ~ time_day, data = df))

a <- df %>%
  pivot_longer(contains("prev")) %>%
  group_by(name) %>% nest() %>%
  mutate(m = map(data, ~glm(value ~ abs(res), data = .x, family = binomial)),
         p = map(m, broom::tidy)) %>%
  unnest(p)

summary(glm(episode_MDD_severe_with_psychosis_bin_prev ~ abs(res), data = df, family = binomial))




