#remotes::install_github("lcpilling/ukbrapR")
library(ukbrapR)
library(tidyr)
library(dplyr)
install.packages("broom")
library(purrr)

ckd <- ukbrapR:::codes_df_ckd

ukbrapR:::codes_df_test

bp <- tribble(~condition, ~vocab_id, ~code,
              # "bp", "ICD10","F31.2",
              # "bp", "ICD9", "296",
              # "bp","Read2", "E11..11",
              # "bp", "CTV3","E1176",
              "sleep", "ICD10", "G47",
              "insomnias", "ICD10", "G470",
              "hypersomnias", "ICD10", "G471",
              "sleep-wake schedule", "ICD10", "G472",
              "apnea", "ICD10", "G473",
              "narcolepsy", "ICD10", "G474",
              "other", "ICD10" , "G478",
              "unspecified", "ICD10", "G479",
              "sleep", "ukb_noncancer", "1123")


# G47.0 Disorders of initiating and maintaining sleep [insomnias]113
# G47.1 Disorders of excessive somnolence [hypersomnias]137
# G47.2 Disorders of the sleep-wake schedule35
# G47.3 Sleep apnoea3847
# G47.4 Narcolepsy and cataplexy24
# G47.8 Other sleep disorders140
# G47.9 Sleep disorder, unspecified

diagnosis_list <- get_diagnoses(bp)
diagnosis_df <- get_df(diagnosis_list, group_by="condition")


df <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  filter(i == 0) %>%
  left_join(diagnosis_df)

df$res <- residuals(lm(pred_mean ~ time_day, data = df))

a <- df %>%
  pivot_longer(ends_with("bin")) %>%
  group_by(name) %>% nest() %>%
  mutate(m = map(data, ~glm(value ~ res, data = .x, family = binomial)),
         p = map(m, broom::tidy),
         n = map_dbl(data, ~sum(.x$value == 1, na.rm = T))) %>%
  unnest(p) %>% select(-data, -m)





