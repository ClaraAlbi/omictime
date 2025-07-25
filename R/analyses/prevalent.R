install.packages("broom")
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

# PREVALENT

covs <- readRDS("/mnt/project/biomarkers/covs.rds")
gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20")
colnames(gen_covs) <- c("eid", paste0("PC", 1:20))

data_b <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(i == 0) %>%
  left_join(covs) %>% left_join(gen_covs)
data_b$res <- residuals(lm(pred_lasso ~ time_day, data = data_b))

fields <- data.table::fread("field.tsv")

outcomes <- tribble(~field_id, ~phen,
                    "130894", "Depressive_episode",
                    "130896", "Recurrent_depression",
                    "130892", "Bipolar",
                    "130874", "Schizophrenia",
                    "130708", "T2D",
                    "130792", "Obesity",
                    "131286", "Hypertension",
                    "131306", "IHD",
                    "131380", "Atherosclerosis")

dis2 <- data.table::fread("/mnt/project/vars_diseases_2.tsv") %>%
  select(eid, contains(outcomes$field_id)) %>%
  inner_join(data_b %>% select(eid, date_bsampling)) %>%
  mutate(across(starts_with("13"), ~ case_when(.x < date_bsampling ~ 1,
                                               is.na(.x) ~ 0), .names = "{.col}_prevalent"),
         across(starts_with("13"), ~ case_when(.x > date_bsampling ~ 1,
                                               is.na(.x) ~ 0), .names = "{.col}_incident")) %>%
  select(-contains("_prevalent_incident")) %>%
  select(eid, contains("prevalent"))



Absgap <- dis2 %>%
  pivot_longer(-eid) %>%
  group_by(name) %>%
  nest() %>%
  ungroup() %>%
  mutate(field_id = str_extract(name, "^\\d+")) %>%
  left_join(outcomes) %>%
  mutate(data = map(data, ~.x %>%
                      inner_join(data_b) %>%
                      filter(!is.na(res)) %>%
                      filter(!is.na(value)) %>%
                      filter(if_all(where(is.factor), ~ . %in% names(which(table(.) >= 50)))) %>%
                      mutate(across(where(is.factor), droplevels))),
         mod = map(data, ~broom::tidy(glm(value ~ abs_res, data = .x %>% mutate(abs_res = abs(res)), family = "binomial"), conf.int = FALSE)),
         n = map_dbl(data, ~sum(.x$value == 1, na.rm = T))) %>%
  select(field_id, phen, mod, n) %>%
  unnest(mod)
