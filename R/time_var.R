library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(lubridate)
library(ggplot2)
library(purrr)

time <- data.table::fread("/mnt/project/blood_sampling.tsv") %>%
  mutate(max_time = pmax(`3166-0.0`, `3166-0.1`, `3166-0.2`, `3166-0.3`, `3166-0.4`, `3166-0.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(fasting = `74-0.0`,
         num_bsamples = `68-0.0`,
         date_bsampling = date)

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
  mutate(gap = time_day - pred_lasso_olink) %>%
  filter(!is.na(pred_lasso_olink)) %>%
  select(eid, time_day, gap, pred_lasso_olink) %>%
  left_join(covs) %>%
  left_join(gen_covs) %>%
  filter(is_white == 1 & rel == 0)


data_b %>%
  mutate(absgap = abs(gap)) %>%
  select(eid, time_day, gap, absgap, sex, age_recruitment, batch, contains("PC")) %>%
  mutate(batch = paste0("b", batch))

