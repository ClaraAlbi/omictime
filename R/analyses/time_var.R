library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(lubridate)
library(ggplot2)
library(purrr)

time <- data.table::fread("/mnt/project/clara/blood_sampling_i0.csv") %>%
  mutate(max_time = pmax(p3166_i0_a0, p3166_i0_a1, p3166_i0_a2, p3166_i0_a3, p3166_i0_a4, p3166_i0_a5, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = "T") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(fasting = p74_i0,
         num_bsamples = p68_i0,
         date_bsampling = date)

covs <- data.table::fread("/mnt/project/clara/covariates.csv") %>%
  mutate(across(c(2,3, 6, 7, 9), as.factor),
         bmi = p21002_i0 / (p50_i0/100)^2,
         smoking = case_when(p20116_i0 == "-3" ~ NA_character_,
                             TRUE ~ p20116_i0)) %>% select(-p50_i0, -p21002_i0, -p20116_i0)

colnames(covs) <- c("eid", "sex", "birth_year",  "age_recruitment",  "assessment_centre", "month_attending", "bmi", "smoking")

gen_covs <- data.table::fread("/mnt/project/clara/genetic_covs.csv") %>%
  select(eid, p22009_a1:p22009_a20, p22006, p22021, p22000)
colnames(gen_covs) <- c("eid", paste0("PC", 1:20), "is_white", "rel", "batch")


saveRDS(time, "time.rds")
saveRDS(covs, "covs.rds")
saveRDS(gen_covs, "gencovs.rds")


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


out <- data_b %>%
  mutate(FID = eid, absgap = abs(gap)) %>%
  select(FID, IID = eid, time_day, gap, absgap, sex, age_recruitment, batch, contains("PC")) %>%
  mutate(batch = paste0("b", batch))

data.table::fwrite(out, "phenotypes_gap.txt")
