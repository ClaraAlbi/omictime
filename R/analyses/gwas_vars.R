

olink_cohort <- data.table::fread("/mnt/project/olink_instance_0.csv", select = "eid")

olink_cohort %>%
  mutate(IID = eid) %>%
  rename(FID = eid) %>%
  data.table::fwrite(., "olink_cohort.txt", sep = "\t", quote = FALSE, row.names = FALSE)

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
  filter(is_white == 1) %>%
  filter(rel == 0)

data_b %>% mutate(FID = eid) %>%
  select(FID, eid, sex, age_recruitment, batch, contains("PC")) %>%
  mutate(across(c(sex, batch), as.factor)) %>%
  rename(IID = eid) %>%
  data.table::fwrite(., "covariates.txt", sep = "\t", quote = FALSE, row.names = FALSE)

data_b %>% mutate(FID = eid, absgap = abs(gap)) %>%
  select(FID, eid, gap, absgap) %>%
  rename(IID = eid) %>%
  data.table::fwrite(., "phenotypes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
