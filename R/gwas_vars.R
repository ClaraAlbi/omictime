

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

