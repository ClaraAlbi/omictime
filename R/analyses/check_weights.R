library(purrr)
library(tidyr)

### LASSO res

data_nmr <- readRDS("/mnt/project/biomarkers_3/covariate_res/res_olink.rds") %>%
  mutate(across(-eid, ~scale(.x)[,1]))

time <- readRDS("/mnt/project/biomarkers/time.rds")

data_f <- data_nmr %>%
  left_join(time %>% select(eid, time_day))

type <- "olink"


subs <- readRDS("/mnt/project/biomarkers_3/cv.olink_subsets.rds")

lasso_weights <- tibble(f = l[str_detect(l, "cv.olink_lasso_")])


genes <- data.table::fread("genes.txt") %>%
  pull(symbol) %>% tolower() %>% unique()

genes[genes %in% colnames(data_nmr)]

l <- list.files("/mnt/project/biomarkers_3", full.names = T)

y_pred_lasso <- predict(lasso_model, s = best_lambda, newx = xtest_imp)[,1]


lasso_top <- tibble(f = l[str_detect(l, "coefs_olink_")]) %>%
  mutate(lasso_mod = map(f, readRDS)) %>%
  unnest(lasso_mod) %>%
  filter(model == "LASSO") %>%
  group_by(f) %>%
  slice_max(abs(weights), n = 20) %>%
  pull(features) %>% unique()

olink_raw_file %>% select(any_of(lasso_top)) %>%
  cor(.,use = "complete.obs")
