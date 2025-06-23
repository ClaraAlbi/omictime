library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("lightgbm")
library(lightgbm)
install.packages("xgboost")
install.packages("glmnet")

l <- list.files("/mnt/project/biomarkers_3", full.names = T)

data_i0 <- tibble(f = l[str_detect(l, "predictions_i0")]) %>%
  mutate(d = map(f, readRDS),
         type = stringr::str_extract(f, "(?<=predictions_)([^_]+)"),
         cv = stringr::str_extract(f, "(?<=cv)\\d+"),
         N = map_dbl(d, ~nrow(.x)),
         lgb = map_dbl(d, ~cor(.x$time_day, .x$pred_lgb)^2),
         xgboost = map_dbl(d, ~cor(.x$time_day, .x$pred_xgboost)^2),
         lasso = map_dbl(d, ~cor(.x$time_day, .x$pred_lasso)^2),
         lassox2 = map_dbl(d, ~cor(.x$time_day, .x$pred_lassox2)^2)) %>%
  select(-d, -f)

# MODELS
lgb1 <- readRDS("/mnt/project/biomarkers_3/cv.i0_lightgbm_cv1.rds")
xgb <- readRDS("/mnt/project/biomarkers_3/cv.i0_xgb_cv1.rds")
lasso <- readRDS("/mnt/project/biomarkers_3/cv.i0_lasso_cv1.rds")
lassox2 <- readRDS("/mnt/project/biomarkers_3/cv.i0_lassox2_cv1.rds")

lightgbm::lgb.save(lgb1, "cv.i0_lightgbm_cv1.rds")
xgboost::xgb.save(xgb, "cv.i0_xgb_cv1.rds")


### validation

i2_data <- data.table::fread("/mnt/project/OLINK_i2.tsv")
i2_meta <- data.table::fread("/mnt/project/blood_sampling_instance2.tsv") %>%
  filter(!is.na(`3166-2.0`)) %>%
  mutate(max_time = pmax(`3166-2.0`,`3166-2.1`,`3166-2.2`,`3166-2.3`,`3166-2.4`, `3166-2.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60)

summary(i2_meta$y[!is.na(i2_meta$time_day)])

i2_imp <- glmnet::makeX(i2_data %>% select(-eid), na.impute = T)
i2_imp_x2 <- glmnet::makeX(i2_data %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

out_i2 <- tibble(y_test = i2_meta$time_day[match(i2_data$eid, i2_meta$eid)],
                 pred_lgb = predict(lgb.load("cv.i0_lightgbm_cv1.rds"), as.matrix(i2_data %>% select(-eid))),
                 pred_xgb = predict(xgboost::xgb.load("cv.i0_xgb_cv1.rds"), as.matrix(i2_data %>% select(-eid))),
                 pred_lasso = predict(lasso, i2_imp)[,1],
                 pred_lassox2 = predict(lassox2, i2_imp_x2)[,1]) %>%
  pivot_longer(-y_test) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value)^2)) %>%
  select(-data)


# i3
i3_data <- data.table::fread("/mnt/project/OLINK_i3.tsv")
i3_meta <- data.table::fread("/mnt/project/blood_sampling_instance3.tsv") %>%
  filter(!is.na(`3166-3.0`)) %>%
  mutate(max_time = pmax(`3166-3.0`,`3166-3.1`,`3166-3.2`,`3166-3.3`,`3166-3.4`, `3166-3.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60)

summary(i3_meta$y[!is.na(i3_meta$time_day)])

i3_imp <- glmnet::makeX(i3_data %>% select(-eid), na.impute = T)
i3_imp_x2 <- glmnet::makeX(i3_data %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

out_i3 <-tibble(y_test = i3_meta$time_day[match(i3_data$eid, i3_meta$eid)],
                pred_lgb = predict(lgb.load("cv.i0_lightgbm_cv1.rds"), as.matrix(i3_data %>% select(-eid))),
                pred_xgb = predict(xgboost::xgb.load("cv.i0_xgb_cv1.rds"), as.matrix(i3_data %>% select(-eid))),
                pred_lasso = predict(lasso, i3_imp)[,1],
                pred_lassox2 = predict(lassox2, i3_imp_x2)[,1]) %>%
  pivot_longer(-y_test) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value)^2)) %>%
  select(-data)

## common 3

ids <- union(i2_data$eid, i3_data$eid)

data_i0_s <- tibble(f = l[str_detect(l, "predictions_i0")]) %>%
  mutate(d = map(f, readRDS),
         type = stringr::str_extract(f, "(?<=predictions_)([^_]+)"),
         cv = stringr::str_extract(f, "(?<=cv)\\d+"),
         d = map(d, ~.x[.x$eid %in% ids,]),
         N = map_dbl(d, ~nrow(.x))) %>%
  select(d) %>%
  unnest(d) %>%
  select(-eid, -cv) %>%
  pivot_longer(-time_day) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$time_day, .x$value)^2)) %>%
  select(-data)

### lambdas
plot(lasso, xvar = "lambda", label = TRUE)
lower_lambdas <- c()
pred_matrix <- predict(lasso, newx = i2_imp, s = lower_lambdas)

