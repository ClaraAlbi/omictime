
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("lightgbm")
library(lightgbm)
install.packages('xgboost', repos = c('https://dmlc.r-universe.dev', 'https://cloud.r-project.org'))
install.packages("glmnet")
library(xgboost)

time_i0 <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  filter(time_day >= 9 & time_day <= 20)

###

lgb1 <- lightgbm::lgb.load("data_share/cv.i0_lightgbm_cv1.rds")
xgb <- xgboost::xgb.load("data_share/cv.finngen_xgb_cv1.rds")
lasso <- readRDS("data_share/cv.i0_lasso_cv1.rds")
lassox2 <- readRDS("data_share/cv.i0_lassox2_cv1.rds")

### PREDS i0

l <- c(list.files("/mnt/project/circadian/results/models",
                  pattern = "predictions", full.names = TRUE))

preds_i0_olink <- tibble(f = l[str_detect(l, "tech_14")]) %>%
  mutate(d = map(f, readRDS)) %>%
  unnest(d) %>%
  rowwise() %>% mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  unnest()

### validation

i2_data <- data.table::fread("/mnt/project/OLINK_i2.tsv") %>%
  select(eid, any_of(readRDS("data_share/olink_panels_1to4_over.rds"))) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

row_na_counts <- rowSums(is.na(i2_data))
i2_data <- i2_data[row_na_counts < ncol(i2_data)/3,]

i2_meta <- data.table::fread("/mnt/project/blood_sampling_instance2.tsv") %>%
  filter(!is.na(`3166-2.0`)) %>%
  filter(eid %in% i2_data$eid) %>%
  mutate(max_time = pmax(`3166-2.0`,`3166-2.1`,`3166-2.2`,`3166-2.3`,`3166-2.4`, `3166-2.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date)

summary(i2_meta$y[!is.na(i2_meta$time_day)])

i2_imp <- glmnet::makeX(i2_data %>% select(-eid), na.impute = T)
i2_imp_x2 <- glmnet::makeX(i2_data %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i2 <- tibble(eid = i2_meta$eid[match(i2_data$eid, i2_meta$eid)],
                   y_test = i2_meta$time_day[match(i2_data$eid, i2_meta$eid)],
                   pred_lgb = predict(lgb1, as.matrix(i2_data %>% select(-eid))),
                   pred_xgboost = predict(xgb, as.matrix(i2_data %>% select(-eid))),
                   pred_lasso = predict(lasso, i2_imp)[,1],
                   pred_lassox2 = predict(lassox2, i2_imp_x2)[,1]) %>%
  left_join(i2_meta %>% select(eid, time_day, date_bsampling))

out_i2 <- preds_i2 %>%
  pivot_longer(-c(eid, y_test, time_day, date_bsampling)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value)^2),
         N = map_dbl(data, ~sum(!is.na(.x$value)))) %>%
  select(-data)




# i3
i3_data <- data.table::fread("/mnt/project/OLINK_i3.tsv") %>%
  select(eid, any_of(readRDS("data_share/olink_panels_1to4_over.rds"))) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

row_na_counts <- rowSums(is.na(i3_data))
i3_data <- i3_data[row_na_counts < ncol(i3_data)/3,]

i3_meta <- data.table::fread("/mnt/project/blood_sampling_instance3.tsv") %>%
  filter(!is.na(`3166-3.0`)) %>%
  #filter(eid %in% i3_data$eid) %>%
  mutate(max_time = pmax(`3166-3.0`,`3166-3.1`,`3166-3.2`,`3166-3.3`,`3166-3.4`, `3166-3.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date)

summary(i3_meta$y[!is.na(i3_meta$time_day)])

i3_imp <- glmnet::makeX(i3_data %>% select(-eid), na.impute = T)
i3_imp_x2 <- glmnet::makeX(i3_data %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i3 <- tibble(eid = i3_meta$eid[match(i3_data$eid, i3_meta$eid)],
                   y_test = i3_meta$time_day[match(i3_data$eid, i3_meta$eid)],
                   pred_lgb = predict(lgb1, as.matrix(i3_data %>% select(-eid))),
                   pred_xgboost = predict(xgb, as.matrix(i3_data %>% select(-eid))),
                   pred_lasso = predict(lasso, i3_imp)[,1],
                   pred_lassox2 = predict(lassox2, i3_imp_x2)[,1]) %>%
  left_join(i3_meta %>% select(eid, time_day, date_bsampling))

out_i3 <- preds_i3 %>%
  pivot_longer(-c(eid, y_test, time_day, date_bsampling)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value)^2),
         N = map_dbl(data, ~sum(!is.na(.x$value)))) %>%
  select(-data)


df <- preds_i0_olink %>% mutate(i = 0) %>%
  left_join(time_i0 %>% select(eid, date_bsampling)) %>%
  select(-f) %>%
  bind_rows(preds_i2 %>% mutate(i = 2) %>% select(-y_test)) %>%
  bind_rows(preds_i3 %>% mutate(i = 3)  %>% select(-y_test)) %>%
  group_by(eid) %>% mutate(n = n())  %>% ungroup()

saveRDS(df, "olink_int_replication.rds")


# df <- preds_i0_olink %>% mutate(i = 0) %>%
#   left_join(time_i0 %>% select(eid, date_bsampling)) %>%
#   bind_rows(readRDS("/mnt/project/olink_int_replication.rds") %>% filter(i!= 0) %>%
#               rowwise() %>% mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
#               unnest())
#
# saveRDS(df, "olink_internal_time_predictions.rds")


# MODELS NMR

time_i1 <- data.table::fread("/mnt/project/blood_sampling_instance1.tsv") %>%
  mutate(max_time = pmax(`3166-1.0`,`3166-1.1`,`3166-1.2`,`3166-1.3`,`3166-1.4`, `3166-1.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date) %>%
  filter(!is.na(time_day))

lgb1_nmr <- lightgbm::lgb.load("data_share/cv.NMR_lightgbm_cv1.rds")
xgb_nmr <- xgboost::xgb.load("data_share/cv.NMR_xgb_cv1.rds")
lasso_nmr <- readRDS("data_share/cv.NMR_lasso_cv1.rds")
lassox2_nmr <- readRDS("data_share/cv.NMR_lassox2_cv1.rds")

summary(time_i1$y[!is.na(time_i1$time_day)])

nmr_i1 <- data.table::fread("/mnt/project/nmr_nightingale_i1.tsv") %>%
  filter(!is.na(`20280-1.0`)) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

i1_imp <- glmnet::makeX(nmr_i1 %>% select(-eid), na.impute = T)
i1_imp_x2 <- glmnet::makeX(nmr_i1 %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i1 <- tibble(eid = time_i1$eid[match(nmr_i1$eid, time_i1$eid)],
                   y_test = time_i1$time_day[match(nmr_i1$eid, time_i1$eid)],
                   pred_lgb = predict(lgb1_nmr, as.matrix(nmr_i1 %>% select(-eid))),
                   pred_xgboost = predict(xgb_nmr, as.matrix(nmr_i1 %>% select(-eid))),
                   pred_lasso = predict(lasso_nmr, i1_imp)[,1],
                   pred_lassox2 = predict(lassox2_nmr, i1_imp_x2)[,1]) %>%
  filter(!is.na(eid)) %>%
  left_join(time_i1 %>% select(eid, time_day, date_bsampling))

out_i1 <- preds_i1 %>%
  pivot_longer(-c(eid, y_test, time_day, date_bsampling)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value, use = "complete.obs")^2),
         N = map_dbl(data, ~sum(!is.na(.x$value)))) %>%
  select(-data)

df_nmr <- preds_i0_nmr %>% mutate(i = 0, time = time_day) %>%
  inner_join(time_i0 %>% select(eid, date_bsampling)) %>%
  select(-f) %>%
  bind_rows(preds_i1 %>% mutate(i = 1) %>% select(-y_test)) %>%
  group_by(eid) %>% mutate(n = n()) %>% ungroup()

saveRDS(df_nmr, "nmr_int_replication.rds")

