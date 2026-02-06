library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(scales)
library(stringr)
install.packages("cowplot")
install.packages("ggpmisc")

files <- list.files("/mnt/project/LASSO_exp/", pattern = "predictions", full.names = TRUE)

# 2. Compute per‐type R² summaries (r2s)
preds <- tibble(file = files) %>%
  mutate(data = map(file, readRDS),
         group = case_when(
           str_detect(file, "14panels")  ~ "14panels",
           str_detect(file, "fem_tech")  ~ "female",
           str_detect(file, "male_tech") ~ "male",
           str_detect(file, "raw") ~ "raw",
           str_detect(file, "tech_14") ~ "tech"
         )) %>% unnest(data) %>%
  rowwise() %>% mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>% ungroup()


preds %>%
  group_by(group) %>%
  nest() %>%
  mutate(pmean = map_dbl(data, ~ cor(.x$time_day, .x$pred_mean)^2),
         lgb   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lgb)^2),
         xgboost = map_dbl(data, ~ cor(.x$time_day, .x$pred_xgboost)^2),
         lasso   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lasso)^2),
         lassox2 = map_dbl(data, ~ cor(.x$time_day, .x$pred_lassox2)^2))

cor(d$`14panels`, d$tech)

d <- preds %>%
  pivot_wider(id_cols = c(eid, time_day), names_from = group, values_from = pred_mean)

d2 <- preds %>% filter(group == "14panels") %>% select(-file, -group)


saveRDS(d2, "olink_int_panels14.rds")
