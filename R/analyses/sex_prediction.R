library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(scales)
install.packages("cowplot")
install.packages("ggpmisc")

files <- list.files("/mnt/project/circadian/results/models/", pattern = "predictions", full.names = TRUE)

# 2. Compute per‐type R² summaries (r2s)
preds <- tibble(file = files) %>%
  mutate(data = map(file, readRDS)) %>% unnest(data)

preds_2 <- preds %>%
  rowwise() %>% mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>% ungroup() %>%
  mutate(sex = str_extract(file, "(?<=_)(fem|male|raw|tech)(?=_)"),
         type = str_extract(file, "(?<=predictions_)[^_]+")) %>%
  filter(sex %in% c("fem", "male"))


preds_2 %>%
  group_by(sex) %>%
  nest() %>%
  mutate(pmean = map_dbl(data, ~ cor(.x$time_day, .x$pred_mean)^2),
         lgb   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lgb)^2),
         xgboost = map_dbl(data, ~ cor(.x$time_day, .x$pred_xgboost)^2),
         lasso   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lasso)^2),
         lassox2 = map_dbl(data, ~ cor(.x$time_day, .x$pred_lassox2)^2))

