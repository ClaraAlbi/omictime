
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

files <- list.files("/mnt/project/biomarkers_3/out/", full.names = T)[59:63]

data <- tibble(file = files) %>%
  mutate(
    data = map(file, readRDS),
    data = map(data, ~.x %>% rowwise() %>% mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)))),
    type = str_extract(file, "(?<=predictions_)[^_]+"),
    cv   = str_extract(file, "(?<=cv)\\d+"),
    N    = map_dbl(data, ~ sum(!is.na(.x[[4]]))),
    pmean = map_dbl(data, ~ cor(.x$time_day, .x$pred_mean)^2),
    lgb   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lgb)^2),
    xgboost = map_dbl(data, ~ cor(.x$time_day, .x$pred_xgboost)^2),
    lasso   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lasso)^2),
    lassox2 = map_dbl(data, ~ cor(.x$time_day, .x$pred_lassox2)^2)
  ) %>%
  select(-data, -file)

# type   cv        N pmean   lgb xgboost lasso lassox2
# <chr>  <chr> <dbl> <dbl> <dbl>   <dbl> <dbl>   <dbl>
#   1 reveal 1     10447 0.540 0.505   0.438 0.525   0.534
# 2 reveal 2     10574 0.534 0.500   0.432 0.520   0.529
# 3 reveal 3     10464 0.548 0.509   0.440 0.535   0.542
# 4 reveal 4     10298 0.550 0.515   0.449 0.537   0.543
# 5 reveal 5     10460 0.541 0.512   0.445 0.521   0.530
