
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

l <- list.files("/mnt/project/biomarkers_3", full.names = T)

data <- tibble(f = l[str_detect(l, "predictions")]) %>%
  mutate(d = map(f, readRDS)) %>%
  mutate(type = stringr::str_extract(f, "(?<=predictions_)([^_]+)"),
         cv = stringr::str_extract(f, "(?<=cv)\\d+"),
         N = map_dbl(d, ~sum(!is.na(.x$pred_lasso))),
         lgb = map_dbl(d, ~cor(.x$time_day, .x$pred_lgb)^2),
         xgboost = map_dbl(d, ~cor(.x$time_day, .x$pred_xgboost)^2),
         lasso = map_dbl(d, ~cor(.x$time_day, .x$pred_lasso)^2),
         lassox2 = map_dbl(d, ~cor(.x$time_day, .x$pred_lassox2)^2)) %>%
  select(-d, -f)

saveRDS(data, "model_benchmark.rds")


data %>%
  mutate(type = factor(type, levels = c("all", "olink", "NMR", "labs", "counts"))) %>%
  filter(!is.na(type)) %>%
  pivot_longer(c(-cv, -type, -N)) %>%
  group_by(type, name) %>%
  mutate(m_r2 = mean(value),
         N_cv = round(mean(N), 0)) %>%
  ggplot(aes(x = type, y = value, fill = name)) +
  geom_col(aes(y = m_r2, fill = name),
           position = position_dodge(width = 0.7),
           width = 0.7, color = NA) +
  geom_jitter(aes(fill = name), color = "black", shape = 21,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
              size = 1, alpha = 0.9) +
  scale_y_continuous(n.breaks = 8) +
  scale_x_discrete(expand = c(0.01, 0)) +
  facet_grid(~type + N_cv, scales = "free") +
  labs(y = "R2") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())
