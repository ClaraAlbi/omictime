library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(scales)
library(stringr)

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


predictions_other <- preds %>%
  group_by(group, cv) %>%
  nest() %>%
  mutate(pmean = map_dbl(data, ~ cor(.x$time_day, .x$pred_mean)^2),
         lgb   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lgb)^2),
         xgboost = map_dbl(data, ~ cor(.x$time_day, .x$pred_xgboost)^2),
         lasso   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lasso)^2),
         lassox2 = map_dbl(data, ~ cor(.x$time_day, .x$pred_lassox2)^2)) %>%
  select(-data)

cor(d$`14panels`, d$tech)

d <- preds %>%
  pivot_wider(id_cols = c(eid, time_day), names_from = group, values_from = pred_mean)

d2 <- preds %>% filter(group == "14panels") %>% select(-file, -group)

#saveRDS(d2, "olink_int_panels14.rds")

#### PLOT

fp <- predictions_other %>%
  pivot_longer(-c(group, cv)) %>%
  group_by(group, name) %>%
  mutate(m_r2 = mean(value),
         sd_r2 = sd(value),
         name = factor(name, levels = c("pmean", "lasso", "lassox2", "lgb", "xgboost")),
         group = factor(group, levels = c("raw", "tech", "female", "male", "14panels"),
                        labels = c("Raw", "Technical", "Female-only", "Male-only", "All"))) %>% ungroup() %>%
  ggplot(aes(x = group, y = value, fill = name)) +
  geom_col(aes(y = m_r2, fill = name),
           position = position_dodge(width = 0.7),
           width = 0.7, alpha = 0.2) +
  geom_jitter(color = "black", shape = 21,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
              size = 1) +
  scale_y_continuous(n.breaks = 8) +
  scale_x_discrete(expand = c(0.01, 0)) +
  facet_grid(~group, scales = "free") +
  scale_fill_viridis_d() +
  labs(y = "R2", fill = "Model") +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        legend.text     = element_text(size = 16),
        legend.title    = element_text(size = 18))


ggsave("plots/F3_model_benchmark_covariates.png", fp, width = 10, height = 5)



predictions_other %>%
  group_by(group) %>%
  summarise(m_r2 = across(pmean:lassox2, mean),
         sd_r2 = across(pmean:lassox2, sd))

# group    m_r2$pmean  $lgb $xgboost $lasso $lassox2 sd_r2$pmean    $lgb $xgboost  $lasso $lassox2
# <chr>         <dbl> <dbl>    <dbl>  <dbl>    <dbl>       <dbl>   <dbl>    <dbl>   <dbl>    <dbl>
#   1 14panels      0.677 0.627    0.642  0.663    0.670     0.00524 0.00616  0.00587 0.00552  0.00468
# 2 female        0.679 0.617    0.630  0.673    0.674     0.0105  0.0120   0.0113  0.00856  0.00984
# 3 male          0.655 0.594    0.606  0.648    0.647     0.0107  0.0103   0.0130  0.0101   0.0111
# 4 raw           0.687 0.637    0.650  0.666    0.674     0.00248 0.00304  0.00381 0.00233  0.00135
# 5 tech          0.677 0.627    0.642  0.664    0.669     0.00248 0.00294  0.00338 0.00267  0.00321
