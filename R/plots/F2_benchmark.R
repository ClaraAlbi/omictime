
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

files <- c(list.files("/mnt/project/biomarkers_3",
                      pattern = "predictions", full.names = TRUE)[-c(31:35, 1:5, 16:20)],
           list.files("/mnt/project/biomarkers_3/covariate_res/MODELS",
                      pattern = "predictions", full.names = TRUE))

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



res <- data %>%
  group_by(type) %>%
  mutate(samples = sum(N)) %>%
  mutate(type = factor(type, levels = c("all", "olink", "NMR", "labs", "counts"), labels = c("All", "Proteomics", "Metabolomics", "Biochemistry", "Cell counts"))) %>%
  filter(!is.na(type)) %>%
  pivot_longer(c(-cv, -type, -samples, -N)) %>%
  group_by(type, name) %>%
  mutate(m_r2 = mean(value),
         sd_r2 = sd(value),
         N_cv = round(mean(samples), 0),
         name = factor(name, levels = c("pmean", "lasso", "lassox2", "lgb", "xgboost")))

saveRDS(res, "data_share/prediction_accuracy.rds")

pbenchmark <- res %>%
  ggplot(aes(x = type, y = value, fill = name)) +
  geom_col(aes(y = m_r2, fill = name),
          position = position_dodge(width = 0.7),
          width = 0.7, alpha = 0.2) +
  geom_jitter(color = "black", shape = 21,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
              size = 1) +
  scale_y_continuous(n.breaks = 8) +
  scale_x_discrete(expand = c(0.01, 0)) +
  facet_grid(~type + samples, scales = "free") +
  scale_fill_viridis_d() +
  labs(y = "R2", fill = "Model") +
  theme_minimal() +
  theme(text = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.key.size = unit(1.2, "lines"),
        legend.text     = element_text(size = 16),
        legend.title    = element_text(size = 18))

ggsave("plots/F3_model_benchmark.png", pbenchmark, width = 10, height = 5)

