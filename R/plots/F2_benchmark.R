
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

saveRDS(data, "data/model_benchmark.rds")

pbenchmark <- data %>%
  group_by(type) %>%
  mutate(samples = sum(N)) %>%
  mutate(type = factor(type, levels = c("all", "olink", "NMR", "labs", "counts"), labels = c("All", "Proteomics", "Metabolomics", "Biochemistry", "Cell counts"))) %>%
  filter(!is.na(type)) %>%
  pivot_longer(c(-cv, -type, -samples, -N)) %>%
  group_by(type, name) %>%
  mutate(m_r2 = mean(value),
         N_cv = round(mean(samples), 0)) %>%

  ggplot(aes(x = type, y = value, fill = name)) +
  geom_col(aes(y = m_r2, fill = name),
          position = position_dodge(width = 0.7),
          width = 0.7, color = NA) +
  geom_jitter(color = "black", shape = 21,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
              size = 2, alpha = 0.9) +
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

