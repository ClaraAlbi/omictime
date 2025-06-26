library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

l <- list.files("/mnt/project/biomarkers_3", full.names = T)

data_ind <- tibble(f = l[str_detect(l, "predictions")]) %>%
  mutate(d = map(f, readRDS)) %>%
  mutate(type = stringr::str_extract(f, "(?<=predictions_)([^_]+)"),
         cv = stringr::str_extract(f, "(?<=cv)\\d+"),
         N = map_dbl(d, ~sum(!is.na(.x[[4]])))) %>%
  select(-cv) %>%
  unnest(d)

r2s <- tibble(f = l[str_detect(l, "predictions")]) %>%
  mutate(d = map(f, readRDS),
         type = stringr::str_extract(f, "(?<=predictions_)([^_]+)"),
         cv = stringr::str_extract(f, "(?<=cv)\\d+"),
         N = map_dbl(d, ~sum(!is.na(.x[[4]]))),
         lgb = map_dbl(d, ~cor(.x$time_day, .x$pred_lgb)^2),
         xgboost = map_dbl(d, ~cor(.x$time_day, .x$pred_xgboost)^2),
         lasso = map_dbl(d, ~cor(.x$time_day, .x$pred_lasso)^2),
         lassox2 = map_dbl(d, ~cor(.x$time_day, .x$pred_lassox2)^2)) %>%
  select(-d, -f) %>%
  group_by(type) %>%
  summarise(across(lgb:lassox2, ~mean(.x), .names = "{col}_mr2"),
            across(lgb:lassox2, ~sd(.x), .names = "{col}_sdr2"),
            N = sum(N)) %>%
  pivot_longer(
    cols = ends_with("mr2") | ends_with("sdr2"),
    names_to = c("model", "metric"),
    names_pattern = "(.*)_(mr2|sdr2)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = metric, values_from = value)


data_long <- data_ind %>%
  filter(type %in% c("all", "olink", "NMR", "labs", "counts")) %>%
  pivot_longer(cols = starts_with("pred_"),
               names_to = "model",
               values_to = "prediction") %>%
  mutate(model = str_remove(model, "pred_"))

best_m <- r2s %>%
  group_by(type) %>%
  slice_max(mr2, n = 1, with_ties = FALSE) %>%
  select(type, best_model = model, top_R2 = mr2, Nt = N)

pl <- data_long %>%
  left_join(best_m, by = "type") %>%
  filter(model == best_model) %>%
  mutate(type = factor(type, levels = c("all", "olink", "NMR", "labs", "counts"),
                       labels = c("All", "Proteomics", "Metabolomics", "Biochemistry", "Cell counts"))) %>%
  ggplot(aes(x = time_day, y = prediction, color = type)) +
  geom_point(alpha = 0.5, size = 0.8) +
  facet_wrap(~type + paste0("N: ", Nt) + paste0("R2: ",round(top_R2, 2)), ncol = 5) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(10, 15, 20)) +
  scale_color_manual(values = c("gray", "#76B041", "#2374AB", "#E85F5C", "#8F3985")) +
  labs(y = "Predicted omic time", x = "Recorded time of day") +
  theme(legend.position = "none", strip.text = element_text(size = 16, hjust = 0))


ggsave("plots/F3_pred.png", pl, width = 10, height = 4)

