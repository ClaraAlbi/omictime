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
            mean_N = mean(N)) %>%
  pivot_longer(
    cols = ends_with("mr2") | ends_with("sdr2"),
    names_to = c("model", "metric"),
    names_pattern = "(.*)_(mr2|sdr2)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = metric, values_from = value)

pl <- data_ind %>%
  filter(type %in% c("all", "olink", "NMR", "labs", "counts")) %>%
  mutate(type = factor(type, levels = c("all", "olink", "NMR", "labs", "counts"))) %>%
  left_join(r2s %>% filter(model == "lgb")) %>%
  ggplot(aes(x = time_day, y = pred_lgb, color = type)) +
  geom_point() +
  facet_wrap(~type , ncol = 5) +
  facet_grid(~type + paste0("N: ", round(mean_N, 0)) + paste0("R2: ",round(mr2, 2), "(", round(sdr2, 2),")")) +
  labs(y = "Predicted omic time", x = "Time day") +
  theme_minimal() +
  scale_x_continuous(breaks = c(10, 15, 20)) +
  scale_color_manual(values = c("gray", "#76B041", "#2374AB", "#E85F5C", "#8F3985")) +
<<<<<<< HEAD
  labs(y = "Predicted omic time", x = "Recorded time of day") +
  theme(legend.position = "none", strip.text = element_text(size = 16, hjust = 0))
=======
  theme(legend.position = "none")
>>>>>>> parent of 704fa62 (figure 3)

ggsave("plot_lgb.png", pl, width = 14, height = 3)

