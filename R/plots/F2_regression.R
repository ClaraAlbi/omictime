library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

files <- list.files("/mnt/project/biomarkers_3",
                    pattern = "predictions", full.names = TRUE)

# 2. Compute per‐type R² summaries (r2s)
r2s <- tibble(file = files) %>%
  mutate(
    data = map(file, readRDS),
    type = str_extract(file, "(?<=predictions_)[^_]+"),
    cv   = str_extract(file, "(?<=cv)\\d+"),
    N    = map_dbl(data, ~ sum(!is.na(.x[[4]]))),
    lgb   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lgb)^2),
    xgboost = map_dbl(data, ~ cor(.x$time_day, .x$pred_xgboost)^2),
    lasso   = map_dbl(data, ~ cor(.x$time_day, .x$pred_lasso)^2),
    lassox2 = map_dbl(data, ~ cor(.x$time_day, .x$pred_lassox2)^2)
  ) %>%
  select(-data, -file) %>%
  group_by(type) %>%
  summarise(
    across(lgb:lassox2, ~ mean(.x), .names = "{col}_mr2"),
    N = sum(N),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = ends_with("mr2"),
    names_to  = "model",
    names_pattern = "(.*)_mr2",
    values_to = "mr2"
  )


data_ind <- tibble(file = files) %>%
  filter(str_detect(file, "predictions")) %>%
  mutate(
    df   = map(file, readRDS),
    type = str_extract(file, "(?<=predictions_)[^_]+")
  ) %>%
  select(type, df) %>%
  unnest(df)


data_long <- data_ind %>%
  filter(type %in% c("all","olink","NMR","labs","counts")) %>%
  pivot_longer(
    cols      = starts_with("pred_"),
    names_to  = "model",
    values_to = "prediction"
  ) %>%
  mutate(model = str_remove(model, "pred_"))

best_m <- r2s %>%
  group_by(type) %>%
  slice_max(mr2, n = 1, with_ties = FALSE) %>%
  select(type,
         best_model = model,
         top_R2     = mr2,
         Nt         = N)

plot_data <- data_long %>%
  left_join(best_m, by = "type") %>%
  filter(model == best_model) %>%
  mutate(type = factor(
    type,
    levels = c("all","olink","NMR","labs","counts"),
    labels = c("All","Proteomics","Metabolomics","Biochemistry","Cell counts")
  ))

# 2. Define the formula object
formula <- y ~ x

# 3. Build the plot, layering on the eq‐labels
pl <- plot_data %>%
  slice(rdunif(100000, nrow(plot_data))) %>%
  ggplot(aes(x = time_day, y = prediction)) +
  # your raw points
  geom_point(aes(color = type), alpha = 0.5, size = 0.8) +
  # the overall regression line
  geom_smooth(
    method  = "lm",
    formula = formula,
    color   = "red",
    size    = 1.2,
    se      = FALSE
  ) +
  # facet by type
  facet_wrap(~ type, ncol = 5) +
  # add R² and adj. R² in top‐left
  ggpmisc::stat_poly_eq(
    aes(label = paste(after_stat(rr.label), sep = "*\", \"*")),
    formula = formula,
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 4,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping    = aes(label = paste("italic(n) ==", Nt)),
    formula = formula,
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 4,
    color = "black"
  ) +
  # scales & styling
  scale_x_continuous(breaks = c(10, 15, 20)) +
  scale_color_manual(values = c(
    "All"          = "gray",
    "Proteomics"   = "#76B041",
    "Metabolomics" = "#2374AB",
    "Biochemistry" = "#E85F5C",
    "Cell counts"  = "#8F3985"
  )) +
  labs(
    x     = "Recorded time of day",
    y     = "Predicted omic time") +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(face = "bold"), legend.position = "none"
  )


ggsave("plots/F3_pred.png", pl, width = 10, height = 3)

