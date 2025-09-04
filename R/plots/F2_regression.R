library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(scales)
install.packages("ggpmisc")

files <- c(list.files("/mnt/project/biomarkers_3",
                      pattern = "predictions", full.names = TRUE)[-c(31:35, 1:5)],
           list.files("/mnt/project/biomarkers_3/covariate_res/MODELS",
                      pattern = "predictions", full.names = TRUE))

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
  filter(type %in% c("all","olink","NMR", "labs", "counts")) %>%
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
    labels = c("All","Proteomics","Metabolomics", "Biochemistry","Cell counts")
  ))

# 2. Define the formula object
formula <- y ~ x

# 3. Build the plot, layering on the eq‐labels
pl <- plot_data %>%
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
    strip.text = element_text(size = 12, face = "bold", hjust = 0),
    axis.title = element_text(face = "bold"), legend.position = "none"
  )


ggsave("plots/F3_pred.png", pl, width = 10, height = 3)



#### Correlations

pred_wide <- plot_data %>%
  select(eid, time_day, type, prediction) %>%
  pivot_wider(
    names_from  = type,
    values_from = prediction
  )

# 2. Compute correlation matrix (pairwise complete.obs)
corr_mat <- pred_wide %>%
  select(-eid, -time_day) %>%
  cor(use = "pairwise.complete.obs")

# 3. Melt into long form for ggplot
types <- colnames(corr_mat)

corr_long <- corr_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("type1") %>%
  pivot_longer(-type1, names_to = "type2", values_to = "r") %>%
  mutate(
    type1 = factor(type1, levels = types),
    type2 = factor(type2, levels = types)
  ) %>%
  filter(as.integer(type1) >= as.integer(type2))

# 5. Plot heatmap of just that half
p_heat <- ggplot(corr_long, aes(x = type2, y = type1, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r)), size = 4) +
  scale_fill_gradient2(
    low      = muted("blue"),
    mid      = "white",
    high     = muted("red"),
    midpoint = 0,
    limits   = c(-1, 1),
    name     = expression(r)
  ) +
  coord_fixed() +
  labs(
    x     = NULL, y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid       = element_blank(),
    axis.ticks       = element_blank(), legend.position = "none"
  )


ggsave("plots/F2_heat.png", p_heat)
