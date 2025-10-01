library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(scales)
install.packages("ggpmisc")

files <- c(list.files("/mnt/project/biomarkers_3",
                      pattern = "predictions", full.names = TRUE)[-c(31:35, 1:5, 16:20)],
           list.files("/mnt/project/biomarkers_3/covariate_res/MODELS",
                      pattern = "predictions", full.names = TRUE))

# 2. Compute per‐type R² summaries (r2s)
r2s <- tibble(file = files) %>%
  mutate(data = map(file, readRDS),
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
  select(-data, -file) %>%
  group_by(type) %>%
  summarise(
    across(pmean:lassox2, ~ mean(.x), .names = "{col}_mr2"),
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
  unnest(df) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)))


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

# plot_data <- data_long %>%
#   left_join(best_m, by = "type") %>%
#   filter(model == best_model) %>%
#   mutate(type = factor(
#     type,
#     levels = c("all","olink","NMR","labs","counts"),
#     labels = c("All","Proteomics","Metabolomics", "Biochemistry","Cell counts")
#   ))

plot_data <- data_long %>%
  filter(model == "mean") %>%
  mutate(type = factor(
    type,
    levels = c("all","olink","NMR","labs","counts"),
    labels = c("All","Proteomics","Metabolomics", "Biochemistry","Cell counts")
  )) %>%
  group_by(type) %>%
  mutate(Nt = n())


# 2. Define the formula object
formula <- y ~ x

# 3. Build the plot, layering on the eq‐labels
pl <- plot_data %>%
  filter(time_day > 9 & time_day < 20) %>%
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
  #scale_y_continuous(breaks = c(9, 12, 15, 18, 21), limits = c(9,21)) +
  scale_y_continuous(breaks = c(10, 15, 20), limits = c(9, 20)) +
  #scale_x_continuous(breaks = c(9, 12, 15, 18, 21), limits = c(9,21)) +
  scale_x_continuous(breaks = c(10, 15, 20), limits = c(9, 20)) +
  scale_color_manual(values = c(
    "All"          = "gray",
    "Proteomics"   = "#76B041",
    "Metabolomics" = "#2374AB",
    "Biochemistry" = "#E85F5C",
    "Cell counts"  = "#8F3985"
  )) +
  labs(x     = "Recorded time of day",
    y     = "Predicted internal time") +
  theme_classic(base_size = 11) +
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
    high     = muted("green"),
    midpoint = 0,
    limits   = c(-1, 1),
    name     = expression(r)
  ) +
  coord_fixed() +
  labs(
    x     = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title = element_blank(), legend.position = "none"
  )


ggsave("plots/F2_heat.png", p_heat, width = 4, height = 4)


### HEAT across models

pred_cols_a <- grep("^pred_", colnames(data_ind), value = TRUE)
pred_cols <- str_remove(pred_cols_a, "pred_")

# Function to compute correlation matrix -> long table
cor_fun <- function(df) {
  mat <- cor(df[, paste0("pred_", pred_cols)], use = "pairwise.complete.obs")
  # melt to long
  df_long <- as.data.frame(as.table(mat)) %>%
    rename(pred1 = Var1, pred2 = Var2, cor = Freq)

  # keep lower triangle only
  df_long <- df_long %>%
    filter(as.numeric(factor(pred1, levels = pred_cols_a)) >=
             as.numeric(factor(pred2, levels = pred_cols_a)))

  return(df_long)
}

cor_by_type <- data_ind %>%
  filter(type %in% c("all","olink","NMR","labs","counts")) %>%
  mutate(type = factor(
    type,
    levels = c("all","olink","NMR","labs","counts"),
    labels = c("All","Proteomics","Metabolomics", "Biochemistry","Cell counts")
  )) %>%
  group_by(type) %>%
  group_modify(~ cor_fun(.x)) %>%
  ungroup()

# plot

heat_mod <- ggplot(cor_by_type, aes(x = pred1, y = pred2, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 2) +
  scale_fill_gradient2(low = "lightblue", mid = "white", high = "pink", midpoint = 0.9) +
  facet_grid(~ type) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title = element_blank(), legend.position = "none"
  )

ggsave("plots/F3_head_mod.png", heat_mod, width = 8, height = 3)

