library(glmnet)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)

panels_in <- readRDS("data_share/olink_panels_1to4_over.rds")

res_olink <- readRDS("/mnt/project/clara/res_olink.rds") %>%
  select(eid, any_of(panels_in))

X_mat  <- glmnet::makeX(res_olink %>% select(-eid), na.impute = TRUE)

time_i0 <- readRDS("/mnt/project/clara/time.rds")
y_obs <- time_i0$time_day[match(res_olink$eid, time_i0$eid)]

thresholds <- c(0, 0.01, 0.05, 0.1)

# ---- Config: adjust if needed ----
model_pattern <- "cv.i0_lasso_cv"   # will look for files like cv.i0_lasso_cv1.rds ... cv.i0_lasso_cv5.rds
model_dir <- "data_share"

rds_files <- list.files(model_dir, pattern = paste0("^", model_pattern, ".*\\.rds$"), full.names = TRUE)
lasso_models <- map(rds_files, readRDS)

# ---- helpers ----
safe_cor2 <- function(y, yhat) {
  if(all(is.na(yhat))) return(NA_real_)
  cv <- suppressWarnings(cor(y, yhat, use = "complete.obs"))
  if(is.na(cv)) NA_real_ else cv^2
}

# get feature names selected in a model at given threshold (exclude intercept)
selected_features <- function(model, thr = 0) {
  cf <- as.numeric(coef(model)[, 1])
  nm <- rownames(coef(model))
  names(cf) <- nm
  feats <- nm[-1][abs(cf[-1]) >= thr]   # >= thr means selected at this threshold
  feats
}

# predict by zeroing out everything except 'keep' features (keep is vector of predictor names)
predict_using_keep <- function(model, X, keep) {
  cf <- as.numeric(coef(model)[,1])
  nm <- rownames(coef(model))
  names(cf) <- nm
  # zero out predictors not in keep
  keep_set <- intersect(names(cf[-1]), keep)  # ensure they exist
  cf[-1][!names(cf[-1]) %in% keep_set] <- 0
  intercept <- cf[1]
  betas <- cf[-1]
  as.numeric(intercept + X %*% betas)
}

# ---- compute intersection per threshold ----
# For each threshold, compute per-model selected features, then the intersection across all models
sel_per_thr <- map(thresholds, function(thr) {
  sel_lists <- map(lasso_models, selected_features, thr = thr)
  common <- reduce(sel_lists, intersect)
  list(threshold = thr, sel_lists = sel_lists, common = common)
})

# ---- compute R2 per model using only the intersection features for each threshold ----
results <- map_dfr(seq_along(sel_per_thr), function(i) {
  thr_info <- sel_per_thr[[i]]
  thr <- thr_info$threshold
  common_feats <- thr_info$common
  n_common <- length(common_feats)

  # For each model compute prediction using only the intersection and its R2
  map_dfr(seq_along(lasso_models), function(midx) {
    mod <- lasso_models[[midx]]
    # If no common features, predictions become intercept-only
    if(n_common == 0) {
      # predict intercept only: get intercept and use that as prediction
      cf_all <- as.numeric(coef(mod)[,1]); pred <- rep(cf_all[1], nrow(X_mat))
    } else {
      pred <- predict_using_keep(mod, X_mat, common_feats)
    }
    tibble(
      model_idx = midx,
      threshold = thr,
      R2 = safe_cor2(y_obs, pred),
      n_common = n_common
    )
  })
})

# prepare factors/labels
results <- results %>%
  mutate(
    threshold_label = ifelse(threshold == 0, "0", as.character(threshold)),
    model_idx = factor(model_idx))

# ---- Plot: grouped bars per threshold, colored by CV fold, label = n_common ----
pos <- position_dodge(width = 0.8)

results_mean <- results %>%
  group_by(threshold_label, n_common) %>%
  summarise(m_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

# --- plot with facet, bars, mean line + mean text ---
ordered_levels <- c( "0", "0.01", "0.05", "0.1") %>% unique()

results <- results %>%
  mutate(threshold_label = factor(threshold_label, levels = ordered_levels),
         n_common = forcats::fct_reorder(factor(n_common), desc(n_common)))

results_mean <- results_mean %>%
  mutate(threshold_label = factor(threshold_label, levels = ordered_levels))

# --- plot ---
p <- ggplot(results, aes(x = threshold_label, y = R2, fill = model_idx)) +
  geom_col(position = pos, width = 0.7) +
  geom_hline(
    data = results_mean,
    mapping = aes(yintercept = m_R2),
    color = "black",
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = results_mean,
    mapping = aes(
      x = threshold_label,
      y = m_R2,
      label = round(m_R2, 2)
    ),
    color = "black",
    vjust = -1.2,
    size = 3.8,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_d(option = "H") +
  labs(subtitle = "Number of proteins:",
    x = "Coefficient threshold",
    y = expression(R^2),
    fill = "CV fold"
  ) +
  facet_grid(~n_common, scales = "free") +
  theme_classic(base_size = 13) +
  ylim(0, max(results$R2, na.rm = TRUE) * 1.25)


ggsave("plots/F3_LASSO_coef.png", p, width = 6, height = 3)
