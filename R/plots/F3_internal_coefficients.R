library(glmnet)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
install.packages("forcats")

data_all <- readRDS("/mnt/project/biomarkers_3/covariate_res/OLINK/res_olink.rds") %>%
  inner_join(readRDS("/mnt/project/biomarkers_3/covariate_res/res_nmr.rds")) %>%
  inner_join(readRDS("/mnt/project/biomarkers_3/covariate_res/res_labs.rds")) %>%
  inner_join(readRDS("/mnt/project/biomarkers_3/covariate_res/res_counts.rds"))

col_na_counts <- colSums(is.na(data_all))/nrow(data_all)
row_na_counts <- rowSums(is.na(data_all))
data_all <- data_all[row_na_counts < 1000, col_na_counts < 0.1]

cv_subsets <- readRDS("data_share/cv.all_subsets.rds")

X_mat  <- glmnet::makeX(data_all %>% select(-eid), na.impute = TRUE)

time_i0 <- readRDS("/mnt/project/biomarkers/time.rds")
y_obs <- time_i0$time_day[match(data_all$eid, time_i0$eid)]

thresholds <- c(0, 0.01, 0.05, 0.1)

# ---- Config: adjust if needed ----
model_pattern <- "cv.all_lasso_cv"
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
  ix <- which(coef(model) != 0)
  cf <- as.numeric(coef(model)[ix, 1])
  nm <- rownames(coef(model))[ix]
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

# ---- compute common features and all feature lists per threshold ----
common_and_features_per_thr <- map(thresholds, function(thr) {
  sel_lists <- map(lasso_models, selected_features, thr = thr)
  common_feats <- reduce(sel_lists, intersect)
  list(
    threshold = thr,
    n_common = length(common_feats),
    common_features = common_feats,
    sel_lists = sel_lists  # store all individual fold feature lists
  )
})

# ---- compute R2 per model using its own selected features at each threshold ----
results <- map_dfr(seq_along(common_and_features_per_thr), function(i) {
  thr_info <- common_and_features_per_thr[[i]]
  thr <- thr_info$threshold
  n_common <- thr_info$n_common

  # For each model (fold) compute prediction using its own selected features at this threshold
  map_dfr(seq_along(lasso_models), function(midx) {
    mod <- lasso_models[[midx]]

    # Get features selected by THIS model at THIS threshold
    selected_feats <- thr_info$sel_lists[[midx]]
    n_features <- length(selected_feats)

    # Get test indices for this fold (cv_subsets already matches data_all)
    test_idx <- which(cv_subsets == midx)

    # If no test samples in this fold, skip
    if(length(test_idx) == 0) {
      return(tibble(
        model_idx = midx,
        threshold = thr,
        R2 = NA_real_,
        n_features = n_features,
        n_common = n_common,
        features = list(selected_feats)
      ))
    }

    # Get test data
    X_test <- X_mat[test_idx, , drop = FALSE]
    y_test <- y_obs[test_idx]

    # If no selected features, predictions become intercept-only
    if(n_features == 0) {
      cf_all <- as.numeric(coef(mod)[,1])
      pred <- rep(cf_all[1], length(test_idx))
    } else {
      pred <- predict_using_keep(mod, X_test, selected_feats)
    }

    tibble(
      model_idx = as.character(midx),
      threshold = thr,
      R2 = safe_cor2(y_test, pred),
      n_features = n_features,
      n_common = n_common,
      features = list(selected_feats)  # store feature names as list column
    )
  })
})

# Save the feature information for later use
feature_info <- list(
  results_with_features = results,
  common_features_per_threshold = map(common_and_features_per_thr, function(x) {
    list(threshold = x$threshold,
         n_common = x$n_common,
         common_features = x$common_features)
  })
)


# ---------------------------
# compute predictions using ONLY common features, evaluated on CV fold 1 (model = fold1)
# ---------------------------
common_preds_cv1 <- map_dfr(common_and_features_per_thr, function(thr_info) {
  thr <- thr_info$threshold
  common_feats <- thr_info$common_features
  n_common <- length(common_feats)

  fold_idx <- 1
  mod <- lasso_models[[fold_idx]]
  test_idx <- which(cv_subsets == fold_idx)

  if (length(test_idx) == 0) {
    return(tibble(
      model_idx = "∩",
      threshold = thr,
      R2 = NA_real_,
      n_features = n_common,
      n_common = n_common,
      features = list(common_feats)
    ))
  }

  X_test <- X_mat[test_idx, , drop = FALSE]
  y_test <- y_obs[test_idx]

  if (n_common == 0) {
    cf_all <- as.numeric(as.matrix(coef(mod))[, 1])
    pred <- rep(cf_all[1], length(test_idx))
  } else {
    pred <- predict_using_keep(mod, X_test, common_feats)
  }

  tibble(
    model_idx = "∩",
    threshold = thr,
    R2 = safe_cor2(y_test, pred),
    n_features = n_common,
    n_common = n_common,
    features = list(common_feats)
  )
})

# append the common_cv1 rows
results <- bind_rows(results, common_preds_cv1)


# prepare factors/labels
results <- results %>%
  mutate(
    threshold_label = ifelse(threshold == 0, "0", as.character(threshold)),
    model_idx = factor(model_idx))

# ---- Plot: grouped bars per threshold, colored by CV fold ----
pos <- position_dodge(width = 0.8)

results_mean <- results %>%
  filter(model_idx != "∩") %>%
  group_by(threshold_label, n_common) %>%
  summarise(m_R2 = mean(R2, na.rm = TRUE),
            mean_n_features = round(mean(n_features, na.rm = TRUE), 1),
            .groups = "drop")

# --- plot with facet by threshold ---
ordered_levels <- c("0", "0.01", "0.05", "0.1") %>% unique()

results <- results %>%
  mutate(threshold_label = factor(threshold_label, levels = ordered_levels))

results_mean <- results_mean %>%
  mutate(threshold_label = factor(threshold_label, levels = ordered_levels))

p <- ggplot(results, aes(x = model_idx, y = R2, fill = model_idx)) +
  geom_col(width = 0.7) +
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
      x = 3.5,
      y = m_R2,
      label = paste0("Mean R² = ", round(m_R2, 2))
    ),
    color = "black",
    vjust = -0.8,
    size = 2.5,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_text(
    aes(label = n_features, y = 0.1),
    color = "white",
    size = 2.5, angle = 90,
    fontface = "bold"
  ) +
  scale_fill_viridis_d(option = "H") +
  labs(
    subtitle = "LASSO threshold:",
    x = "CV fold",
    y = expression(R^2),
    fill = "CV fold"
  ) +
  facet_wrap(~threshold_label,
             nrow = 1) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  ylim(0, max(results$R2, na.rm = TRUE) * 1.05)

ggsave("plots/F3_LASSO_coef.png", p, width = 6, height = 3)




