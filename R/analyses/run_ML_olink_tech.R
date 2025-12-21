# ---- deps
library(dplyr)
library(tidyr)
library(glue)
library(Matrix)   # sparse-friendly ops if needed
library(lightgbm)
library(xgboost)
library(glmnet)

type <- "olink_tech"

data_all <- readRDS("/mnt/project/res_olink_tech.rds")
stopifnot(all(c("eid") %in% names(data_all)))

time <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  select(eid, time_day)

dat <- data_all %>%
  left_join(time, by = "eid") %>%
  relocate(eid, time_day)

# feature names once
feat <- setdiff(names(dat), c("eid", "time_day"))

# ---- CV splits (fixed & saved for reproducibility)
cv_subsets <- sample(1:5, nrow(dat), replace = TRUE)
saveRDS(cv_subsets, glue("cv.{type}_subsets.rds"))

# ---- helpers (fast + leak-free)
# train-centered/scaled imputation, applied to train/test consistently
train_scaler <- function(X) {
  mu <- colMeans(X, na.rm = TRUE)
  # robust sd: fallback to 1 when 0/NA
  s  <- suppressWarnings(apply(X, 2, sd, na.rm = TRUE))
  s[!is.finite(s) | s == 0] <- 1
  list(mu = mu, s = s)
}

apply_scaler <- function(X, scaler) {
  # impute NA with train mean, then scale
  for (j in seq_along(scaler$mu)) {
    xj <- X[, j]
    nas <- is.na(xj)
    if (any(nas)) xj[nas] <- scaler$mu[j]
    X[, j] <- (xj - scaler$mu[j]) / scaler$s[j]
  }
  X
}

mk_poly2 <- function(X) {
  X2 <- X ^ 2
  colnames(X2) <- paste0(colnames(X), "_sq")
  cbind(X, X2)
}

# ---- main loop
system.time({
  for (cv in 1:5) {
    idx_train <- which(cv_subsets != cv)
    idx_test  <- which(cv_subsets == cv)

    # y
    y_train <- dat$time_day[idx_train]
    y_test  <- dat$time_day[idx_test]

    # X (dense matrix once; rely on impute/scale per fold)
    X_train_raw <- as.matrix(dat[idx_train, feat, drop = FALSE])
    X_test_raw  <- as.matrix(dat[idx_test , feat, drop = FALSE])

    # scale/impute using TRAIN stats
    scaler <- train_scaler(X_train_raw)
    X_train <- apply_scaler(X_train_raw, scaler)
    X_test  <- apply_scaler(X_test_raw , scaler)

    # ===================== LIGHTGBM =====================
    # Use early stopping on a small validation slice (no heavy lgb.cv inside outer CV)
    ntr <- length(y_train)
    nval <- max(1L, floor(0.2 * ntr))
    val_ids <- seq_len(nval)
    tr_ids  <- (nval + 1L):ntr

    dtrain_lgb <- lgb.Dataset(data = X_train[tr_ids, , drop = FALSE], label = y_train[tr_ids])
    dvalid_lgb <- lgb.Dataset.create.valid(dtrain_lgb,
                                           data = X_train[val_ids, , drop = FALSE],
                                           label = y_train[val_ids])

    lgb_params <- list(
      objective = "regression",
      metric = "l2",
      feature_fraction = 0.8,
      bagging_fraction = 0.8,
      bagging_freq = 5,
      min_data_in_leaf = 50,
      min_sum_hessian_in_leaf = 1,
      lambda_l1 = 5,
      lambda_l2 = 2,
      min_gain_to_split = 8,
      num_threads = max(1L, parallel::detectCores() - 1L)
    )

    lgb_model <- lgb.train(
      params = lgb_params,
      data = dtrain_lgb,
      valids = list(valid = dvalid_lgb),
      nrounds = 5000,
      early_stopping_rounds = 100,
      verbose = -1
    )
    saveRDS(lgb_model, glue("cv.{type}_lightgbm_cv{cv}.rds"))

    # ===================== XGBOOST ======================
    dtrain_xgb <- xgb.DMatrix(X_train, label = y_train)
    dvalid_xgb <- xgb.DMatrix(X_train[val_ids, , drop = FALSE], label = y_train[val_ids])

    xgb_params <- list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      eta = 0.05,
      max_depth = 6,
      subsample = 0.8,
      colsample_bytree = 0.8,
      nthread = max(1L, parallel::detectCores() - 1L)
    )

    xgb_model <- xgb.train(
      params = xgb_params,
      data = dtrain_xgb,
      nrounds = 5000,
      watchlist = list(valid = dvalid_xgb),
      early_stopping_rounds = 100,
      verbose = 0
    )
    saveRDS(xgb_model, glue("cv.{type}_xgb_cv{cv}.rds"))

    # ===================== LASSO ========================
    # glmnet standardizes by default; we already scaled → turn it off
    lasso_cv <- cv.glmnet(
      x = X_train, y = y_train,
      alpha = 1, family = "gaussian",
      standardize = FALSE,
      nfolds = 5, grouped = TRUE,
      intercept = TRUE
    )
    saveRDS(lasso_cv, glue("cv.{type}_lasso_cv{cv}.rds"))

    # ===================== LASSO (with squares) =========
    X_train_x2 <- mk_poly2(X_train)
    X_test_x2  <- mk_poly2(X_test)

    lasso_x2_cv <- cv.glmnet(
      x = X_train_x2, y = y_train,
      alpha = 1, family = "gaussian",
      standardize = FALSE,
      nfolds = 5, grouped = TRUE,
      intercept = TRUE
    )
    saveRDS(lasso_x2_cv, glue("cv.{type}_lassox2_cv{cv}.rds"))

    # ===================== PREDICT ======================
    # LightGBM / XGBoost
    y_pred_lgb <- predict(lgb_model, X_test, raw = FALSE)
    y_pred_xgb <- predict(xgb_model, xgb.DMatrix(X_test))

    # LASSO
    y_pred_lasso   <- as.numeric(predict(lasso_cv,   s = "lambda.min", newx = X_test))
    y_pred_lassox2 <- as.numeric(predict(lasso_x2_cv, s = "lambda.min", newx = X_test_x2))

    out_pred <- tibble(
      eid        = dat$eid[idx_test],
      cv         = cv,
      time_day   = y_test,
      pred_lgb   = y_pred_lgb,
      pred_xgboost = y_pred_xgb,
      pred_lasso   = y_pred_lasso,
      pred_lassox2 = y_pred_lassox2
    )
    saveRDS(out_pred, glue("predictions_{type}_cv{cv}.rds"))

    # ===================== FEATURE IMPORTANCES ==========
    # LASSO coefs
    coefs_lasso   <- as.matrix(coef(lasso_cv, s = "lambda.min"))
    feats_lasso   <- rownames(coefs_lasso)[coefs_lasso[,1] != 0]
    w_lasso       <- coefs_lasso[coefs_lasso[,1] != 0, 1]

    coefs_lassox2 <- as.matrix(coef(lasso_x2_cv, s = "lambda.min"))
    feats_lassox2 <- rownames(coefs_lassox2)[coefs_lassox2[,1] != 0]
    w_lassox2     <- coefs_lassox2[coefs_lassox2[,1] != 0, 1]

    # XGB importance
    imp_xgb <- xgb.importance(model = xgb_model)
    feats_xgb <- imp_xgb$Feature
    w_xgb     <- imp_xgb$Gain

    # LGB importance
    imp_lgb <- lgb.importance(lgb_model)
    feats_lgb <- imp_lgb$Feature
    w_lgb     <- imp_lgb$Gain

    out_coefs <- tibble(
      model   = c(rep("LASSO", length(feats_lasso)),
                  rep("LASSOx2", length(feats_lassox2)),
                  rep("XGB", length(feats_xgb)),
                  rep("LGB", length(feats_lgb))),
      feature = c(feats_lasso, feats_lassox2, feats_xgb, feats_lgb),
      weight  = c(w_lasso, w_lassox2, w_xgb, w_lgb)
    )

    saveRDS(out_coefs, glue("coefs_{type}_cv{cv}.rds"))

    #quick sanity (optional):
    cat(glue("CV{cv}: R2 LGB={round(cor(y_test, y_pred_lgb)^2,3)}  ",
             "XGB={round(cor(y_test, y_pred_xgb)^2,3)}  ",
             "LASSO={round(cor(y_test, y_pred_lasso)^2,3)}  ",
             "LASSOx2={round(cor(y_test, y_pred_lassox2)^2,3)}\n"))
    if (cv %% 1 == 0) gc(FALSE)
  }
})
