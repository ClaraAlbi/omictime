library(dplyr)
library(tidyr)
library(glue)

data_nmr <- readRDS("/mnt/project/biomarkers_3/covariate_res/res_olink.rds") %>%
  left_join(readRDS("/mnt/project/biomarkers_3/covariate_res/res_nmr.rds")) %>%
  left_join(readRDS("/mnt/project/biomarkers_3/covariate_res/res_labs.rds")) %>%
  left_join(readRDS("/mnt/project/biomarkers_3/covariate_res/res_counts.rds")) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

time <- readRDS("/mnt/project/biomarkers/time.rds")

data_f <- data_nmr %>%
  left_join(time %>% select(eid, time_day))

type <- "all"

cv_subsets <- sample(1:5, nrow(data_f), replace = T)
saveRDS(cv_subsets, glue("cv.{type}_subsets.rds"))
#system(glue("dx upload cv.{type}_subsets.rds --destination project-Gp18X60Jk8K49jYjyKqXYY9b:/biomarkers_3/out/"))

system.time(for (cv in 1:5) {
  ids <- which(cv_subsets != cv)
  data_cv <- data_f[ids,]
  x <- as.matrix(data_cv %>% select(-time_day, -eid))
  y <- as.matrix(data_cv$time_day)
  
  ### LIGHTGBM
  
  lgb_grid <- list(
    objective = "regression",
    metric = "l2",
    min_sum_hessian_in_leaf = 1,
    feature_fraction = 0.8,  # Optimized
    bagging_fraction = 0.8,  # Optimized
    bagging_freq = 5,
    min_data = 50,  # Optimized
    lambda_l1 = 5,  # Optimized
    lambda_l2 = 2,  # Optimized
    min_gain_to_split = 8,  # Optimized
    is_unbalance = TRUE
  )
  
  data_gbm <- lightgbm::lgb.Dataset(data = x, label = y)
  
  # Cross-validation
  m <- lightgbm::lgb.cv(
    params = lgb_grid,
    data = data_gbm,
    nrounds = 1000,
    early_stopping_rounds = 50,
    eval_freq = 20,
    nfold = 5,
    stratified = TRUE
  )
  
  best_iter <- m$best_iter
  
  # Train final model
  lgb_model <- lightgbm::lgb.train(
    params = lgb_grid,
    data = data_gbm,
    nrounds = best_iter,
    eval_freq = 20
  )
  
  saveRDS(lgb_model, glue("cv.{type}_lightgbm_cv{cv}.rds"))
  
  #XGBOOST
  dtrain <- xgboost::xgb.DMatrix(data = x, label = y)
  
  # Define parameters
  xgb_params <- list(
    objective = "reg:squarederror",
    eta = 0.1,
    max_depth = 6,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  
  # Train model
  xgb_model <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = 100)
  saveRDS(xgb_model, glue("cv.{type}_xgb_cv{cv}.rds"))
  
  
  ### LASSO
  x_imp <- glmnet::makeX(data_cv %>% select(-time_day, -eid), na.impute = T)
  lasso_model <- glmnet::cv.glmnet(x_imp, y, alpha = 1)
  
  # Get best lambda and selected variables
  best_lambda <- lasso_model$lambda.min
  
  saveRDS(lasso_model, glue("cv.{type}_lasso_cv{cv}.rds"))
  
  ### LASSO X2
  x_imp_x2 <- glmnet::makeX(data_cv %>% select(-time_day, -eid) %>% 
                              mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T) 
  
  lasso_model_x2 <- glmnet::cv.glmnet(x_imp_x2, y, alpha = 1)
  
  # Get best lambda and selected variables
  best_lambda_x2 <- lasso_model_x2$lambda.min
  
  saveRDS(lasso_model_x2, glue("cv.{type}_lassox2_cv{cv}.rds"))
  
  # TEST MODELS
  
  y_test <- data_f$time_day[which(cv_subsets == cv)]
  y_pred_lgb <- predict(lgb_model, newdata = as.matrix(data_f[which(cv_subsets == cv),] %>% select(-time_day, -eid)))
  y_pred_xgboost <- predict(xgb_model, newdata = as.matrix(data_f[which(cv_subsets == cv),] %>% select(-time_day, -eid)))
  xtest_imp <- glmnet::makeX(data_f[which(cv_subsets == cv),] %>% select(-time_day, -eid), na.impute = T)
  y_pred_lasso <- predict(lasso_model, s = best_lambda, newx = xtest_imp)[,1]
  
  xtest_imp_x2 <- glmnet::makeX(data_f[which(cv_subsets == cv),] %>% select(-time_day, -eid) %>%
                                  mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)
  y_pred_lasso_x2 <- predict(lasso_model_x2, s = best_lambda_x2, newx = xtest_imp_x2)[,1]
  
  out <- tibble(eid = data_f$eid[which(cv_subsets == cv)],
                cv = cv,
                time_day = y_test,
                pred_lgb = y_pred_lgb,
                pred_xgboost = y_pred_xgboost,
                pred_lasso = y_pred_lasso,
                pred_lassox2 = y_pred_lasso_x2)
  
  saveRDS(out, glue("predictions_{type}_cv{cv}.rds"))
  
  # SAVE FEATURES
  coef_lasso <- coef(lasso_model, s = best_lambda)
  
  selected_features_lasso <- rownames(coef_lasso)[coef_lasso[, 1] != 0]
  selected_weights_lasso <- as.numeric(coef_lasso[coef_lasso[, 1] != 0, ])
  
  coef_lasso_x2 <- coef(lasso_model_x2, s = best_lambda)
  
  selected_features_lassox2 <- rownames(coef_lasso_x2)[coef_lasso_x2[, 1] != 0]
  selected_weights_lassox2 <- as.numeric(coef_lasso_x2[coef_lasso_x2[, 1] != 0, ])
  
  xgb_importance <- xgboost::xgb.importance(model = xgb_model)
  selected_features_xgboost <- xgb_importance$Feature
  selected_weights_xgboost <- xgb_importance$Gain 
  
  lgb_importance <- lightgbm::lgb.importance(lgb_model)
  selected_features_lgl <- lgb_importance$Feature
  selected_weights_lgl <- lgb_importance$Gain  # Importance score
  
  out_coefs <- tribble(~model, ~features, ~weights,
                       "LASSO", selected_features_lasso, selected_weights_lasso, 
                       "LASSOx2", selected_features_lassox2, selected_weights_lassox2, 
                       "XGB", selected_features_xgboost, selected_weights_xgboost, 
                       "LGB", selected_features_lgl, selected_weights_lgl) %>%
    unnest()
  
  saveRDS(out_coefs, glue("coefs_{type}_cv{cv}.rds"))
  
  cor(y_test, y_pred_lgb)^2
  cor(y_test, y_pred_xgboost)^2
  cor(y_test, y_pred_lasso)^2
  cor(y_test, y_pred_lasso_x2)^2
  
})



