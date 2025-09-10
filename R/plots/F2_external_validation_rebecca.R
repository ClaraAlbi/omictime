
library(tidyr)
library(dplyr)
library(ggplot2)
library(lightgbm)
library(xgboost)
library(glmnet)

lgb1 <- lightgbm::lgb.load("cv.olink_lightgbm_cv1.rds")
xgb <- xgboost::xgb.load("cv.olink_xgb_cv1.rds")
lasso <- readRDS("cv.olink_lasso_cv1.rds")
lassox2 <- readRDS("cv.olink_lassox2_cv1.rds")

# Raw protein file format columns
# eid, time_day, prot1, prot2, prot3 etc

olink_raw_file <- data.table::fread("/mnt/project/olink_instance_0.csv")

# Keep only proteins in
prot_order <- lasso$glmnet.fit$beta@Dimnames[[1]]

olink_scaled <- olink_raw_file %>%
  select(eid, time_day, any_of(tolower(prot_order))) %>%
  mutate(across(-c(eid, time_day),  ~scale(.x)[,1]))

# Impute values needed to use glmnet
olink_scaled_imp <- glmnet::makeX(olink_scaled %>% select(-eid, -time_day), na.impute = T)
olink_scaled_imp_x2 <- glmnet::makeX(olink_scaled %>% select(-eid, -time_day) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

# Individual-lebel predictions
out_pred <- tibble(eid = olink_scaled$eid,
                   time_day = olink_scaled$time_day,
                   pred_lgb = predict(lgb1, as.matrix(olink_scaled %>% select(-eid, -time_day))),
                   pred_xgb = predict(xgb, as.matrix(olink_scaled %>% select(-eid, -time_day))),
                   pred_lasso = predict(lasso, olink_scaled_imp)[,1],
                   pred_lassox2 = predict(lassox2, olink_scaled_imp_x2)[,1]) %>%
  rowwise() %>%
  mutate(pred_mean =  mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)),
         mod_sd = sd(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)))

saveRDS(out_pred, "predictions_internal_time.rds")


# Plot time progression
out_pred %>%
  ggplot(aes(x = time_day, y = pred_mean, color = eid)) + geom_point()

# Prediction accuracy estimation
out_pred %>%
  pivot_longer(c(-time_day, -eid)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor.test(.x$time_day, .x$value, use = "complete.obs")^2)) %>%
  select(-data)




