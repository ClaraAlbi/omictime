lasso_coefs <- coef(lasso, s = "lambda.min")
lasso_vars <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
lasso_vars <- lasso_vars[lasso_vars != "(Intercept)"]

# LASSO x2 features (non-zero coefficients)
lassox2_coefs <- coef(lassox2, s = "lambda.min")
lassox2_vars <- rownames(lassox2_coefs)[which(lassox2_coefs != 0)]
lassox2_vars <- lassox2_vars[lassox2_vars != "(Intercept)"]

# Separate original and squared variables for lassox2
lassox2_original <- lassox2_vars[!grepl("_sq$", lassox2_vars)]
lassox2_squared <- lassox2_vars[grepl("_sq$", lassox2_vars)]
lassox2_base_for_squared <- gsub("_sq$", "", lassox2_squared)

# XGBoost features
xgb_vars <- xgb.feature_names(xgb)

# LightGBM features
lgb_vars <- lightgbm::lgb.feature_names(lgb1)

# ============================================
# Get available features in new data
# ============================================
available_features <- tolower(names(olink_raw_file)[names(olink_raw_file) != "eid"])

# ============================================
# Find common features for each model
# ============================================
common_lasso <- intersect(tolower(lasso_vars), available_features)
common_xgb <- intersect(tolower(xgb_vars), available_features)
common_lgb <- intersect(tolower(lgb_vars), available_features)

# For lassox2, find common original and base variables for squared terms
common_lassox2_original <- intersect(tolower(lassox2_original), available_features)
common_lassox2_base <- intersect(tolower(lassox2_base_for_squared), available_features)
