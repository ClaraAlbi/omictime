# Load models
lgb1 <- lightgbm::lgb.load("data_share/cv.finngen_lightgbm_cv1.rds")
xgb <- xgboost::xgb.load("data_share/cv.finngen_xgb_cv1.rds")
lasso <- readRDS("data_share/cv.finngen_lasso_cv1.rds")
lassox2 <- readRDS("data_share/cv.finngen_lassox2_cv1.rds")

# File with time of blood sampling info
time_file <- data.table::fread("/mnt/project/blood_sampling.tsv")
# Raw proteomics file
olink_raw_file <- data.table::fread("/mnt/project/olink_instance_0.csv")

# ============================================
# Get features from each model
# ============================================
# LASSO features (non-zero coefficients)
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

# ============================================
# Print summary of feature availability
# ============================================
cat("Feature availability summary:\n")
cat("LASSO: ", length(lasso_vars), "trained ->", length(common_lasso), "available\n")
cat("LASSO x²: ", length(lassox2_vars), "trained ->",
    length(common_lassox2_original) + length(common_lassox2_base), "base available\n")
cat("XGBoost: ", length(xgb_vars), "trained ->", length(common_xgb), "available\n")
cat("LightGBM: ", length(lgb_vars), "trained ->", length(common_lgb), "available\n\n")

cat("Missing features:\n")
cat("LASSO missing:", setdiff(tolower(lasso_vars), available_features), "\n")
cat("XGBoost missing:", setdiff(tolower(xgb_vars), available_features), "\n")
cat("LightGBM missing:", setdiff(tolower(lgb_vars), available_features), "\n\n")

# ============================================
# Prepare data with only common features
# ============================================
# For LASSO: select only common features
finngen_olink_lasso <- olink_raw_file %>%
  select(eid, any_of(common_lasso)) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

# For LASSO x²: select common original features + create squared terms
finngen_olink_lassox2 <- olink_raw_file %>%
  select(eid, any_of(c(common_lassox2_original, common_lassox2_base))) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

# For XGBoost and LightGBM
finngen_olink_xgb <- olink_raw_file %>%
  select(eid, any_of(common_xgb)) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

finngen_olink_lgb <- olink_raw_file %>%
  select(eid, any_of(common_lgb)) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

# ============================================
# Process time data
# ============================================
time_day <- time_file %>%
  mutate(max_time = pmax(`3166-0.0`, `3166-0.1`, `3166-0.2`, `3166-0.3`, `3166-0.4`, `3166-0.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  filter(time_day >= 9 & time_day <= 20)

# ============================================
# Impute and prepare matrices for prediction
# ============================================
# LASSO: impute with only common features
finngen_olink_imp_lasso <- glmnet::makeX(finngen_olink_lasso %>% select(-eid), na.impute = T)

# Keep only the columns that were in the training and are available
available_lasso_cols <- intersect(colnames(finngen_olink_imp_lasso), lasso_vars)
finngen_olink_imp_lasso <- finngen_olink_imp_lasso[, available_lasso_cols, drop = FALSE]

# LASSO x²: create squared terms for common base variables
finngen_olink_lassox2_prepared <- finngen_olink_lassox2 %>%
  select(-eid) %>%
  mutate(across(any_of(common_lassox2_base), list(sq = ~ .^2), .names = "{.col}_sq"))

finngen_olink_imp_lassox2 <- glmnet::makeX(finngen_olink_lassox2_prepared, na.impute = T)

# Keep only the columns that were in the training and are available
available_lassox2_cols <- intersect(colnames(finngen_olink_imp_lassox2), lassox2_vars)
finngen_olink_imp_lassox2 <- finngen_olink_imp_lassox2[, available_lassox2_cols, drop = FALSE]

# ============================================
# Make predictions
# ============================================
out_finngen <- tibble(
  eid = time_day$eid[match(finngen_olink_lasso$eid, time_day$eid)],
  y_test = time_day$time_day[match(finngen_olink_lasso$eid, time_day$eid)],
  pred_lgb = predict(lgb1, as.matrix(finngen_olink_lgb %>% select(-eid))),
  pred_xgb = predict(xgb, xgboost::xgb.DMatrix(data = as.matrix(finngen_olink_xgb %>% select(-eid)))),
  pred_lasso = predict(lasso, finngen_olink_imp_lasso)[,1],
  pred_lassox2 = predict(lassox2, finngen_olink_imp_lassox2)[,1]
)

# ============================================
# Prediction accuracy estimation
# ============================================
pred <- out_finngen %>%
  pivot_longer(c(-y_test, -eid)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value, use = "complete.obs")^2)) %>%
  select(-data)

print(pred)
