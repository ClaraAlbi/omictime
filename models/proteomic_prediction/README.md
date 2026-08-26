# Proteomic prediction model weights

This directory contains fitted models from cross-validation fold 1 of the 14-panel Olink proteomic time-of-day prediction analysis.

## Files

| File | Model | Format |
|---|---|---|
| `cv.olink_14_lasso_cv1.rds` | LASSO regression using the original protein predictors 
| `cv.olink_14_lassox2_cv1.rds` | LASSO regression using original and squared protein predictors 
| `cv.olink_14_lightgbm_cv1.lgb` | LightGBM regression model 
| `cv.olink_14_xgb_cv1.ubj` | XGBoost regression model 

## Loading the models

```r
library(glmnet)
library(lightgbm)
library(xgboost)

lasso_model <- readRDS("cv.olink_14_lasso_cv1.rds")
lasso_x2_model <- readRDS("cv.olink_14_lassox2_cv1.rds")
lightgbm_model <- lightgbm::lgb.load("cv.olink_14_lightgbm_cv1.lgb")
xgboost_model <- xgboost::xgb.load("cv.olink_14_xgb_cv1.ubj")
```

## Obtain predictions

The following example assumes that `X_new` is a numeric matrix containing the new samples in rows and the Olink protein predictors in columns. The matrix must use the same proteins, column order, transformations, imputation, and training-set scaling as the model input.

```r
# X_new: preprocessed numeric matrix with one row per participant
X_new <- as.matrix(new_proteomic_data)

# Extract the training feature names from the fitted LASSO model.
expected_features <- rownames(lasso_model$glmnet.fit$beta)

# Add any missing proteins as zero-filled columns, then place all columns in
# exactly the same order used to fit the model. Additional columns are dropped.
missing_features <- setdiff(expected_features, colnames(X_new))
if (length(missing_features) > 0) {
  X_new <- cbind(
    X_new,
    matrix(
      0,
      nrow = nrow(X_new),
      ncol = length(missing_features),
      dimnames = list(NULL, missing_features)
    )
  )
}
X_new <- X_new[, expected_features, drop = FALSE]

# The squared-feature LASSO expects the original predictors followed by
# their element-wise squares.
X_new_sq <- X_new ^ 2
colnames(X_new_sq) <- paste0(colnames(X_new), "_sq")
X_new_x2 <- cbind(X_new, X_new_sq)

predictions <- data.frame(
  lasso = as.numeric(
    predict(lasso_model, newx = X_new, s = "lambda.min")
  ),
  lasso_x2 = as.numeric(
    predict(lasso_x2_model, newx = X_new_x2, s = "lambda.min")
  ),
  lightgbm = as.numeric(
    predict(lightgbm_model, X_new, raw = FALSE)
  ),
  xgboost = as.numeric(
    predict(xgboost_model, xgboost::xgb.DMatrix(X_new))
  )
)

# Ensemble prediction used in the main analysis: mean of the
# four model predictions for each participant.
predictions$pred_mean <- rowMeans(
  predictions[, c("lasso", "lasso_x2", "lightgbm", "xgboost")]
)

head(predictions)
```

