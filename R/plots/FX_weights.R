# === packages ===
library(tidyr)
library(dplyr)
library(stringr)
library(data.table)
library(viridis)
library(UpSetR)
library(lightgbm)
library(xgboost)
library(glmnet)

# If using LightGBM objects, ensure lightgbm is installed in your environment
# install.packages("lightgbm")  # only if available in your environment; skip if not
# If not available, the script will try fallbacks.

# === helper readers / extractors ===
safe_readRDS <- function(path) {
  tryCatch(readRDS(path), error = function(e) {
    warning("Failed to read ", path, ": ", conditionMessage(e))
    NULL
  })
}

# glmnet-style coefficient extractor (works for cv.glmnet, glmnet, or objects with coef())
extract_glmnet_coefs <- function(obj, tol = 1e-8) {
  # try coef(obj); convert to matrix and take first column if multiple
  cf <- tryCatch({
    as.matrix(coef(obj))
  }, error = function(e) {
    # attempt to find 'glmnet' like element inside object
    possible <- list()
    for (nm in names(obj)) {
      if (inherits(obj[[nm]], "dgCMatrix") || inherits(obj[[nm]], "matrix") ||
          inherits(obj[[nm]], "matrix")) {
        possible[[nm]] <- tryCatch(as.matrix(obj[[nm]]), error = function(e) NULL)
      }
    }
    # pick first reasonable matrix
    possible <- compact(possible)
    if (length(possible) >= 1) possible[[1]] else stop("no coef matrix found")
  })
  if (ncol(cf) > 1) cf <- cf[, 1, drop = FALSE]
  vals <- as.numeric(cf[, 1])
  nm <- rownames(cf)
  tibble(feature = nm, weight = vals) %>%
    # exclude intercept in caller, but keep here for completeness
    mutate(absw = abs(weight)) %>%
    filter(absw > tol | feature == "(Intercept)") %>%
    select(feature, weight)
}

# LightGBM extractor using lightgbm::lgb.importance when available
extract_lightgbm_importance <- function(obj) {
  # Accept either lgb.Booster or list containing booster
  booster <- NULL
  if (!is.null(obj)) {
    if (inherits(obj, "lgb.Booster")) booster <- obj
    else {
      # search for element that is lgb.Booster
      for (nm in names(obj)) {
        if (inherits(obj[[nm]], "lgb.Booster")) {
          booster <- obj[[nm]]; break
        }
      }
    }
  }
  if (is.null(booster)) {
    # try to find a data.frame named feature_importance or similar
    if (is.list(obj)) {
      for (nm in names(obj)) {
        el <- obj[[nm]]
        if (is.data.frame(el) && all(c("Feature", "Gain") %in% colnames(el))) {
          return(tibble(feature = as.character(el$Feature),
                        weight = as.numeric(el$Gain)))
        }
      }
    }
    stop("No LightGBM booster or importance table found")
  }
  # use lgb.importance
  imp <- tryCatch({
    lightgbm::lgb.importance(booster)
  }, error = function(e) {
    stop("lgb.importance failed: ", conditionMessage(e))
  })
  tibble(feature = as.character(imp$Feature),
         weight = as.numeric(imp$Gain))
}

# xgboost extractor (if any xgb.Booster present)
extract_xgboost_importance <- function(obj) {
  booster <- NULL
  if (!is.null(obj)) {
    if (inherits(obj, "xgb.Booster")) booster <- obj
    else {
      for (nm in names(obj)) {
        if (inherits(obj[[nm]], "xgb.Booster")) {
          booster <- obj[[nm]]; break
        }
      }
    }
  }
  if (is.null(booster)) stop("No xgboost booster found")
  imp <- xgboost::xgb.importance(model = booster)
  # xgb.importance returns data.frame with Feature and Gain
  tibble(feature = as.character(imp$Feature),
         weight = as.numeric(imp$Gain))
}

# generic model -> weights wrapper: tries to identify type by filename and object class
extract_weights_from_model <- function(path, obj) {
  fname <- basename(path)
  model_label <- NA_character_

  # determine model label from filename
  if (str_detect(fname, regex("lassox2", ignore_case = TRUE))) model_label <- "LASSO_X2"
  else if (str_detect(fname, regex("lasso", ignore_case = TRUE))) model_label <- "LASSO"
  else if (str_detect(fname, regex("lightgbm|lgb", ignore_case = TRUE))) model_label <- "LIGHTGBM"
  else if (str_detect(fname, regex("xgboost|xgb", ignore_case = TRUE))) model_label <- "XGBOOST"
  else model_label <- "UNKNOWN"

  # try extraction based on model_label or object class
  res <- tryCatch({
    if (model_label %in% c("LASSO", "LASSO_X2")) {
      tbl <- extract_glmnet_coefs(obj)
      tbl %>% mutate(model = model_label, source = fname) %>% select(model, source, everything())
    } else if (model_label == "LIGHTGBM") {
      tbl <- extract_lightgbm_importance(obj)
      tbl %>% mutate(model = model_label, source = fname) %>% select(model, source, everything())
    } else if (model_label == "XGBOOST") {
      tbl <- extract_xgboost_importance(obj)
      tbl %>% mutate(model = model_label, source = fname) %>% select(model, source, everything())
    } else {
      # fallback heuristics by object class
      cls <- class(obj)
      if ("lgb.Booster" %in% cls) {
        tbl <- extract_lightgbm_importance(obj)
        tbl %>% mutate(model = "LIGHTGBM", source = fname) %>% select(model, source, everything())
      } else if (any(str_detect(cls, "xgb"))) {
        tbl <- extract_xgboost_importance(obj)
        tbl %>% mutate(model = "XGBOOST", source = fname) %>% select(model, source, everything())
      } else {
        # try coef fallback
        tbl <- extract_glmnet_coefs(obj)
        tbl %>% mutate(model = "LASSO_like", source = fname) %>% select(model, source, everything())
      }
    }
  }, error = function(e) {
    warning("Failed to extract from ", fname, ": ", conditionMessage(e))
    NULL
  })
  res
}

# === locate your model files ===
model_dir <- "/mnt/project/biomarkers_3/covariate_res/MODELS"
model_files <- list.files(model_dir, pattern = "^cv\\.all_.*\\.rds$", full.names = TRUE)[-16]

# optionally, restrict to the five folds of lasso/lassox2/lightgbm user supplied --
# but we can process all found files
message("Found model files: ", length(model_files))

# If you want to use the exact list you posted, you could set:
# model_files <- c(
#   "/mnt/project/biomarkers_3/covariate_res/MODELS/cv.all_lasso_cv1.rds", ...
# )

# === read & extract ===
weights_list <- list()

for (f in model_files) {
  message("Reading: ", f)
  obj <- safe_readRDS(f)
  if (is.null(obj)) next

  # extract cv index from filename, e.g. cv1
  cv <- str_extract(basename(f), "cv[0-9]+") %>% str_remove("^cv") %>% as.character()
  if (is.na(cv)) cv <- NA_character_

  wt_tbl <- extract_weights_from_model(f, obj)
  if (is.null(wt_tbl)) next

  # attach cv as column and reorder to match your expected layout
  wt_tbl <- wt_tbl %>%
    mutate(cv = cv) %>%
    rename(feature = feature, weight = weight, model = model) %>%
    select(model, cv, feature, weight, source)

  weights_list[[f]] <- wt_tbl
}

weights <- bind_rows(weights_list)

# drop intercept rows if present
weights <- weights %>% filter(feature != "(Intercept)" & feature != "Intercept")


weights %>% group_by(model, cv) %>% count() %>%
  group_by(model) %>% summarise(m = mean(n))

# apply numerical tolerance for coefficient zeroing (for LASSO-like models)
tol <- 1e-8
weights_filtered <- weights %>%
  group_by(model, cv) %>%
  mutate(abs_weight = abs(weight)) %>%
  ungroup() %>%
  filter(abs_weight > tol | model %in% c("LIGHTGBM","XGBOOST"))
# For boosting models we don't threshold here because Gains are non-negative and meaningful.

# compute per-model-feature summaries across folds (mean, sd, n)
weights_summary <- weights_filtered %>%
  group_by(model, feature) %>%
  summarise(mean_weight = mean(weight, na.rm = TRUE),
            sd_weight   = sd(weight, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  filter(n >= 5)   # keep features observed in at least 1 fold; we'll later require n == number_of_folds if needed

# If you only want features that appear in all 5 folds, require n == 5
num_folds <- weights %>% pull(cv) %>% unique() %>% length()
weights_summary_allfolds <- weights_summary %>% filter(n == num_folds)

# scaled importance per model (within-model scaling)
weights_summary_scaled <- weights_summary_allfolds %>%
  mutate(abs_mean = abs(mean_weight)) %>%
  group_by(model) %>%
  mutate(scaled_importance = abs_mean / max(abs_mean, na.rm = TRUE)) %>%
  ungroup()

# === join with fields metadata for prettier names ===
fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

imp <- weights_summary_scaled %>%
  left_join(fields %>% select(field_id, title, main_category) %>% mutate(field_id = as.character(field_id)),
            by = c("feature" = "field_id")) %>%
  mutate(phen = if_else(is.na(title), feature, title)) %>%
  select(model, feature, phen,main_category, mean_weight, sd_weight, n, scaled_importance)

imp %>% group_by(main_category, model) %>% count() %>%
  mutate(main_category = case_when(is.na(main_category) ~ "Proteomics",
                                   main_category == "220" ~ "Metabolomics",
                                   main_category == "100081" ~ "Cell_counts",
                                   main_category == "17518" ~ "Biochemistry")) %>%
  pivot_wider(names_from = main_category, values_from = n) %>%
  mutate(Total = sum(Proteomics, Metabolomics, Cell_counts,  Biochemistry))

# === UpSet: prepare sets of top features per model ===
top_sets <- imp %>%
  group_by(model) %>%
  mutate(rank = rank(-scaled_importance, ties.method = "first")) %>%
  # choose top K to include in set comparisons; use top_n or all features that were in all folds
  filter(rank <= 500) %>%
  summarise(features = list(unique(feature)), .groups = "drop")

feat_list <- setNames(top_sets$features, top_sets$model)

# convert to binary presence matrix
upset_data <- fromList(feat_list)

# drop features appearing in only 1 model
upset_data_filtered <- upset_data[rowSums(upset_data) > 1, ]

# save upset plot
png("plots/FS_weights_upset.png", width = 2000, height = 1500, res = 300)
upset(upset_data_filtered,
      nsets = ncol(upset_data_filtered),
      order.by = "freq",
      mainbar.y.label = "Intersection size",
      sets.x.label = "Set size")
dev.off()

