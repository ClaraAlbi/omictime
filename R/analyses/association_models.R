library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
install.packages("broom")
install.packages("survminer")
library(broom)
library(survminer)
library(glue)

args <- commandArgs(trailingOnly = TRUE)
exposure <- args[1]

# --- 0. Prep: cohort, biomarkers & covariates, diseases ----------------

# administrative censor date
censor_date <- as.Date("2024-09-01")

# base cohort: death + baseline date
base_cohort <- fread("/mnt/project/clara/death.csv") %>%
  transmute(
    eid,
    death_date    = as.Date(p40000_i0)
  ) %>%
  left_join(
    readRDS("/mnt/project/biomarkers/time.rds") %>%
      mutate(date_bsampling = as.Date(date_bsampling)) %>%
      filter(date_bsampling > "1900-01-01") %>%
      select(eid, date_bsampling),
    by = "eid"
  )

# biomarker residuals + covs (sex, age, PCs, BMI, smoking, fasting)
time <- readRDS("/mnt/project/biomarkers/time.rds")
bio_covs <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(i == 0) %>%
  left_join(readRDS("/mnt/project/biomarkers/covs.rds"), by = "eid") %>%
  left_join(
    fread("/mnt/project/genetic_covs.tsv") %>%
      select(eid, `22009-0.1`:`22009-0.20`) %>%
      setnames(c("eid", paste0("PC", 1:20))),
    by = "eid"
  ) %>%
  left_join(time %>% select(eid, fasting)) %>%
  mutate(
    res        = residuals(lm(pred_lasso ~ time_day, data = cur_data())),
    res_abs    = abs(res),
    ares_q     = factor(ntile(res_abs, 5), levels = 1:5),
    gap        = pred_lasso - time_day,
    gap_abs    = abs(gap),
  ) %>%
  # ensure covariates are correct type
  mutate(
    sex           = factor(sex),
    smoking       = factor(smoking),
    age_recruitment = age_recruitment,
    BMI           = weight/((height/100)^2)
  ) %>% select(-date_bsampling) %>% filter(fasting < 24)

# diagnosis table
dis2 <- readRDS("/mnt/project/diseases_circadian.rds")

run_all_models <- function(disease_field) {
  print(disease_field)

  # ----- Shared setup -----
  df <- base_cohort %>%
    left_join(dis2 %>% transmute(eid, diag_date = .data[[disease_field]]),
              by = "eid") %>%
    left_join(bio_covs, by = "eid") %>%
    filter(!is.na(gap_abs)) %>%
    mutate(
      prev_case = as.integer(!is.na(diag_date) & diag_date < date_bsampling),
      end_date  = as.Date(pmin(diag_date, death_date, censor_date, na.rm = TRUE)),
      event     = as.integer(!is.na(diag_date) & diag_date <= end_date),
      time_yrs  = as.numeric(end_date - date_bsampling) / 365.25
    )

  # ----- Covariate setup -----
  base_covars   <- c("age_recruitment", paste0("PC", 1:20))
  extra_covars  <- c("BMI", "smoking", "fasting")
  if (disease_field != "date_breast") {
    base_covars <- c("sex", base_covars)
  }

  f_prev1 <- as.formula(paste("prev_case ~", exposure, " + ", paste(base_covars, collapse = " + ")))
  f_prev2 <- update(f_prev1, paste(". ~ . +", paste(extra_covars, collapse = " + ")))

  f_cox1 <- as.formula(paste("Surv(time_yrs, event) ~ ", exposure, " + ",  paste(base_covars, collapse = " + ")))
  f_cox2 <- update(f_cox1, paste(". ~ . +", paste(extra_covars, collapse = " + ")))

  results <- list()

  # ----- Prevalent Model 1 -----
  df_prev <- df %>% filter(!is.na(prev_case))
  if (sum(df_prev$prev_case, na.rm = TRUE) > 0) {
    m_prev1 <- glm(f_prev1, data = df_prev, family = binomial)
    res_prev1 <- broom::tidy(m_prev1, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == exposure) %>%
      transmute(
        disease = disease_field,
        model   = "Prevalent_Model1",
        type    = "logistic",
        term,
        effect  = estimate,
        lo95    = conf.low,
        hi95    = conf.high,
        p.value,
        cases   = sum(df_prev$prev_case == 1, na.rm = TRUE),
        n       = nrow(df_prev)
      )
    results <- append(results, list(res_prev1))
  }

  # ----- Prevalent Model 2 -----
  if (sum(df_prev$prev_case, na.rm = TRUE) > 0) {
    m_prev2 <- glm(f_prev2, data = df_prev, family = binomial)
    res_prev2 <- broom::tidy(m_prev2, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == exposure) %>%
      transmute(
        disease = disease_field,
        model   = "Prevalent_Model2",
        type    = "logistic",
        term,
        effect  = estimate,
        lo95    = conf.low,
        hi95    = conf.high,
        p.value,
        cases   = sum(df_prev$prev_case == 1, na.rm = TRUE),
        n       = nrow(df_prev)
      )
    results <- append(results, list(res_prev2))
  }

  # ----- Incident Model 1 (Cox) -----
  df_cox <- df %>%
    filter(is.na(death_date) | death_date > date_bsampling) %>%
    filter(is.na(diag_date)  | diag_date  > date_bsampling) %>%
    filter(time_yrs > 1)

  if (sum(df_cox$event, na.rm = TRUE) > 0) {
    m_cox1 <- survival::coxph(f_cox1, data = df_cox, x = FALSE, y = FALSE)
    res_cox1 <- broom::tidy(m_cox1, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == exposure) %>%
      transmute(
        disease = disease_field,
        model   = "Incident_Model1",
        type    = "cox",
        term,
        effect  = estimate,
        lo95    = conf.low,
        hi95    = conf.high,
        p.value,
        cases   = sum(df_cox$event == 1, na.rm = TRUE),
        n       = nrow(df_cox)
      )
    results <- append(results, list(res_cox1))
  }

  # ----- Incident Model 2 (Cox fully adjusted) -----
  if (sum(df_cox$event, na.rm = TRUE) > 0) {
    m_cox2 <- survival::coxph(f_cox2, data = df_cox, x = FALSE, y = FALSE)
    res_cox2 <- broom::tidy(m_cox2, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == exposure) %>%
      transmute(
        disease = disease_field,
        model   = "Incident_Model2",
        type    = "cox",
        term,
        effect  = estimate,
        lo95    = conf.low,
        hi95    = conf.high,
        p.value,
        cases   = sum(df_cox$event == 1, na.rm = TRUE),
        n       = nrow(df_cox)
      )
    results <- append(results, list(res_cox2))
  }

  bind_rows(results)
}

# select the diseases you want to plot
disease_cols <- setdiff(names(dis2), "eid")

all_results <- purrr::map_dfr(disease_cols, run_all_models)

saveRDS(all_results, glue("results_{exposure}_diseases.rds"))
