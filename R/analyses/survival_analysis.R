library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
install.packages("broom")
install.packages("survminer")
library(broom)
library(survminer)

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
bio_covs <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(i == 0) %>%
  left_join(readRDS("/mnt/project/biomarkers/covs.rds"), by = "eid") %>%
  left_join(
    fread("/mnt/project/genetic_covs.tsv") %>%
      select(eid, `22009-0.1`:`22009-0.20`) %>%
      setnames(c("eid", paste0("PC", 1:20))),
    by = "eid"
  ) %>%
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
  ) %>% select(-date_bsampling)

# diagnosis table
dis2 <- fread("/mnt/project/top_diseases_IEFG.csv") %>%
  mutate(across(-eid, ~ ifelse(.x > "1902-02-02", as.Date(.x), NA)))  # keep only real dates

# select the diseases you want to plot
disease_cols <- setdiff(names(dis2), "eid")

fit_hr_cont <- function(disease_field) {
  df <- base_cohort %>%
    left_join(dis2 %>% transmute(eid, diag_date = .data[[disease_field]]),
              by = "eid") %>%
    left_join(bio_covs, by = "eid") %>%
    filter(
      is.na(death_date) | death_date > date_bsampling,
      is.na(diag_date)  | diag_date  > date_bsampling
    ) %>%
    mutate(
      end_date = as.Date(    # ← force Date class
        pmin(diag_date, death_date, censor_date, na.rm = TRUE)
      ),
      event    = as.integer(!is.na(diag_date) & diag_date <= end_date),
      time_yrs = as.numeric(end_date - date_bsampling) / 365.25
    ) %>%
    filter(time_yrs > 1)

  if (sum(df$event)==0) {
    return(tibble(
      disease = disease_field,
      model   = c("Model1","Model2"),
      term    = "gap_abs",
      hr      = NA_real_,
      lo95    = NA_real_,
      hi95    = NA_real_,
      p.value = NA_real_,
      events  = 0,
      n       = nrow(df)
    ))
  }
  # formulas
  f1 <- as.formula(paste0(
    "Surv(time_yrs,event) ~ gap_abs + sex + age_recruitment + ",
    paste0("PC",1:20, collapse = "+")
  ))
  f2 <- update(f1, . ~ . + BMI + smoking )

  m1 <- coxph(f1, data = df, x = FALSE, y = FALSE)
  m2 <- coxph(f2, data = df, x = FALSE, y = FALSE)

  tbl1 <- tidy(m1, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term=="gap_abs") %>%
    transmute(
      disease = disease_field,
      model   = "Model1",
      term,
      hr      = estimate,
      lo95    = conf.low,
      hi95    = conf.high,
      p.value = p.value,
      events  = sum(df$event),
      n       = nrow(df)
    )

  tbl2 <- tidy(m2, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term=="gap_abs") %>%
    transmute(
      disease = disease_field,
      model   = "Model2",
      term,
      hr      = estimate,
      lo95    = conf.low,
      hi95    = conf.high,
      p.value = p.value,
      events  = sum(df$event),
      n       = nrow(df)
    )

  bind_rows(tbl1, tbl2)
}

hr_continuous_table <- map_dfr(disease_cols, fit_hr_cont)
saveRDS(hr_continuous_table, "hr_results_abs_gap.rds")

