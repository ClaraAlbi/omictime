library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
install.packages("broom")
library(broom)
# -------------------------------------------------------------------
# 0.  House-keeping helpers ----
censor_date <- as.Date("2025-09-01")

# base cohort: death + baseline
base_cohort <- fread("/mnt/project/clara/death.csv") %>%
  transmute(eid, death_date = as.Date(p40000_i0)) %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds")   %>%
              mutate(date_bsampling = as.Date(date_bsampling)) %>%
              filter(date_bsampling > "1900-01-01") %>%
              select(eid, date_bsampling),
            by = "eid")

# biomarker & covariates (gap, sex, age, PCs…)
bio_covs <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(i == 0) %>%
  left_join(readRDS("/mnt/project/biomarkers/covs.rds"), by = "eid") %>%
  left_join(fread("/mnt/project/genetic_covs.tsv") %>%
              select(eid, "22009-0.1":"22009-0.20") %>%
              setnames(c("eid", paste0("PC", 1:20))),
            by = "eid") %>%
  select(-date_bsampling) %>%
  mutate(res       = residuals(lm(pred_lasso ~ time_day, data = cur_data())),
         gap       = pred_lasso - time_day,
         abs_gap   = abs(gap),
         ares_q    = ntile(abs_gap, 5) |> factor(levels = 1:5))

# -------------------------------------------------------------------
# 1.  Diagnosis table: one column per disease ----
dis2 <- fread("top_diseases_IEFG.csv") %>%
  mutate(across(-eid, ~ ifelse(.x > "1902-02-02", as.Date(.x), NA)))  # keep only real dates

disease_cols <- setdiff(names(dis2), "eid")   # every diagnosis field
disease_cols <- c("130894-0.0", "130896-0.0", "130892-0.0", "130874-0.0", "130708-0.0", "130792-0.0", "131286-0.0", "131306-0.0")
# -------------------------------------------------------------------
# 2.  Function to build dataset + run Cox for ONE disease ----
fit_one <- function(disease_field) {

  # Merge diagnosis column and rename to diag_date
  df <- base_cohort %>%
    left_join(dis2 %>% select(eid, !!disease_field), by="eid") %>%
    rename(diag_date = !!sym(disease_field)) %>%
    left_join(bio_covs, by="eid") %>%
    filter(is.na(death_date) | death_date > date_bsampling,
           is.na(diag_date)  | diag_date  > date_bsampling) %>%
    mutate(
      end_date = as.Date(
        pmin(diag_date, death_date, censor_date, na.rm = TRUE)
      ),   # ← re-cast to Date here
      event    = as.integer(
        !is.na(diag_date) &
          diag_date <= pmin(death_date, censor_date, na.rm = TRUE)
      ),
      time_yrs = as.numeric(end_date - date_bsampling) / 365.25
    ) %>%
    filter(time_yrs > 1)

  # If no events, return NA rows so tibble still binds
  if (sum(df$event) == 0)
    return(tibble(disease = disease_field,
                  term = "abs(gap)",
                  hr = NA_real_, lo95 = NA_real_, hi95 = NA_real_,
                  p = NA_real_, events = 0, n = nrow(df),
                  model = list(NULL)))

  # Fit Cox
  cox_m <- coxph(Surv(time_yrs, event) ~ abs_gap + sex + age_recruitment,
                 data = df, x = FALSE, y = FALSE)

  # Extract HR for abs_gap
  hr_tbl <- broom::tidy(cox_m,
                        exponentiate = TRUE,
                        conf.int     = TRUE) %>%    # <— add this
    filter(term == "abs_gap") %>%
    transmute(
      disease = disease_field,
      term,
      hr      = estimate,
      lo95    = conf.low,
      hi95    = conf.high,
      p       = p.value,
      events  = sum(df$event),
      n       = nrow(df)
    )

  hr_tbl$model <- list(cox_m)   # nested list-col
  hr_tbl
}

# -------------------------------------------------------------------
# 3.  Run across all diseases in parallel (map) ----
hr_results <- map_dfr(disease_cols, fit_one)

# -------------------------------------------------------------------
# 4.  Inspect
hr_results %>%
  select(disease, hr, lo95, hi95, p, events, n) %>%
  arrange(p)
