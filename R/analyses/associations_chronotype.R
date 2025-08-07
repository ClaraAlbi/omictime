library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)
library(stringr)
install.packages("broom")

chronotype <- fread("/mnt/project/chronotype.tsv")

# administrative censor date
censor_date <- as.Date("2022-05-31")

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
  left_join(chronotype %>% select(eid, chronotype = `1180-0.0`) %>%
              mutate(chrono = case_when(chronotype == "Definitely a 'morning' person" ~ 5,
                                        chronotype == "More a 'morning' than 'evening' person" ~ 4,
                                        chronotype == "Do not know" ~ 3,
                                        chronotype == "More an 'evening' than a 'morning' person" ~ 2,
                                        chronotype == "Definitely an 'evening' person" ~ 1,
                                        TRUE ~ NA_integer_
                                        ),
                     chronotype = factor(chronotype, levels = c("Definitely a 'morning' person",
                                                                "More a 'morning' than 'evening' person",
                                                                "Do not know",
                                                                "More an 'evening' than a 'morning' person",
                                                                "Definitely an 'evening' person")))) %>%
  mutate(
    res        = residuals(lm(pred_lasso ~ time_day, data = cur_data())),
    res_abs    = abs(res),
    res_q     = ntile(res, 10),
    gap        = pred_lasso - time_day,
    gap_abs    = abs(gap),
    gap_q     = ntile(gap, 5)
  ) %>%
  # ensure covariates are correct type
  mutate(
    sex           = factor(sex),
    smoking       = factor(smoking),
    age_recruitment = age_recruitment,
    BMI           = weight/((height/100)^2)
  ) %>% select(-date_bsampling) %>% filter(fasting < 24)


acc_chr <- bio_covs %>%
  filter(!is.na(chrono)) %>%
  group_by(chrono) %>%
  mutate(total_c = n()) %>%
  group_by(chrono, chronotype, gap_q, total_c) %>%
  summarise(prop = n() / unique(total_c), .groups = "drop",
            mean_gap = mean(gap))

acc_chr %>% ggplot(aes(x = gap_q, y = prop, fill = chronotype)) +
  geom_col(position = "dodge") +
  labs(x = "Acceleration defined chronotype")



# diagnosis table
dis2 <- readRDS("/mnt/project/diseases_circadian.rds")

run_all_models <- function(disease_field) {
  print(disease_field)


  df <- base_cohort %>%
    left_join(dis2 %>% transmute(eid, diag_date = .data[[disease_field]]),
              by = "eid") %>%
    left_join(bio_covs, by = "eid") %>%
    filter(!is.na(gap_abs)) %>%
    mutate(
      prev_case = as.integer(!is.na(diag_date) & diag_date < date_bsampling),
      end_date  = as.Date(pmin(diag_date, death_date, censor_date, na.rm = TRUE)),
      event     = as.integer(!is.na(diag_date) & diag_date <= end_date),
      time_yrs  = as.numeric(end_date - date_bsampling) / 365.25,
      chrono    = factor(chrono, levels = 1:5),
      gap_q    = factor(gap_q, levels = 1:5)
    )

  df$chrono <- relevel(df$chrono, ref = 5)
  df$gap_q <- relevel(df$gap_q, ref = 5)

  if (disease_field == "date_breast") {
    df <- df %>% filter(sex == 0 | sex == "Female")  # adjust based on your sex coding
  }

  # ----- Covariate setup -----
  base_covars   <- c("age_recruitment", paste0("PC", 1:20))
  extra_covars  <- c("BMI", "smoking", "fasting")
  if (disease_field != "date_breast") {
    base_covars <- c("sex", base_covars)
  }

  f_prev1 <- as.formula(paste("prev_case ~", "chrono + ", paste(base_covars, collapse = " + ")))
  f_prev2 <- as.formula(paste("prev_case ~", "gap_q + ", paste(base_covars, collapse = " + ")))
  f_cox1 <- as.formula(paste("Surv(time_yrs, event) ~ ", "chrono + ",  paste(base_covars, collapse = " + ")))
  f_cox2 <- as.formula(paste("Surv(time_yrs, event) ~ ", "gap_q + ",  paste(base_covars, collapse = " + ")))

  results <- list()

  # ----- Prevalent Model 1 -----
  df_prev <- df %>% filter(!is.na(prev_case))
  if (sum(df_prev$prev_case, na.rm = TRUE) > 0) {
    m_prev1 <- glm(f_prev1, data = df_prev, family = binomial)
    res_prev1 <- broom::tidy(m_prev1, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(str_detect(term, "chrono")) %>%
      transmute(
        disease = disease_field,
        model   = "Model1",
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

  if (sum(df_prev$prev_case, na.rm = TRUE) > 0) {
    m_prev2 <- glm(f_prev2, data = df_prev, family = binomial)
    res_prev2<- broom::tidy(m_prev2, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(str_detect(term, "gap")) %>%
      transmute(
        disease = disease_field,
        model   = "Model2",
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

  df_cox <- df %>%
    filter(is.na(death_date) | death_date > date_bsampling) %>%
    filter(is.na(diag_date)  | diag_date  > date_bsampling) %>%
    filter(time_yrs > 1)

  if (sum(df_cox$event, na.rm = TRUE) > 0) {
    m_cox1 <- survival::coxph(f_cox1, data = df_cox, x = FALSE, y = FALSE)
    res_cox1 <- broom::tidy(m_cox1, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(str_detect(term, "chrono")) %>%
      transmute(
        disease = disease_field,
        model   = "Model1",
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

  if (sum(df_cox$event, na.rm = TRUE) > 0) {
    m_cox2 <- survival::coxph(f_cox2, data = df_cox, x = FALSE, y = FALSE)
    res_cox2 <- broom::tidy(m_cox2, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(str_detect(term, "gap")) %>%
      transmute(
        disease = disease_field,
        model   = "Model2",
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

disease_cols <- setdiff(names(dis2), "eid")
all_results <- purrr::map_dfr(disease_cols, run_all_models)


saveRDS(all_results, "results_chronotype_diseases_1y.rds")
