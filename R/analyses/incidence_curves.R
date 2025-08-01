library(dplyr)
library(purrr)
library(survival)
#install.packages("survminer")
library(survminer)
library(ggpubr)


exposure <- "gap_abs"

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


# 2. Function to build the dataset and return a survminer plot
make_ci_plot <- function(disease_field) {
  df_plot <- base_cohort %>%
    left_join(dis2 %>% select(eid, !!sym(disease_field))
              %>% rename(diag_date = !!sym(disease_field)),
              by = "eid") %>%
    left_join(bio_covs %>% select(eid, gap_abs, ares_q), by = "eid") %>%
    filter(is.na(death_date) | death_date > date_bsampling,
           is.na(diag_date)  | diag_date  > date_bsampling) %>%
    mutate(
      end_date = as.Date(pmin(diag_date, death_date, censor_date, na.rm = TRUE)),
      event    = as.integer(!is.na(diag_date) &
                              diag_date <= pmin(death_date, censor_date, na.rm = TRUE)),
      time_yrs = as.numeric(end_date - date_bsampling) / 365.25,
      ares_q   = factor(ares_q, levels = 1:5)
    ) %>%
    filter(time_yrs > 1, ares_q %in% c("1","3","5"))

  # Fit survival curves
  fit <- survfit(Surv(time_yrs, event) ~ ares_q, data = df_plot)

  # Plot cumulative incidence = 1 – survival

  df_km <- survminer::surv_summary(fit, data = df_plot) %>%
    mutate(cuminc = 1 - surv)

  ggplot(df_km, aes(time, cuminc, color = strata)) +
    #geom_point(alpha = 0.3) +
    geom_errorbar(aes(ymin = 1 - lower, ymax = 1 - upper)) +

    labs(
      title = paste("Disease:", disease_field),
      x = "Years since baseline",
      y = "Cumulative incidence",
      color = "res_abs quintile"
    ) +
    scale_color_manual(values = c("steelblue", "forestgreen", "firebrick")) +
    theme_minimal()

}

# 3. Generate one plot per disease
ci_plots <- map(indicated_diseases, make_ci_plot)

# 4. Display them
# If you’re in an interactive session:
walk(ci_plots, print)

# Or, to arrange them in a grid:
library(gridExtra)
grid.arrange(grobs = map(ci_plots, `[[`, "plot"), ncol = 1)
