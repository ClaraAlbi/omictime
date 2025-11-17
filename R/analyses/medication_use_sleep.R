library(tidyr)
library(lubridate)
library(dplyr)
library(glue)
library(purrr)
library(ggplot2)
install.packages("forcats")
library(forcats)
install.packages("janitor")
install.packages("broom")
library(stringr)

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  filter(smoking != "-3") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
  )

pcs <- data.table::fread("/mnt/project/covariates.txt") %>%
  select(eid = 1, contains("PC"))


meds <- data.table::fread("/mnt/project/medications.tsv")

# https://www.sciencedirect.com/science/article/abs/pii/S1389945717302757?via%3Dihub
# Current medications were self-reported to the research nurse, and the participants were dichotomised according to whether they were taking:
#   sleep medication (sedatives and hypnotics),
#   any other psychotropic medication (mood stabilisers, antidepressants, and antipsychotics) or
#   antihypertensive medication (ACE inhibitors, angiotensin II antagonists, beta blockers, calcium channel blockers, and diuretics).

cohort <- meds %>%
  #slice(1:1000) %>%
  pivot_longer(-eid) %>%
  filter(!is.na(value)) %>%
  left_join(X41467_2019_9572_MOESM3_ESM %>% select(Category, `Medication ATC code`, `Coding a`, `Drug name`), by = c("value" = "Coding a")) %>%
  filter(str_detect(`Medication ATC code`, "N05C") | #sleep medication (sedatives and hypnotics)
           str_detect(`Medication ATC code`, "N03AG01") | # mood stabiliser
           str_detect(`Medication ATC code`, "N03AX09") |  # mood stabiliser
           str_detect(`Medication ATC code`, "N03AF01") |  # mood stabiliser
           str_detect(`Medication ATC code`, "N03AF02") |  # mood stabiliser
           str_detect(`Medication ATC code`, "N05AN01") |  # lithium

           str_detect(`Medication ATC code`, "N06A") | # N06A ANTIDEPRESSANTS


           str_detect(`Medication ATC code`, "C09A") | # ACE inhibitors, plain
           str_detect(`Medication ATC code`, "C09C") | # ANGIOTENSIN II RECEPTOR BLOCKERS (ARBs), PLAIN
           str_detect(`Medication ATC code`, "C07") | # BETA BLOCKING AGENTS
           str_detect(`Medication ATC code`, "C08") | #CALCIUM CHANNEL BLOCKERS
           str_detect(`Medication ATC code`, "C03")) #C03 DIURETICS

cohort_f<- cohort %>%
  mutate(
    N05C = str_detect(`Medication ATC code`, "^N05C"),
    N05A = str_detect(`Medication ATC code`, "^N05A"),
    N03 = str_detect(`Medication ATC code`, "^N03"),
    N06 = str_detect(`Medication ATC code`, "^N06"),
    C0 = str_detect(`Medication ATC code`, "^C0")
  ) %>%
  group_by(eid) %>%
  summarise(
    across(c(N05C,N05A, N03, N06, C0), ~ any(.x, na.rm = TRUE)),
    .groups = "drop"
  )

data <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  inner_join(readRDS("/mnt/project/olink_int_replication.rds") %>% select(-date_bsampling)) %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  left_join(covs) %>%
  left_join(pcs) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  left_join(cohort_f) %>%
  mutate(has_prescription = ifelse(eid %in% cohort$eid, 1, 0),
         antihypertensive = case_when(isTRUE(C0) ~ 1, is.na(C0) | isFALSE(C0) ~ 0),
         sleep_medication = case_when(isTRUE(N05C) ~ 1, is.na(N05C) | isFALSE(N05C) ~ 0),
         antidepressants = case_when(isTRUE(N06) ~ 1, is.na(N06) | isFALSE(N06) ~ 0),
         mood_stabiliser = case_when(isTRUE(N03) ~ 1, is.na(N03) | isFALSE(N03) ~ 0),
         lithium = case_when(isTRUE(N05A) ~ 1, is.na(N05A) | isFALSE(N05A) ~ 0),
         across(c(has_prescription, antihypertensive, sleep_medication, antidepressants, mood_stabiliser, lithium), as.factor))

data$res <- residuals(lm(pred_mean ~ time_day, data = data))

table(data$has_prescription)
# 0     1
# 34868 17678

data %>%

vars <- c("has_prescription",
          "antihypertensive",
          "sleep_medication",
          "antidepressants",
          "mood_stabiliser",
          "lithium")

covars <- c("sex", "age_recruitment", paste0("PC", 1:10))

results <- map_dfr(vars, function(v) {
  adj_vars <- if (v %in% c("sex", "age_recruitment")) paste0("PC", 1:10) else covars

  # Combine predictor + covariates safely
  rhs <- paste(c(v, adj_vars), collapse = " + ")
  f <- as.formula(paste("abs(res) ~ ", rhs))

  fit <- lm(f, data = data)

  broom::tidy(fit) %>%
    filter(str_detect(term, paste0("^", v))) %>%
    mutate(predictor = v, reference = FALSE) %>%
    {
      if (is.factor(data[[v]])) {
        ref <- tibble(
          term = paste0(v, levels(data[[v]])[1]),
          estimate = 0, std.error = NA, statistic = NA, p.value = NA,
          predictor = v, reference = TRUE
        )
        bind_rows(ref, .)
      } else .
    }
})

saveRDS(results, "data_share/results_associations_medication_CM.rds")



library(table1)
my_render_cont <- function(x){
  with(
    stats.apply.rounding(stats.default(x)),
    c(
      "",
      `Mean (SD)` = sprintf("%s (%s)", MEAN, SD),
      `Median [Q1, Q3]` = sprintf("%s [%s, %s]",
                                  MEDIAN, Q1, Q3)
    )
  )
}


tab_desc <- table1::table1(~ age_recruitment + sex + has_prescription + antihypertensive + sleep_medication + antidepressants + mood_stabiliser + lithium,
                           data = data,
                           render.cont = my_render_cont)


