library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
library(ggplot2)
install.packages("broom")


rs <- data.table::fread("data_share/supplementary_data1.csv") %>%
  mutate(pfdr = p.adjust(p.value_time_day)) %>%
  filter(pfdr < 0.05) %>%
  filter(r2_time_day > 0.01)

olink <- readRDS("/mnt/project/biomarkers_res/res_olink_tech.rds") %>%
  select(eid, any_of(rs$FID))

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight / (height/100)^2,
         age_group = case_when(
           age_recruitment >= 40 & age_recruitment < 50 ~ "40-50",
           age_recruitment >= 50 & age_recruitment < 60 ~ "50-60",
           age_recruitment >= 60 & age_recruitment < 70 ~ "60-70",
      TRUE ~ NA_character_
    ),
    age_group = factor(age_group, levels = c("40-50", "50-60", "60-70")),
    sex = factor(sex),
    smoking = factor(smoking),
    month_attending = factor(month_attending)) %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds")) %>%
  mutate(cos_time = cos(2 * pi * time_day / 24),
         sin_time = sin(2 * pi * time_day / 24))


d <- olink %>%
  left_join(covs)


###

run_time_interaction <- function(data, outcome, interaction_var, alpha = 0.05) {

  f_base <- base_formula(outcome)
  print(outcome)

  f_int <- as.formula(paste0(
    outcome,
    " ~ (cos_time + sin_time) * ", interaction_var,
    " + age_group + sex + bmi + smoking + fasting + month_attending"
  ))

  m0 <- lm(f_base, data = data)
  m1 <- lm(f_int, data = data)

  anova_res <- anova(m0, m1) %>%
    as.data.frame() %>%
    slice(2) %>%
    transmute(
      outcome = outcome,
      interaction = interaction_var,
      df = Df,
      f_stat = F,
      p_value = `Pr(>F)`
    )

  use_interaction <- anova_res$p_value < alpha

  est_res <- if (use_interaction) {
    broom::tidy(m1) %>%
      mutate(
        outcome = outcome,
        model = paste0("time_x_", interaction_var)
      )
  } else {
    broom::tidy(m0) %>%
      mutate(
        outcome = outcome,
        model = "baseline"
      )
  }

  list(
    anova = anova_res,
    estimates = est_res,
    interaction_kept = use_interaction
  )
}


biomarkers <- colnames(olink)[-1]

base_formula <- function(outcome) {
  as.formula(paste0(
    outcome,
    " ~ cos_time + sin_time + age_group + sex + ",
    "bmi + smoking + fasting + month_attending"
  ))
}

interaction_vars <- c("sex")

results <- map(
  biomarkers,
  function(biomarker) {
    map(
      interaction_vars,
      function(int_var) {
        run_time_interaction(
          data = d,
          outcome = biomarker,
          interaction_var = int_var)})})

anova_table <- map_dfr(
  results,
  ~ map_dfr(.x, "anova")
) %>% mutate(pfdr = p.adjust(p_value)) %>%
  filter(pfdr < 0.05)

estimate_table <- map_dfr(
  results,
  ~ map_dfr(.x, "estimates")
) %>% filter(outcome %in% anova_table$outcome)






