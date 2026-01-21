library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)
library(stringr)
library(scales)
library(lubridate)
library(data.table)
install.packages("broom")
install.packages("forcats")
library(broom)
library(forcats)

d <- readRDS("/mnt/project/diseases_circadian.rds")

outcomes <- colnames(d)[-1]

outcome_labels <- c(
  ra_date = "Rheumatoid arthritis",
  hyp     = "Hypertension",
  dys     = "Dyslipidaemia",
  dep_ep  = "Depression (episode)",
  chd     = "Chronic heart disease",
  obe     = "Obesity",
  t2d     = "Type 2 diabetes",
  dem     = "Dementia",
  dep_rec = "Depression (recurrent)",
  bip     = "Bipolar disorder",
  scz     = "Schizophrenia",
  liv     = "Chronic liver disease"
)


# count how many diagnoses happened >15 years ago (per row)
d_long <- d %>%
  mutate(across(all_of(outcomes), as.Date)) %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds") %>% select(eid, date_bsampling)) %>%
  mutate(
    date_bsampling = ymd(date_bsampling),
    across(all_of(outcomes), ymd)
  ) %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "diagnosis",
    values_to = "earliest_dx_date"
  )

d_long <- d_long %>%
  mutate(
    case_15y = case_when(
      !is.na(earliest_dx_date) &
        earliest_dx_date >= date_bsampling %m-% years(15) &
        earliest_dx_date <= date_bsampling ~ 1L,

      is.na(earliest_dx_date) ~ 0L,

      TRUE ~ NA_integer_
    )
  )

d_15_cc <- d_long %>%
  select(eid, date_bsampling, diagnosis, case_15y) %>%
  pivot_wider(
    id_cols = c(eid, date_bsampling),
    names_from = diagnosis,
    values_from = case_15y
  )



##############


covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         smoking = case_when(smoking == -3 ~NA, TRUE ~ smoking),
         smoking = as.factor(smoking),
         sex = as.factor(sex),
         assessment_centre = as.factor(assessment_centre))

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20")
colnames(gen_covs) <- c("eid", paste0("PC", 1:20))

chrono <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid, chrono = `1180-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ "Definitely morning",
    chrono == 2 ~ "Rather morning",
    chrono == -1~ "Don't know",
    chrono == 3 ~ "Rather evening",
    chrono == 4 ~ "Definitely evening",
    TRUE ~ NA_character_))

df <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  filter(i == 0) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate( season = case_when(
    m %in% c("12", "01", "02") ~ "Winter",
    m %in% c("03", "04", "05") ~ "Spring",
    m %in% c("06", "07", "08") ~ "Summer",
    TRUE ~ "Fall"
  ))
df$res <- residuals(lm(pred_mean ~ time_day, data = df))

data <- df %>%
  left_join(covs) %>%
  left_join(gen_covs) %>%
  left_join(chrono)

data1 <- data %>%
  select(eid, res, sex, age_recruitment,assessment_centre,season, smoking, bmi, chrono, any_of(paste0("PC", 1:20))) %>%
  left_join(d_15_cc) %>%
  mutate(
    res_sd_group = case_when(
      res >  2  ~ "Accelerated",
      res < -2  ~ "Delayed",
      TRUE        ~ "Middle"
    ),
    res_sd_group = factor(res_sd_group, levels = c("Middle", "Accelerated", "Delayed"))
  )



covariates <- c("sex", "age_recruitment", "chrono", "season", "assessment_centre", paste0("PC", 1:20))

res_continuous <- map_dfr(outcomes, function(outcome) {
  fit <- glm(
    as.formula(
      paste(
        outcome, "~ res +",
        paste(covariates, collapse = " + ")
      )
    ),
    data1, family ="binomial"
  )


  tidy(fit, exponentiate = TRUE) %>%
    filter(term == "res") %>%
    mutate(outcome = outcome)
}) %>%
  mutate(
    conf.low  = exp(log(estimate) - 1.96 * std.error),
    conf.high = exp(log(estimate) + 1.96 * std.error),
    pfdr      = p.adjust(p.value)
  )



middle_grp <- "Middle"
extreme_grps <- c("Accelerated", "Delayed")

res_sd_extremes <- map_dfr(outcomes, function(outcome) {

  map_dfr(extreme_grps, function(g) {

    df <- data1 %>%
      filter(res_sd_group %in% c(g, middle_grp)) %>%
      mutate(g_case = if_else(res_sd_group == g, 1L, 0L))

    fit <- glm(
      as.formula(
        paste(
          outcome, "~ g_case +",
          paste(covariates, collapse = " + ")
        )
      ),
      data = df,
      family = "binomial"
    )

    tidy(fit, exponentiate = TRUE) %>%
      filter(term == "g_case") %>%
      mutate(
        outcome  = outcome,
        contrast = paste0(g, "_vs_Middle"),
        group    = g
      )
  })
}) %>%
  mutate(
    conf.low  = exp(log(estimate) - 1.96 * std.error),
    conf.high = exp(log(estimate) + 1.96 * std.error),
    pfdr      = p.adjust(p.value)
  )


# Continuous var
#
ns <- d_long %>%
  filter(eid %in% data1$eid) %>%
  rename(outcome = diagnosis) %>%
  group_by(outcome) %>%
  summarise(
    cases = sum(case_15y == 1, na.rm = TRUE),
    controls = sum(case_15y == 0, na.rm = TRUE),
    excluded = sum(is.na(case_15y))
  )


#### PLOT

case_counts_q <- data1 %>%
  select(res_sd_group, all_of(outcomes)) %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "outcome",
    values_to = "case"
  ) %>%
  group_by(outcome, res_sd_group) %>%
  summarise(
    cases = sum(case == 1, na.rm = TRUE),
    .groups = "drop"
  )



plot_df <-
  bind_rows(
  res_sd_extremes %>%
  left_join(
    case_counts_q,
    by = c("outcome", "group" = "res_sd_group")
  ),
  res_continuous %>%
  left_join(ns) %>% mutate(group = "Continuous")) %>%
  filter(cases > 10) %>%
  mutate(
    group = factor(group, levels = c("Continuous",  "Accelerated", "Delayed")),
    outcome_full = recode(outcome, !!!outcome_labels),
    count_label = paste0("n=", cases),
    fill_fdr = pfdr < 0.05)

p_qs <- ggplot(plot_df, aes(x = estimate, y = outcome_full, color = fct_rev(group),  alpha = fill_fdr, fill = fct_rev(group))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0,
    position = position_dodge(width = 0.65)
  ) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  geom_point(
    shape = 21,
    size = 1.5,
    stroke = 1,
    position = position_dodge(width = 0.65)
  ) +
  geom_text(
    aes(x = conf.high * 1.18, label = count_label),
    size = 3,
    position = position_dodge(width = 0.65),
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_fill_manual(
    values = c(
      "Delayed" = "#b2182b",
      "Accelerated" = "#2166ac",
      "Continuous" = "black"
    ), label = c("Delayed (CA < -2)",
                 "Accelerated (CA > 2)",
                 "Continuous CA")
  ) +
  scale_color_manual(
    values = c(
      "Continuous" = "black",
      "Accelerated" = "#2166ac",
      "Delayed" = "#b2182b"
    ), guide = "none"
  ) +
  labs(
    x = "Odds Ratio",
    y = NULL,
    color = NULL,
    fill = "Risk factor",
    alpha = "FDR < 5%"
  ) +
  guides(
    alpha = guide_legend(order = 2, nrow = 2, byrow = TRUE),
    fill  = guide_legend(order = 1, nrow = 3, byrow = TRUE, reverse  = TRUE, override.aes = list(size = 4))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 30, 10, 10),
    legend.position = c(0.6, 0.68),
    legend.justification = c(0, 0)
  )

ggsave("plots/FX_disease_q.png", p_qs, width = 7, height = 7)

