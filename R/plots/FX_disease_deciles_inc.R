# Required packages (load if not already)
library(dplyr)
library(tidyr)
library(lubridate)
library(broom)
library(purrr)
library(forcats)
library(ggplot2)
library(scales)
library(table1)

# ---- 1) Read input disease file and create long table of earliest diagnosis dates ----
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

# Convert outcome columns to Date and pivot longer
d_long_future <- d %>%
  mutate(across(all_of(outcomes), as.Date)) %>%
  # join sampling date
  left_join(readRDS("/mnt/project/biomarkers/time.rds") %>% select(eid, date_bsampling), by = "eid") %>%
  mutate(
    date_bsampling = ymd(date_bsampling),
    across(all_of(outcomes), ymd)
  ) %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "diagnosis",
    values_to = "earliest_dx_date"
  )

# ---- 2) Create future-15-year case indicator ----
d_long_future <- d_long_future %>%
  mutate(
    # case_15y_future: 1 if diagnosis occurs AFTER baseline and within 15 years,
    # 0 if no date (kept same semantics as original code) -- adjust if you want NA for censored
    case_15y_future = case_when(
      !is.na(earliest_dx_date) &
        earliest_dx_date > date_bsampling &
        earliest_dx_date <= date_bsampling %m+% years(15) ~ 1L,

      is.na(earliest_dx_date) ~ 0L,

      TRUE ~ NA_integer_
    )
  )

# Pivot back to wide so each eid has columns per diagnosis (like d_15_cc previously)
d_15_future_cc <- d_long_future %>%
  select(eid, date_bsampling, diagnosis, case_15y_future) %>%
  pivot_wider(
    id_cols = c(eid, date_bsampling),
    names_from = diagnosis,
    values_from = case_15y_future
  )

# ---- 3) Existing data processing: join covariates and create medication flags etc. ----
# NOTE: I assume df, covs, a, job_vars, sleep, gen_covs, dep, vars_join, cohort_f, cohort exist in environment as before
data <- df %>%
  left_join(covs) %>%
  left_join(a) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  left_join(gen_covs) %>%
  left_join(dep) %>%
  left_join(vars_join %>% select(eid, chrono_Nightshift)) %>%
  # left_join(cohort_f) %>%   # keep if you previously used it
  mutate(has_prescription = ifelse(eid %in% cohort$eid, 1, 0),
         antihypertensive = ifelse(!is.na(C0), as.integer(C0 == TRUE), 0),
         sleep_medication = ifelse(!is.na(N05C), as.integer(N05C == TRUE), 0),
         antidepressants = ifelse(!is.na(N06), as.integer(N06 == TRUE), 0),
         mood_stabiliser = ifelse(!is.na(N03), as.integer(N03 == TRUE), 0),
         lithium = ifelse(!is.na(N05A), as.integer(N05A == TRUE), 0),
         across(c(has_prescription, antihypertensive, sleep_medication, antidepressants, mood_stabiliser, lithium), as.factor))

# residuals as before
data$res <- residuals(lm(pred_mean ~ time_day, data = data))

# ---- 4) Join future-case table into analysis dataset and create res_sd_group ----
data1 <- data %>%
  left_join(d_15_future_cc, by = c("eid", "date_bsampling")) %>%   # <-- replaced d_15_cc with d_15_future_cc
  mutate(
    res_sd_group = case_when(
      res >  2  ~ "Accelerated",
      res < -2  ~ "Delayed",
      TRUE        ~ "Middle"
    ),
    res_sd_group = factor(res_sd_group, levels = c("Middle", "Accelerated", "Delayed"))
  )


# ---- 6) Models ----
covariates <- c("sex", "age_recruitment", "chrono", "season", "TDI", "bmi", "smoking", "assessment_centre", paste0("PC", 1:20))

# Continuous res effect (per unit res)
res_continuous <- map_dfr(outcomes, function(outcome) {
  fit <- glm(
    as.formula(
      paste(
        # use the future-coded wide columns created earlier: these are named same as outcomes
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

# Extremes vs Middle
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

# ---- 7) Counts: number of incident cases (future 15y) among analysis sample ----
ns <- d_long_future %>%
  filter(eid %in% data1$eid) %>%
  rename(outcome = diagnosis) %>%
  group_by(outcome) %>%
  summarise(
    cases = sum(case_15y_future == 1, na.rm = TRUE),
    controls = sum(case_15y_future == 0, na.rm = TRUE),
    excluded = sum(is.na(case_15y_future)),
    .groups = "drop"
  )

# ---- 8) Prepare plotting dataframe (same logic as original, but now based on future cases) ----
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
    fill_sig = if_else(pfdr < 0.05, group, NA_character_)
  )

# Save results (changed filename to indicate future 15y)
saveRDS(plot_df, "data_share/associations_results_disease_CA_future15y.rds")


plot_df_inc <- readRDS("data_share/associations_results_disease_CA_future15y.rds")


# ---- 9) Plot (unchanged style) ----
p_qs_inc <- ggplot(plot_df_inc, aes(x = estimate, y = outcome_full, color = fct_rev(group),  fill = fill_sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey70") +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0,
    position = position_dodge(width = 0.65)
  ) +
  geom_point(
    shape = 21,
    size = 2,
    stroke = 1,
    position = position_dodge(width = 0.65)
  ) +
  geom_text(
    aes(x = conf.high * 1.15, label = count_label),
    size = 3,
    position = position_dodge(width = 0.65),
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  scale_color_manual(
    values = c(
      "Delayed" = "#b2182b",
      "Accelerated" = "#2166ac",
      "Continuous" = "black"
    ), label = c("Delayed (CA < -2)",
                 "Accelerated (CA > 2)",
                 "Continuous CA")
  ) +
  scale_fill_manual(
    values = c(
      "Delayed"     = "#b2182b",
      "Accelerated" = "#2166ac",
      "Continuous"  = "black"
    ),
    na.value = "white",
    guide = "none"
  ) +
  labs(title = "Incident cases (<15y)",
    x = "Odds Ratio",
    y = NULL,
    fill = "FDR < 5%",
    color = "Risk factor"
  ) +
  guides(
    color  = guide_legend(order = 1, nrow = 3, byrow = TRUE, reverse  = TRUE, override.aes = list(size = 4))
  ) +
  theme_classic(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 30, 10, 10),
    #legend.position = c(0.55, 0.85),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill  = "white",
      color = "black",
      linewidth = 0.4
    ),
    legend.position = "none"
  )

ggsave("plots/FX_disease_q_future15y.png", p_qs_inc, width = 6, height = 8)


install.packages("patchwork")

p_qs + p_qs_inc

