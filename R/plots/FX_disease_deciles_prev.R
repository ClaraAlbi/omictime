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
  filter(date_bsampling != "1900-01-01") %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "diagnosis",
    values_to = "earliest_dx_date"
  )

d_long <- d_long %>%
  left_join(covs %>% select(eid, birth_year)) %>%
  mutate(age = year(earliest_dx_date) - birth_year) %>%
  filter(is.na(age) | age >= 10) %>%
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
  select(eid, date_bsampling, diagnosis, case_15y, age) %>%
  pivot_wider(
    id_cols = c(eid, date_bsampling, age),
    names_from = diagnosis,
    values_from = case_15y
  )



##############


data <- df %>%
  left_join(covs) %>%
  left_join(a) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  left_join(gen_covs) %>%
  left_join(dep) %>%
  left_join(vars_join %>% select(eid, chrono_Nightshift)) #%>%
  left_join(cohort_f) %>%
  mutate(has_prescription = ifelse(eid %in% cohort$eid, 1, 0),
         antihypertensive = ifelse(!is.na(C0), as.integer(C0 == TRUE), 0),
         sleep_medication = ifelse(!is.na(N05C), as.integer(N05C == TRUE), 0),
         antidepressants = ifelse(!is.na(N06), as.integer(N06 == TRUE), 0),
         mood_stabiliser = ifelse(!is.na(N03), as.integer(N03 == TRUE), 0),
         lithium = ifelse(!is.na(N05A), as.integer(N05A == TRUE), 0),
         across(c(has_prescription, antihypertensive, sleep_medication, antidepressants, mood_stabiliser, lithium), as.factor))


data$res <- residuals(lm(pred_mean ~ time_day, data = data))


data1 <- data %>%
  #select(eid, res, sex, age_recruitment,assessment_centre,season, smoking, bmi, chrono, any_of(paste0("PC", 1:20))) %>%
  left_join(d_15_cc) %>%
  mutate(
    res_sd_group = case_when(
      res >  2  ~ "Accelerated",
      res < -2  ~ "Delayed",
      TRUE        ~ "Middle"
    ),
    res_sd_group = factor(res_sd_group, levels = c("Middle", "Accelerated", "Delayed"))
  )


library(table1)
my_render_cont <- function(x){
  with(
    stats.apply.rounding(stats.default(x)),
    c(
      "",
      `Mean (SD)` = sprintf("%s (%s)", MEAN, SD),
      `Median [Q1, Q3]` = sprintf("%s [%s, %s]", MEDIAN, Q1, Q3),
      `Min, Max` = sprintf("%s, %s", MIN, MAX)
    )
  )
}

pvalue_AM <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(names(x), times = sapply(x, length)))
  keep <- g %in% c("Middle", "Accelerated")
  y <- y[keep]
  g <- droplevels(g[keep])

  if (is.numeric(y)) {
    p <- tryCatch(
      t.test(y ~ g)$p.value,
      error = function(e) wilcox.test(y ~ g)$p.value
    )
  } else {
    tbl <- table(y, g)
    p <- if (any(tbl < 5)) fisher.test(tbl)$p.value else chisq.test(tbl)$p.value
  }

  format.pval(p, digits = 3, eps = 0.001)
}

pvalue_DM <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(names(x), times = sapply(x, length)))
  keep <- g %in% c("Middle", "Delayed")
  y <- y[keep]
  g <- droplevels(g[keep])

  if (is.numeric(y)) {
    p <- tryCatch(
      t.test(y ~ g)$p.value,
      error = function(e) wilcox.test(y ~ g)$p.value
    )
  } else {
    tbl <- table(y, g)
    p <- if (any(tbl < 5)) fisher.test(tbl)$p.value else chisq.test(tbl)$p.value
  }

  format.pval(p, digits = 3, eps = 0.001)
}


tab_desc <- table1::table1(
  ~ age_recruitment + sex + p30079 + TDI + bmi + smoking + season +
    day_type + fri_sun + autumnDST + springDST + chrono +
    h_sleep + wakeup + ever_insomnia + shift_work + night_shift +
    chrono_Nightshift #+ has_prescription + antihypertensive + leep_medication + antidepressants + mood_stabiliser + lithium
  | res_sd_group,
  data = data1,
  overall = FALSE,
  render.cont = my_render_cont,
  topclass = "Rtable1-grid",
  extra.col = list(
    `P (Acc vs Mid)` = pvalue_AM,
    `P (Del vs Mid)` = pvalue_DM
  )
)




covariates <- c("sex", "age_recruitment", "chrono", "season", "TDI", "bmi", "smoking", "assessment_centre", paste0("PC", 1:20))

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
    fill_sig = if_else(pfdr < 0.05, group, NA_character_))

saveRDS(plot_df, "data_share/associations_results_disease_CA.rds")


plot_df <- readRDS("data_share/associations_results_disease_CA.rds")


p_qs <- ggplot(plot_df, aes(x = estimate, y = outcome_full, color = fct_rev(group),  fill = fill_sig)) +
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
    aes(x = conf.high * 1.3, label = count_label),
    size = 3,
    position = position_dodge(width = 0.65),
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2, 4),
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
    na.value = "white",   # ← makes non-significant points hollow
    guide = "none"
  ) +
  labs(title = "Prevalent cases (<15y)",
    x = "Odds Ratio",
    y = NULL,
    fill = "FDR < 5%",
    color = "Circadian acceleration"
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
    legend.position = c(0.55, 0.85),
    legend.justification = c(0, 0),
    legend.background = element_rect(
      fill  = "white",
      color = "black",
      linewidth = 0.4
    )
  )

ggsave("plots/FX_disease_q_prev.png", p_qs, width = 6, height = 8)



