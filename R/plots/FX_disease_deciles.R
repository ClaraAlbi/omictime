library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)
library(stringr)
library(lubridate)
library(data.table)
install.packages("broom")
install.packages("forcats")
library(broom)
library(forcats)

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
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T)
df$res <- residuals(lm(pred_mean ~ time_day, data = df))

data <- df %>%
  left_join(covs) %>%
  left_join(gen_covs) %>%
  left_join(chrono)



data1 <- data %>%
  select(eid, res, sex, age_recruitment,assessment_centre, smoking, bmi, chrono, any_of(paste0("PC", 1:20))) %>%
  left_join(d) %>%
  mutate(
    across(29:40,
           ~ if_else(
             !is.na(na_if(as.character(.x), "")) &
               as.Date(na_if(as.character(.x), "")) <= date_bsampling,
             1L,
             0L
           ),
           .names = "{.col}_prevalent"
    )
  ) 



data1 <- data1 %>%
  mutate(
    ca_quintile = ntile(res, 5),
    ca_quintile_lab = paste0("Q", ca_quintile)
  )

case_counts_vs_mid <- data1 %>%
  mutate(
    contrast = case_when(
      ca_quintile == 1 ~ "Q1_vs_Q3",
      ca_quintile == 2 ~ "Q2_vs_Q3",
      ca_quintile == 3 ~ "Q3",
      ca_quintile == 4 ~ "Q4_vs_Q3",
      ca_quintile == 5 ~ "Q5_vs_Q3"
    )
  ) %>%
  select(contrast, all_of(outcomes)) %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "outcome",
    values_to = "case"
  ) %>%
  group_by(outcome, contrast) %>%
  summarise(
    n_cases = sum(case == 1, na.rm = TRUE),
    .groups = "drop"
  )

case_counts_total <- data1 %>%
  select(all_of(outcomes)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "outcome",
    values_to = "case"
  ) %>%
  group_by(outcome) %>%
  summarise(
    n_cases = sum(case == 1, na.rm = TRUE)
  ) %>%
  mutate(contrast = "CA_per_SD")

middle_quintile <- 3

outcomes <- names(data1) %>%
  grep("_prevalent$", ., value = TRUE)


covariates <- c(
  "sex",
  "age_recruitment",
  "bmi",
  "smoking",
  "assessment_centre",
  paste0("PC", 1:20)
)

quintiles <- 1:5

res_quintile_vs_mid <- map_dfr(quintiles, function(q) {
  
  if (q == middle_quintile) return(NULL)
  
  map_dfr(outcomes, function(outcome) {
    
    df <- data1 %>%
      filter(ca_quintile %in% c(q, middle_quintile)) %>%
      mutate(q_case = if_else(ca_quintile == q, 1L, 0L))
    
    fit <- run_logistic(
      as.formula(
        paste(outcome, "~ q_case +", paste(covariates, collapse = " + "))
      ),
      df
    )
    
    if (is.null(fit)) return(NULL)
    
    tidy(fit, exponentiate = TRUE) %>%
      filter(term == "q_case") %>%
      mutate(
        outcome  = outcome,
        contrast = paste0("Q", q, "_vs_Q3"),
        quintile = q
      )
  })
}) %>%
  mutate(
    conf.low  = exp(log(estimate) - 1.96 * std.error),
    conf.high = exp(log(estimate) + 1.96 * std.error)
  )


data1 <- data1 %>%
  mutate(res_z = as.numeric(scale(res)))

res_quant <- map_dfr(outcomes, function(outcome) {
  
  fit <- run_logistic(
    as.formula(
      paste(outcome, "~ res_z +", paste(covariates, collapse = " + "))
    ),
    data1
  )
  
  if (is.null(fit)) return(NULL)
  
  tidy(fit, exponentiate = TRUE) %>%
    filter(term == "res_z") %>%
    mutate(
      outcome  = outcome,
      contrast = "CA_per_SD",
      quintile = NA_integer_
    )
}) %>%
  mutate(
    conf.low  = exp(log(estimate) - 1.96 * std.error),
    conf.high = exp(log(estimate) + 1.96 * std.error)
  )

res_all <- bind_rows(res_quintile_vs_mid, res_quant)

res_all <- res_all %>%
  left_join(case_counts_vs_mid, by = c("outcome", "contrast")) %>%
  left_join(
    case_counts_total,
    by = c("outcome", "contrast"),
    suffix = c("", "_total")
  ) %>%
  mutate(
    n_cases = coalesce(n_cases, n_cases_total)
  ) %>%
  select(-n_cases_total)

plot_df <- res_all %>%
  mutate(
    count_label = if_else(
      !is.na(n_cases),
      paste0("n=", n_cases),
      NA_character_
    )
  )

res_all %>%
  mutate(
    count_label = paste0("n=", n_cases),
    outcome = fct_reorder(outcome, estimate),
    contrast = factor(
      contrast,
      levels = c("Q1_vs_Q3", "Q2_vs_Q3", "Q4_vs_Q3", "Q5_vs_Q3", "CA_per_SD")
    )
  ) %>%
  ggplot(aes(x = estimate, y = outcome, color = contrast)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.25,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    size = 2.2,
    position = position_dodge(width = 0.6)
  ) +
  geom_text(
    aes(
      x = conf.high * 1.15,
      label = count_label
    ),
    size = 3,
    position = position_dodge(width = 0.6),
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    x = "Odds Ratio (log scale)",
    y = NULL,
    color = NULL
  ) +
  scale_color_manual(
    values = c(
      "Q1_vs_Q3"  = "#2166ac",  # blue
      "Q2_vs_Q3"  = "#67a9cf",  # light blue
      "Q4_vs_Q3"  = "#f4a582",  # light red
      "Q5_vs_Q3"  = "#b2182b",  # red
      "CA_per_SD" = "black"     # continuous → clearly distinct
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 14)
  )

d <- readRDS("diseases_circadian.rds") %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds") %>% select(eid, date_bsampling))


data1 <- data %>%
  select(eid, res, sex, age_recruitment,assessment_centre, smoking, bmi, chrono, any_of(paste0("PC", 1:20))) %>%
  left_join(d) %>%
  mutate(
    across(29:40,
           ~ if_else(
             !is.na(na_if(as.character(.x), "")) &
          as.Date(na_if(as.character(.x), "")) <= date_bsampling,
        1L,
        0L
      ),
      .names = "{.col}_prevalent"
    )
  ) 



  mutate(
    ca_decile = ntile(res, 10),
    ca_top10 = if_else(ca_decile == 10, 1L, 0L),
    ca_bottom10 = if_else(ca_decile == 1, 1L, 0L),
    ca_extreme = case_when(
      ca_decile == 10 ~ "Top10",
      ca_decile == 1  ~ "Bottom10",
      TRUE            ~ "Middle"
    )
  )

data_extremes <- data1 %>%
  filter(ca_decile %in% c(1, 10)) %>%
  mutate(
    ca_group = if_else(ca_decile == 10, "Top10", "Bottom10")
  )

q10 <- quantile(data1$res, probs = 0.10, na.rm = TRUE)
q90 <- quantile(data1$res, probs = 0.90, na.rm = TRUE)

data1 <- data1 %>%
  mutate(
    ca_extreme = case_when(
      res <= q10 ~ "Bottom10",
      res >= q90 ~ "Top10",
      TRUE       ~ "Middle"
    ),
    ca_top10 = if_else(res >= q90, 1L, 0L),
    ca_bottom10 = if_else(res <= q10, 1L, 0L)
  )



outcomes <- names(data1) %>%
  grep("_prevalent$", ., value = TRUE)


case_counts <- data1 %>%
  select(ca_extreme, all_of(outcomes)) %>%
  pivot_longer(
    cols = all_of(outcomes),
    names_to = "outcome",
    values_to = "case"
  ) %>%
  group_by(outcome, ca_extreme) %>%
  summarise(
    n_cases = sum(case == 1, na.rm = TRUE),
    .groups = "drop"
  )


covariates <- c(
  "sex",
  "age_recruitment",
  "bmi",
  "smoking",
#  "chrono",
  paste0("PC", 1:20)
)

run_logistic <- function(formula, data) {
  tryCatch(
    glm(formula, data = data, family = binomial),
    error = function(e) NULL
  )
}

res_top_vs_mid <- map_dfr(outcomes, function(outcome) {

  df <- data1 %>%
    filter(ca_extreme %in% c("Top10", "Middle")) %>%
    mutate(ca_top = if_else(ca_extreme == "Top10", 1L, 0L))

  fit <- run_logistic(
    as.formula(
      paste(outcome, "~ ca_top +", paste(covariates, collapse = " + "))
    ),
    df
  )

  if (is.null(fit)) return(NULL)

  tidy(fit, conf.int = FALSE, exponentiate = TRUE) %>%
    filter(term == "ca_top") %>%
    mutate(
      outcome = outcome,
      contrast = "Top10_vs_Middle"
    )
})

res_bottom_vs_mid <- map_dfr(outcomes, function(outcome) {

  df <- data1 %>%
    filter(ca_extreme %in% c("Bottom10", "Middle")) %>%
    mutate(ca_bottom = if_else(ca_extreme == "Bottom10", 1L, 0L))

  fit <- run_logistic(
    as.formula(
      paste(outcome, "~ ca_bottom +", paste(covariates, collapse = " + "))
    ),
    df
  )

  if (is.null(fit)) return(NULL)

  tidy(fit, conf.int = FALSE, exponentiate = TRUE) %>%
    filter(term == "ca_bottom") %>%
    mutate(
      outcome = outcome,
      contrast = "Bottom10_vs_Middle"
    )
})



res_all <- bind_rows(res_top_vs_mid, res_bottom_vs_mid) %>%
  left_join(
    case_counts %>%
      mutate(term = case_when(
        ca_extreme == "Bottom10" ~ "ca_bottom",
        ca_extreme == "Top10"    ~ "ca_top"
      )),
    by = c("outcome", "term")
  ) %>%
  filter(n_cases > 10) %>%
  mutate(
    conf.low  = exp(log(estimate) - 1.96 * std.error),
    conf.high = exp(log(estimate) + 1.96 * std.error)
  )



res_all %>%
  mutate(
    outcome = forcats::fct_reorder(outcome, estimate)
  ) %>%
  ggplot(aes(x = estimate, y = outcome, color = contrast)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.25,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    size = 2.2,
    position = position_dodge(width = 0.6)
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.5, 2),
    labels = label_number(accuracy = 0.01)
  ) +
  labs(
    x = "Odds Ratio (log scale)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    axis.text.y = element_text(size = 12)
  )
