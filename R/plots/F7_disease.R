#remotes::install_github("lcpilling/ukbrapR")
library(ukbrapR)
library(tidyr)
library(dplyr)
install.packages("broom")
library(purrr)

ckd <- ukbrapR:::codes_df_ckd

ukbrapR:::codes_df_test

bp <- tribble(~condition, ~vocab_id, ~code,
              # "bp", "ICD10","F31.2",
              # "bp", "ICD9", "296",
              # "bp","Read2", "E11..11",
              # "bp", "CTV3","E1176",
              "sleep", "ICD10", "G47",
              "insomnias", "ICD10", "G470",
              "hypersomnias", "ICD10", "G471",
              "sleep-wake schedule", "ICD10", "G472",
              "apnea", "ICD10", "G473",
              "narcolepsy", "ICD10", "G474",
              "other", "ICD10" , "G478",
              "unspecified", "ICD10", "G479",
              "sleep", "ukb_noncancer", "1123")


# G47.0 Disorders of initiating and maintaining sleep [insomnias]113
# G47.1 Disorders of excessive somnolence [hypersomnias]137
# G47.2 Disorders of the sleep-wake schedule35
# G47.3 Sleep apnoea3847
# G47.4 Narcolepsy and cataplexy24
# G47.8 Other sleep disorders140
# G47.9 Sleep disorder, unspecified

diagnosis_list <- get_diagnoses(bp)
diagnosis_df <- get_df(diagnosis_list, group_by="condition")


df <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  filter(i == 0) %>%
  left_join(diagnosis_df)

df$res <- residuals(lm(pred_mean ~ time_day, data = df))

a <- df %>%
  pivot_longer(ends_with("bin")) %>%
  group_by(name) %>% nest() %>%
  mutate(m = map(data, ~glm(value ~ res, data = .x, family = binomial)),
         p = map(m, broom::tidy),
         n = map_dbl(data, ~sum(.x$value == 1, na.rm = T))) %>%
  unnest(p) %>% select(-data, -m)


data2 <- data %>%
  left_join(diagnosis_df) %>%
  mutate(
  ca_decile = ntile(res, 5),
  ca_top10 = if_else(ca_decile == 5, 1L, 0L),
  ca_bottom10 = if_else(ca_decile == 1, 1L, 0L),
  ca_extreme = case_when(
    ca_decile == 5 ~ "Top25",
    ca_decile == 1  ~ "Bottom25",
    TRUE            ~ "Middle"
  )
)
data_extremes <- data2 %>%
  filter(ca_decile %in% c(1, 5)) %>%
  mutate(
    ca_group = if_else(ca_decile == 5, "Top25", "Bottom25")
  )

q10 <- quantile(data2$res, probs = 0.25, na.rm = TRUE)
q90 <- quantile(data2$res, probs = 0.85, na.rm = TRUE)

data2 <- data2 %>%
  mutate(
    ca_extreme = case_when(
      res <= q10 ~ "Bottom25",
      res >= q90 ~ "Top25",
      TRUE       ~ "Middle"
    ),
    ca_top25 = if_else(res >= q90, 1L, 0L),
    ca_bottom25 = if_else(res <= q10, 1L, 0L)
  )



outcomes <- names(data2) %>%
  grep("_prev$", ., value = TRUE)


case_counts <- data2 %>%
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
  "chrono",
  paste0("PC", 1:20)
)

run_logistic <- function(formula, data) {
  tryCatch(
    glm(formula, data = data, family = binomial),
    error = function(e) NULL
  )
}

res_top_vs_mid <- map_dfr(outcomes, function(outcome) {

  df <- data2 %>%
    filter(ca_extreme %in% c("Top25", "Middle")) %>%
    mutate(ca_top = if_else(ca_extreme == "Top25", 1L, 0L))

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
      contrast = "Top25_vs_Middle"
    )
})

res_bottom_vs_mid <- map_dfr(outcomes, function(outcome) {

  df <- data2 %>%
    filter(ca_extreme %in% c("Bottom25", "Middle")) %>%
    mutate(ca_bottom = if_else(ca_extreme == "Bottom25", 1L, 0L))

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
      contrast = "Bottom25_vs_Middle"
    )
})



res_all <- bind_rows(res_top_vs_mid, res_bottom_vs_mid) %>%
  left_join(case_counts %>% mutate(term = case_when(ca_extreme == "Bottom25" ~ "ca_bottom",
                                                    ca_extreme == "Top25" ~ "ca_top"))) %>%
  #filter(n_cases > 10) %>%
  mutate(
    conf.low  = exp(log(estimate) - 1.96 * std.error),
    conf.high = exp(log(estimate) + 1.96 * std.error))

