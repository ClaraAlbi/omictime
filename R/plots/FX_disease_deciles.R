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

pcs <- data.table::fread("/mnt/project/covariates.txt") %>%
  select(eid = 1, contains("PC"))

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
  left_join(pcs) %>%
  left_join(chrono)

d <- bind_cols(fread("/mnt/project/top_diseases_IEFG.csv"), fread("/mnt/project/immune.csv") %>% select(-eid)) %>%
  mutate(across(-eid, as_date))

data1 <- data %>%
  select(eid, res, sex, age_recruitment,assessment_centre, smoking, bmi, chrono, any_of(paste0("PC", 1:20))) %>%
  left_join(d
            #%>% select(eid, p130892, p130874, p130894, p130896)
            ) %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds") %>% select(eid, date_bsampling)) %>%
  mutate(
    across(
      contains("p13"),
      ~ if_else(
        !is.na(na_if(as.character(.x), "")) &
          as.Date(na_if(as.character(.x), "")) <= date_bsampling,
        1L,
        0L
      ),
      .names = "{.col}_prevalent"
    )
  ) %>%
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


fields <- data.table::fread("/mnt/project/Showcase metadata//field.tsv")

res_all <- bind_rows(res_top_vs_mid, res_bottom_vs_mid) %>%
  left_join(case_counts %>% mutate(term = case_when(ca_extreme == "Bottom10" ~ "ca_bottom",
                                                    ca_extreme == "Top10" ~ "ca_top"))) %>%
  filter(n_cases > 10) %>%
  mutate(
    conf.low  = exp(log(estimate) - 1.96 * std.error),
    conf.high = exp(log(estimate) + 1.96 * std.error),
    outcome = str_remove(outcome, "p"),
    field_id = as.numeric(str_remove(outcome, "_prevalent"))) %>%
  left_join(fields %>% select(field_id, title)) %>%
  mutate(family = str_extract(title, "(?<=Date )\\S+(?= first reported)"),
         family = str_sub(family, 1, 2),
         disorder = sub(".*\\((.*)\\).*", "\\1", title),
         disorder = str_remove(disorder, "\\)"),
         disorder = str_to_sentence(disorder)) #%>%
#filter(!field_id  %in% c(130898, 130902, 130932, 130944, 130852)) %>%
filter(!family %in% c("F5", "F6", "F1", "FO" ))

res_all

res_all %>%
  #filter(field_id == 130892) %>%
  mutate(disorder = fct_reorder(disorder, abs(estimate))) %>%
  ggplot(aes(x = estimate, y = disorder, color = contrast)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.2
  ) +
  geom_text(aes(label = n_cases)) +
  scale_x_log10() +
  facet_wrap(~family, scales = "free") +
  labs(
    x = "Odds Ratio (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 8)
