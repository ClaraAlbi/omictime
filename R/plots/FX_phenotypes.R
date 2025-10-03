library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("broom")
install.packages("table1")
install.packages("forcats")
library(broom)

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = as.factor(sex))

job_vars <- data.table::fread("/mnt/project/job_vars.tsv") %>%
  mutate(night_shift = case_when(`3426-0.0` == 1 ~ "Never",
                                 `3426-0.0` == 2 ~ "Sometimes",
                                 `3426-0.0` == 3 ~ "Usually",
                                 `3426-0.0` == 4 ~ "Always")) %>%
  filter(`3426-0.0` %in% 1:4) %>%
  mutate(night_shift = as.factor(night_shift),
         night_shift = relevel(night_shift, ref = "Never"))

pa <- data.table::fread("/mnt/project/sun_exposure.csv") %>%
  mutate(time_outdoors_s = case_when(as.integer(p1050_i0) > 0 & as.integer(p1050_i0) < 4  ~ "< 4h",
                                     as.integer(p1050_i0) > 4 ~ "> 4h"),
         time_outdoors_w = case_when(as.integer(p1060_i0) > 0 & as.integer(p1050_i0) < 2  ~ "< 2h",
                                     as.integer(p1060_i0) > 2 ~ "> 2h"))

sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid,
         h_sleep = `1160-0.0`,
         chrono = `1180-0.0`,
         ever_insomnia = `1200-0.0`,
         wakeup = `1170-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ "Definitely morning",
    chrono == 2 ~ "Rather morning",
    chrono == -1~ "Don't know",
    chrono == 3 ~ "Rather evening",
    chrono == 4 ~ "Definitely evening",
    TRUE ~ NA_character_),
    chrono = factor(chrono, levels = c("Definitely morning", "Rather morning", "Don't know", "Rather evening", "Definitely evening")),
    ever_insomnia = case_when(ever_insomnia == 1 ~ "Never/rarely",
                              ever_insomnia == 2 ~ "Sometimes",
                              ever_insomnia == 3 ~ "Usually"),
    ever_insomnia = factor(ever_insomnia, levels = c("Never/rarely", "Sometimes", "Usually")))

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  #filter(time_day > 12 & time_day < 18) %>%
  filter(i == 0) %>%
  #filter(cv %in% 1:5) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  ungroup() %>%
  mutate(season = case_when(m %in% c("12", "01", "02") ~ "Winter",
                            m %in% c("03", "04", "05") ~ "Spring",
                            m %in% c("06", "07", "08") ~ "Summer",
                            m %in% c("09", "10", "11") ~ "Fall"),
         season = relevel(as.factor(season), ref = "Winter"))

data <- df %>%
  left_join(covs) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  filter(!is.na(chrono)) %>%
  mutate(h = round(time_day, 0)) %>%
  filter(age_recruitment > 39)

data$res <- residuals(lm(pred_mean ~ time_day, data = data))
data$res_q <- ntile(data$res, 5)


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


tab_desc <- table1::table1(~ time_day + age_recruitment + factor(sex) + chrono + h_sleep + season + night_shift + ever_insomnia,
                           data = data,
                           render.cont = my_render_cont)



# Predictor list
vars <- c("time_day", "age_recruitment", "sex",
          "chrono", "h_sleep", "ever_insomnia",
          "season", "night_shift", "smoking", "bmi")

# Loop over predictors
results <- map_dfr(vars, function(v) {
  f <- as.formula(paste("abs(res) ~", v))
  fit <- lm(f, data = data)

  res <- tidy(fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(predictor = v, reference = FALSE)

  # Add reference row if factor
  if (is.factor(data[[v]])) {
    ref_level <- levels(data[[v]])[1]
    ref_row <- tibble(
      term = paste0(v, ref_level),
      estimate = 0, std.error = NA, statistic = NA, p.value = NA,
      predictor = v, reference = TRUE
    )
    res <- bind_rows(ref_row, res)
  }
  res
})

# Add ORs, CIs, and categories
res <- results %>%
  mutate(
    OR    = exp(estimate),
    lower = ifelse(reference, 1, exp(estimate - 1.96*std.error)),
    upper = ifelse(reference, 1, exp(estimate + 1.96*std.error)),
    Category = case_when(
      predictor %in% c("age_recruitment", "sex", "bmi", "smoking") ~ "Demographics",
      predictor %in% c("chrono", "h_sleep", "ever_insomnia") ~ "Sleep",
      predictor == "season" ~ "Season",
      predictor == "night_shift" ~ "Job",
      predictor == "time_day" ~ "Other",
      TRUE ~ "Other"
    )
  )


factor_lookup <- map_dfr(vars, function(v) {
  if (is.factor(data[[v]])) {
    tibble(predictor = v,
           levels_list = list(levels(data[[v]])))
  } else {
    tibble(predictor = v,
           levels_list = list(NULL))
  }
})

# --- add labels first ---
res <- res %>%
  mutate(
    level = gsub(predictor, "", term),
    level = ifelse(reference, paste0(level, " (ref)"), level),
    display_term = paste0(predictor, ": ", level)
  ) %>%
  left_join(factor_lookup, by = "predictor")

# --- reorder display_term using lookup ---
res <- res %>%
  rowwise() %>%
  mutate(display_term = if (!is.null(levels_list)) {
    lvls <- levels_list
    factor(display_term,
           levels = c(paste0(predictor, ": ", lvls[1], " (ref)"),
                      paste0(predictor, ": ", lvls[-1])))
  } else {
    factor(display_term, levels = unique(display_term))
  }) %>%
  ungroup()

# Pretty labels for predictors
pretty_predictor <- c(
  time_day = "Time of Day",
  age_recruitment = "Age at Recruitment",
  sex = "Sex",
  chrono = "Chronotype",
  h_sleep = "Sleep Duration",
  ever_insomnia = "Insomnia",
  season = "Season",
  night_shift = "Night Shift",
  smoking = "Smoking",
  bmi = "BMI"
)

  res_plot <- res %>%
  mutate(
    predictor_label = pretty_predictor[predictor],
    predictor_label = factor(predictor_label, levels = pretty_predictor)
  )

res_plot <- res %>%
  mutate(
    level_label = gsub(paste0("^", predictor), "", term),
    level_label = ifelse(reference, paste0(level_label, " (ref)"), level_label),
    predictor_label = pretty_predictor[predictor],
    predictor_label = factor(predictor_label, levels = pretty_predictor)
  )

# Plot
ggplot(res_plot,
       aes(x = fct_rev(level_label), y = OR,
           color = Category, shape = reference)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  facet_grid(rows = vars(Category, predictor_label),
             scales = "free_y", space = "free_y") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 19),
                     labels = c("FALSE" = "Estimate", "TRUE" = "Reference")) +
  coord_flip() +
  theme_bw(base_size = 20) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = "Odds Ratio (95% CI)",
       color = "Domain", shape = "")


