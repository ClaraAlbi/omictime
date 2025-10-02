library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("broom")
install.packages("table1")
install.packages("forcats")

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
                              ever_insomnia == 3 ~ "Usually"))

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
  left_join(pa) %>%
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

vars <- c("time_day", "age_recruitment", "sex",
          "chrono", "h_sleep", "ever_insomnia",
          "season", "night_shift", "smoking", "bmi")

results <- map_dfr(vars, function(v) {
  f <- as.formula(paste("abs(res) ~", v))
  fit <- lm(f, data = data)

  res <- broom::tidy(fit) %>%
    filter(term != "(Intercept)") %>%
    mutate(predictor = v, reference = FALSE)

  # Add reference row if predictor is factor
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

res <- res %>%
  mutate(
    # Strip out the variable name prefix from "term"
    level = ifelse(reference,
                   paste0("Reference: ", gsub(predictor, "", term)),
                   gsub(predictor, "", term)),
    display_term = paste0(predictor, ": ", level)
  )

# Now order within each predictor:
res <- res %>%
  group_by(predictor) %>%
  mutate(display_term = factor(display_term,
                               levels = unique(display_term))) %>%
  ungroup()


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
    # remove predictor prefix from term
    level_label = gsub(paste0("^", predictor), "", term),
    # mark reference rows
    level_label = ifelse(reference, paste0(level_label, " (ref)"), level_label),
    # pretty predictor labels
    predictor_label = pretty_predictor[predictor],
    predictor_label = factor(predictor_label, levels = pretty_predictor)
  )


library(forcats)

ggplot(res_plot,
       aes(x = fct_rev(level_label), y = OR,
           color = Category, shape = reference)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, na.rm = TRUE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  facet_grid(rows = vars(Category, predictor_label),
             scales = "free_y", space = "free_y") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c("TRUE" = 21, "FALSE" = 19),
                     labels = c("FALSE" = "Estimate", "TRUE" = "Reference")) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.y = element_text(angle = 0, hjust = 0),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = "Odds Ratio (95% CI)",
       color = "Domain", shape = "")
# Run regressions separately
results <- map_dfr(vars, function(v) {
  f <- as.formula(paste("res ~", v))
  fit <- lm(f, data = data)
  broom::tidy(fit) %>% filter(term != "(Intercept)") %>% mutate(predictor = v)
})

res <- results %>%
  mutate(
    OR = exp(estimate),
    lower = exp(estimate - 1.96*std.error),
    upper = exp(estimate + 1.96*std.error),
    Category = case_when(
      term %in% c("age_recruitment", "sex", "bmi", "smoking") ~ "Demographics",
      grepl("^chrono|h_sleep|ever_insomnia", term) ~ "Sleep",
      grepl("^season", term) ~ "Season",
      grepl("^night_shift", term) ~ "Job",
      term == "time_day" ~ "Other",
      TRUE ~ "Other"
    )
  )

ggplot(res,
       aes(x = reorder(term, OR), y = OR, color = Category)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  facet_grid(rows = vars(Category, predictor) , scales = "free", space = "free") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(title = "Odds Ratios by Predictor",
       x = "Predictor",
       y = "Odds Ratio (95% CI)",
       color = "Category")



# Step 1: rescale predictions to match min/max
true_min <- min(a$h, na.rm = TRUE)
true_max <- max(a$h, na.rm = TRUE)
pred_min <- min(a$pred_mean, na.rm = TRUE)
pred_max <- max(a$pred_mean, na.rm = TRUE)

scale_factor <- (true_max - true_min) / (pred_max - pred_min)
shift <- true_min - pred_min * scale_factor

a <- a %>%
  mutate(pred_rescaled = pred_mean * scale_factor + shift)

# Step 2: spline calibration on rescaled predictions
fit_spline <- lm(h ~ ns(pred_rescaled, df = 4), data = a)
a <- a %>%
  ungroup() %>%
  mutate(pred_calib_spline = predict(fit_spline, newdata = a))

calib_df <- a %>%
  group_by(h) %>%
  summarise(
    raw      = mean(pred_mean, na.rm = TRUE),
    spline   = mean(pred_calib_spline, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(raw, spline),
                      names_to = "method", values_to = "mean_pred") %>%
  mutate(method = recode(method,
                         raw = "Raw",
                         spline = "Spline-calibrated"))

ggplot(calib_df, aes(x = h, y = mean_pred, color = method, linetype = method)) +
  geom_point(aes(size = n), alpha = 0.6) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  coord_equal(xlim = c(8,21), ylim = c(8,21)) +
  theme_classic(base_size = 14) +
  labs(x = "Recorded time", y = "Mean predicted time",
       color = "Method", linetype = "Method", size = "Sample size")




