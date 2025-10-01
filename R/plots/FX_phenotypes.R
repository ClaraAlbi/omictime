library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("broom")
install.packages("table1")

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex))

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

a <- df %>%  left_join(covs) %>%
  left_join(job_vars) %>%
  left_join(pa) %>%
  left_join(sleep) %>%
  filter(!is.na(chrono)) %>%
  mutate(h = round(time_day, 0)) %>%
  filter(age_recruitment > 39)

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



install.packages("table1")
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

a$res <- residuals(lm(pred_mean ~ time_day, data = a))
a$res_q <- ntile(a$res, 5)

tab_desc <- table1::table1(~ time_day + age_recruitment + factor(sex) + chrono + h_sleep + season + night_shift + ever_insomnia,
                           data = a,
                           render.cont = my_render_cont)


vars <- c("time_day", "age_recruitment", "sex",
          "chrono", "h_sleep", "ever_insomnia",
          "season", "night_shift", "smoking", "bmi")

# Run regressions separately
results <- map_dfr(vars, function(v) {
  f <- as.formula(paste("res ~", v))
  fit <- lm(f, data = a)
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

ggplot(results %>% filter(term != "(Intercept)"),
       aes(x = reorder(term, OR), y = OR, color = Category)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(title = "Odds Ratios by Predictor",
       x = "Predictor",
       y = "Odds Ratio (95% CI)",
       color = "Category")

t.test(res ~ factor(sex), data = a)

ggplot(a, aes(x = factor(h), y = pred_mean, fill = factor(chrono))) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Recorded time", y = "Mean prediction", fill = "Chronotype")

ggplot(a, aes(x = sex, y = pred_mean, fill = factor(sex))) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Recorded time", y = "Mean prediction", fill = "Sex")

ggplot(a, aes(x = factor(age_recruitment), y = pred_mean)) +
  geom_boxplot() +
  geom_hline(yintercept = mean(a$pred_mean), linetype = 2, color = "red") +
  theme_classic() +
  labs(x = "Recorded time", y = "Mean prediction", fill = "Chronotype")





