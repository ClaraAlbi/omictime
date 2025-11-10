library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("broom")
install.packages("table1")
install.packages("forcats")
library(broom)
library(forcats)
library(lubridate)

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  filter(smoking != "-3") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
         )

job_vars <- data.table::fread("/mnt/project/job_vars.tsv") %>%
  mutate(night_shift = case_when(`3426-0.0` == 1 ~ "Never",
                                 `3426-0.0` == 2 ~ "Sometimes",
                                 `3426-0.0` == 3 ~ "Usually",
                                 `3426-0.0` == 4 ~ "Always"),
         shift_work = case_when(`826-0.0` == 1 ~ "Never/rarely",
                                 `826-0.0` == 2 ~ "Sometimes",
                                 `826-0.0` == 3 ~ "Usually",
                                 `826-0.0` == 4 ~ "Always")) %>%
  #filter(`3426-0.0` %in% 1:4) %>%
  mutate(night_shift = factor(night_shift, levels = c("Never", "Sometimes", "Usually", "Always")),
         shift_work = factor(shift_work, levels = c("Never/rarely", "Sometimes", "Usually", "Always")))

pcs <- data.table::fread("/mnt/project/covariates.txt") %>%
  select(eid = 1, contains("PC"))

# MH <- data.table::fread("/mnt/project/psychosocial_MH.csv") %>%
#   select(eid, contains("p20"), contains("p19")) %>%
#   select(-contains("p2012")) %>%
#   mutate(across(contains("p"), as.factor))

dep <- data.table::fread("/mnt/project/other_covs.tsv") %>%
  select(eid, TDI = 2, eth = 3)


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
    wakeup = case_when(
      wakeup == 1 ~ "Not at all easy",
      wakeup == 2~ "Not very easy",
      wakeup == 3 ~ "Fairly easy",
      wakeup == 4 ~"Very easy",
      wakeup == -1 ~ NA_character_),
      wakeup = factor(wakeup, levels = c("Very easy", "Fairly easy",  "Not very easy", "Not at all easy")),
    chrono = factor(chrono, levels = c("Definitely morning", "Rather morning", "Don't know", "Rather evening", "Definitely evening")),
    ever_insomnia = case_when(ever_insomnia == 1 ~ "Never/rarely",
                              ever_insomnia == 2 ~ "Sometimes",
                              ever_insomnia == 3 ~ "Usually"),
    ever_insomnia = factor(ever_insomnia, levels = c("Never/rarely", "Sometimes", "Usually")))

df_temp <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  inner_join(readRDS("/mnt/project/olink_int_replication.rds") %>% select(-date_bsampling)) %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  mutate(date = as.POSIXct(date_bsampling, tz = "Europe/London"))

years <- unique(year(df_temp$date))

dst_transitions <- map_df(years, ~{
  test_dates <- seq(as.Date(paste0(.x, "-01-01")),
                    as.Date(paste0(.x, "-12-31")),
                    by = "day")

  is_dst <- dst(force_tz(as.POSIXct(test_dates), "Europe/London"))
  dst_changes <- which(diff(is_dst) != 0)

  if (length(dst_changes) >= 2) {
    tibble(
      year = .x,
      spring_dst = test_dates[dst_changes[1] + 1],
      fall_dst = test_dates[dst_changes[2] + 1]
    )
  }
})

df <- df_temp %>%
  mutate(
    date = as.POSIXct(date_bsampling, tz = "Europe/London"),
    is_dst = as.logical(as.POSIXlt(date, tz = "Europe/London")$isdst),
    is_dst = factor(is_dst, levels = c(FALSE, TRUE), labels = c("No", "Yes")),
    date_only = as.Date(date)
  ) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = FALSE) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  ungroup() %>%
  mutate(
    season = case_when(
      m %in% c("12", "01", "02") ~ "Winter",
      m %in% c("03", "04", "05") ~ "Spring",
      m %in% c("06", "07", "08") ~ "Summer",
      m %in% c("09", "10", "11") ~ "Fall"
    ),
    season = relevel(as.factor(season), ref = "Winter"),

    # Weekday/weekend classification
    day_of_week = wday(date_bsampling, label=TRUE, week_start=1),
    day_of_week = relevel(factor(day_of_week, ordered = FALSE), ref = "Wed"),
    is_weekend = wday(date_bsampling, week_start = 1) %in% c(6),
    day_type = if_else(is_weekend, "Weekend", "Weekday"),

    # Add year for joining
    year = year(date_bsampling)
  ) %>%
  left_join(dst_transitions, by = "year") %>%
  mutate(
    # DST classification
    dst_category = case_when(
      date_bsampling >= (spring_dst - 2) & date_bsampling <= (spring_dst - 1) ~ "before_spring_DST",
      date_bsampling >= (spring_dst + 1) & date_bsampling <= (spring_dst + 2) ~ "after_spring_DST",
      date_bsampling >= (fall_dst - 2) & date_bsampling <= (fall_dst - 1) ~ "before_fall_DST",
      date_bsampling >= (fall_dst + 1) & date_bsampling <= (fall_dst + 2) ~ "after_fall_DST",
      TRUE ~ "normal"
    ),
    dst_category = relevel(factor(dst_category), ref = "normal")
  ) #%>%
  select(-year, -spring_dst, -fall_dst, -date_only)


#labs <- readRDS("/mnt/project/biomarkers_3/covariate_res/res_labs.rds")

#fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

#colnames(labs) <- c("eid", fields %>% filter(field_id %in% as.numeric(colnames(labs))) %>% pull(title))

data <- df %>%
  left_join(covs) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  #left_join(labs) %>%
  #left_join(MH) %>%
  left_join(pcs) %>%
  left_join(dep) %>%
  filter(!is.na(chrono)) %>%
  mutate(h = round(time_day, 0)) %>%
  filter(age_recruitment > 39)

data$res <- residuals(lm(pred_mean ~ time_day, data = data))
data$res_q <- ntile(data$res, 5)


a <- broom::tidy(lm(res ~ dst_category, data = data))

library(ggplot2)


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


# --- 1. Predictor list ---
vars <- c("time_day", "age_recruitment", "sex", "chrono", "h_sleep", "ever_insomnia",
          "season", "night_shift", "smoking", "bmi", "is_dst", "wakeup", "shift_work", "TDI")

#vars <- colnames(MH)[-1]

covars <- c("sex", "age_recruitment", paste0("PC", 1:10))

#data <- data %>% filter(eth == 1001)

# --- 2. Fit models and extract results ---
results <- map_dfr(vars, function(v) {
  adj_vars <- if (v %in% c("sex", "age_recruitment")) character(0) else covars

  # Combine predictor + covariates safely
  rhs <- paste(c(v, adj_vars), collapse = " + ")
  f <- as.formula(paste("res ~", rhs))

  fit <- lm(f, data = data)

  tidy(fit) %>%
    filter(str_detect(term, paste0("^", v))) %>%
    mutate(predictor = v, reference = FALSE) %>%
    {
      if (is.factor(data[[v]])) {
        ref <- tibble(
          term = paste0(v, levels(data[[v]])[1]),
          estimate = 0, std.error = NA, statistic = NA, p.value = NA,
          predictor = v, reference = TRUE
        )
        bind_rows(ref, .)
      } else .
    }
})


results <- readRDS("data_share/results_associations_phenotypes.rds")


# --- 3. Compute ORs, CIs, and domain categories ---
res <- results %>%
  mutate(term = case_when(str_detect(term, "night_shift") ~ paste0(term, "."),
                          TRUE ~ term)) %>%
  mutate(
    OR    = exp(estimate),
    lower = ifelse(reference, 1, exp(estimate - 1.96 * std.error)),
    upper = ifelse(reference, 1, exp(estimate + 1.96 * std.error)),
    Category = case_when(
      predictor %in% c("age_recruitment", "sex", "bmi", "smoking", "TDI") ~ "Demographics",
      predictor %in% c("time_day") ~ "Time",
      predictor %in% c("chrono", "h_sleep", "ever_insomnia", "wakeup") ~ "Sleep",
      predictor %in% c("season", "is_dst") ~ "Season",
      predictor %in% c("shift_work", "night_shift") ~ "Job",
      TRUE ~ "Other"
    )
  )

#saveRDS(res, "data_share/results_associations_phenotypes.rds")

# --- 4. Factor level lookup for consistent ordering ---
factor_lookup <- map_dfr(vars, \(v)
                         tibble(predictor = v,
                                levels_list = list(if (is.factor(data[[v]])) levels(data[[v]]) else NULL))
)

# --- 5. Clean display terms ---
res2 <- res %>%
  mutate(level = str_remove(term, paste0("^", predictor)),
         display_term = ifelse(reference, paste0(level, " (ref)"),
                               ifelse(level == "", predictor, level))) %>%
  group_by(predictor) %>% mutate(n = n()) %>%
  mutate(display_term = case_when(n == 1 ~ "",
                                  TRUE ~ display_term)) %>%
  left_join(factor_lookup, by = "predictor") %>%
  group_by(predictor) %>%
  mutate(display_term = {
    lvls <- levels_list[[1]]
    if (is.null(lvls)) factor(display_term, levels = unique(display_term))
    else factor(display_term, levels = c(paste0(lvls[1], " (ref)"), lvls[-1]))
  }) %>%
  ungroup()

# --- 6. Pretty labels and domain colors ---
pretty_predictor <- c(
  time_day = "Time of Day", age_recruitment = "Age at Recruitment", sex = "Sex",
  chrono = "Chronotype", h_sleep = "Sleep Duration", ever_insomnia = "Insomnia",
  season = "Season", is_dst = "Daylight Savings", night_shift = "Night Shift",
  smoking = "Smoking", wakeup = "Waking Easiness", bmi = "BMI",
  shift_work = "Shift Work", TDI  = "Townsend DI"
)

domain_colors <- c(
  "Demographics" = "#d62728", "Job" = "#1f77b4", "Time" = "#2ca02c",
  "Season" = "#9467bd", "Sleep" = "#ff7f0e"
)

res_plot <- res2 %>%
  mutate(
    predictor_label = factor(pretty_predictor[predictor], levels = pretty_predictor),
    Domain = factor(Category, levels = c("Time", "Demographics", "Sleep", "Season", "Job"))
  )

# --- 7. Plot ---
p_d <- res_plot %>%
  filter(Domain == "Demographics") %>%
  ggplot(aes(x = OR, y = fct_rev(display_term), color = Domain)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(aes(shape = reference), size = 3, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21),
                     labels = c(`FALSE` = "Estimate", `TRUE` = "Reference")) +
  scale_color_manual(values = domain_colors) +
  scale_x_continuous(limits = c(0.5, 1.2)) +
  facet_grid(rows = vars(predictor_label), scales = "free", space = "free") +
  labs(x = "Odds Ratio (95% CI)", y = NULL, shape = NULL, color = "Domain") +
  theme_bw(base_size = 14) +
  theme(
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    panel.grid.major.y = element_blank(),
    strip.background = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

res_plot %>%
  filter(Domain %in% c("Sleep")) %>%
  ggplot(aes(x = OR, y = fct_rev(display_term), color = Domain)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(aes(shape = reference), size = 3, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21),
                     labels = c(`FALSE` = "Estimate", `TRUE` = "Reference")) +
  scale_color_manual(values = domain_colors) +
  scale_x_continuous(limits = c(0.5, 1.2)) +
  facet_grid(rows = vars(predictor_label), scales = "free", space = "free") +
  labs(x = "Odds Ratio (95% CI)", y = NULL, shape = NULL, color = "Domain") +
  theme_bw(base_size = 14) +
  theme(
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    panel.grid.major.y = element_blank(),
    strip.background = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

p
# ggsave("plots/FX_phenotypes.png", p, width = 12, height = 7)












b <- res %>%
  filter(predictor == "night_shift")

a <- d2 %>%
  filter(predictor_label == "Night Shift")
  ggplot(aes(x = OR, y = display_term, color = Domain)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(aes(shape = reference), size = 3, fill = "white")
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21),
                     labels = c(`FALSE` = "Estimate", `TRUE` = "Reference")) +
  scale_color_manual(values = domain_colors) +
  facet_grid(rows = vars(Domain, predictor_label), scales = "free", space = "free")


#### biomr


results <- map_dfr(colnames(labs)[-1], function(v) {
  f <- as.formula(paste0("res ~ ", "`",v, "`", " + time_day + age_recruitment + sex"))
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


res <- results %>%
  filter(term != "time_day") %>%
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
      predictor %in% colnames(labs) ~ "Labs",
      TRUE ~ "Other"
    )
  )

res %>%
  filter(Category == "Labs") %>%
  ggplot(aes(x = OR, y = fct_rev(term))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(aes(shape = reference), size = 3, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21),
                     labels = c(`FALSE` = "Estimate", `TRUE` = "Reference")) +
  #scale_color_manual(values = domain_colors) +
  #facet_grid(rows = vars(Domain), scales = "free", space = "free") +
  #facet_grid(~ Domain + predictor_label, scales = "free_y", ncol = 1, strip.position = "right") +
  labs(
    x = "Odds Ratio (95% CI)",
    y = NULL,
    shape = NULL,
    color = "Domain"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.text.y.right = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 16),
    legend.position = "bottom"
  )


ggplot(data, aes(x = age_recruitment, y = res)) +
  geom_smooth()
  geom_point()
