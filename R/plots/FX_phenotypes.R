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
                                 `3426-0.0` == 4 ~ "Always")) %>%
  #filter(`3426-0.0` %in% 1:4) %>%
  mutate(night_shift = factor(night_shift, levels = c("Never", "Sometimes", "Usually", "Always")))

# pa <- data.table::fread("/mnt/project/sun_exposure.csv") %>%
#   mutate(time_outdoors_s = case_when(as.integer(p1050_i0) > 0 & as.integer(p1050_i0) < 4  ~ "< 4h",
#                                      as.integer(p1050_i0) > 4 ~ "> 4h"),
#          time_outdoors_w = case_when(as.integer(p1060_i0) > 0 & as.integer(p1050_i0) < 2  ~ "< 2h",
#                                      as.integer(p1060_i0) > 2 ~ "> 2h"))

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

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  #filter(time_day > 12 & time_day < 18) %>%
  filter(i == 0) %>%
  #filter(cv %in% 1:5) %>%
  mutate(date = as.POSIXct(
    date_bsampling,
    tz = "Europe/London"),
    is_dst = as.logical(as.POSIXlt(date, tz = "Europe/London")$isdst),
    is_dst = factor(is_dst, levels = c(FALSE, TRUE), labels = c("No", "Yes"))) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  ungroup() %>%
  mutate(season = case_when(m %in% c("12", "01", "02") ~ "Winter",
                            m %in% c("03", "04", "05") ~ "Spring",
                            m %in% c("06", "07", "08") ~ "Summer",
                            m %in% c("09", "10", "11") ~ "Fall"),
         season = relevel(as.factor(season), ref = "Winter"))

labs <- readRDS("/mnt/project/biomarkers_3/covariate_res/res_labs.rds")

fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

colnames(labs) <- c("eid", fields %>% filter(field_id %in% as.numeric(colnames(labs))) %>% pull(title))

data <- df %>%
  left_join(covs) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  #left_join(labs) %>%
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
library(tidyverse)
library(broom)

vars <- c("time_day", "age_recruitment", "sex",
          "chrono", "h_sleep", "ever_insomnia",
          "season", "night_shift", "smoking", "bmi", "is_dst", "wakeup")

# Loop over predictors
results <- map_dfr(vars, function(v) {
  f <- as.formula(paste0("res ~ ", "`",v, "`"))
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
      predictor %in% c("chrono", "h_sleep", "ever_insomnia", "wakeup") ~ "Sleep",
      predictor %in% c("season", "is_dst") ~ "Season",
      predictor == "night_shift" ~ "Job",
      predictor == "time_day" ~ "Other",
      predictor %in% colnames(labs) ~ "Labs",
      TRUE ~ "Other"
    )
  )

# Level lookup for ordering
factor_lookup <- map_dfr(vars, function(v) {
  if (is.factor(data[[v]])) {
    tibble(predictor = v,
           levels_list = list(levels(data[[v]])))
  } else {
    tibble(predictor = v,
           levels_list = list(NULL))
  }
})



# Add display labels and reorder
res2 <- res %>%
  mutate(
    level = str_remove(term, paste0("^", predictor)),
    level_r = ifelse(reference, paste0(level, " (ref)"), level),
    display_term = ifelse(level == "", predictor, level_r)
  ) %>%
  left_join(factor_lookup, by = "predictor") %>%
  group_by(predictor) %>%
  mutate(
    display_term = {
      lvls <- levels_list[[1]]
      if (!is.null(lvls)) {
        lvls <- unique(lvls)
        factor(as.character(display_term),
               levels = c(paste0(lvls[1], " (ref)"),
                          paste0(lvls[-1])))
      } else {
        factor(as.character(display_term),
               levels = unique(as.character(display_term)))
      }
    }
  )





# Pretty predictor labels
pretty_predictor <- c(
  time_day = "Time of Day",
  age_recruitment = "Age at Recruitment",
  sex = "Sex",
  chrono = "Chronotype",
  h_sleep = "Sleep Duration",
  ever_insomnia = "Insomnia",
  season = "Season",
  is_dst = "Daylight savings",
  night_shift = "Night Shift",
  smoking = "Smoking",
  wakeup = "Waking easiness",
  bmi = "BMI"
)

res_plot <- res2 %>%
  mutate(
    predictor_label = pretty_predictor[predictor],
    predictor_label = factor(predictor_label,
                             levels = pretty_predictor)
  )

# --- Plot ---


domain_colors <- c(
  "Demographics" = "#d62728",
  "Job" = "#1f77b4",
  "Other" = "#2ca02c",
  "Season" = "#9467bd",
  "Sleep" = "#ff7f0e"
)


d2 <- res_plot %>%
  mutate(
    Domain = factor(Category, levels = c("Other", "Demographics", "Sleep", "Season", "Job"))
  )


p <- ggplot(d2, aes(x = OR, y = fct_rev(display_term), color = Domain)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(aes(shape = reference), size = 3, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21),
                     labels = c(`FALSE` = "Estimate", `TRUE` = "Reference")) +
  scale_color_manual(values = domain_colors) +
  facet_grid(rows = vars(Domain, predictor_label), scales = "free", space = "free") +
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

ggsave("plots/FX_phenotypes.png", p, width = 12, height = 7)



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
