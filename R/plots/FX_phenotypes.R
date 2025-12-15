library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("broom")
install.packages("table1")
install.packages("forcats")
install.packages("ggh4x")
library(broom)
library(forcats)
library(lubridate)
library(ggh4x)

pred <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  inner_join(readRDS("/mnt/project/biomarkers/time.rds")) %>%
  filter(i == 0) %>%
  filter(!is.na(time_day)) %>%
  mutate(date = as.POSIXct(date_bsampling, tz = "Europe/London"))

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
         assessment_centre = as.factor(assessment_centre)
  )

a <- data.table::fread("/mnt/project/ancestry_new.csv") %>%
  mutate(p30079 = case_when(p30079 == "" ~ NA_character_,
                            TRUE ~ p30079),
         p30079 = relevel(as.factor(p30079), ref = "European ancestry (EUR)"))

job_vars <- data.table::fread("/mnt/project/job_vars.tsv") %>%
  mutate(night_shift = case_when(`3426-0.0` == 1 ~ "Never",
                                 `3426-0.0` == 2 ~ "Sometimes",
                                 `3426-0.0` == 3 ~ "Usually",
                                 `3426-0.0` == 4 ~ "Always"),
         shift_work = case_when(`826-0.0` == 1 ~ "Never/rarely",
                                `826-0.0` == 2 ~ "Sometimes",
                                `826-0.0` == 3 ~ "Usually",
                                `826-0.0` == 4 ~ "Always")) %>%
  mutate(night_shift = factor(night_shift, levels = c("Never", "Sometimes", "Usually", "Always")),
         shift_work = factor(shift_work, levels = c("Never/rarely", "Sometimes", "Usually", "Always")))

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20", `22006-0.0`, `22021-0.0`, `22000-0.0`)
colnames(gen_covs) <- c("eid", paste0("PC", 1:20), "is_white", "rel", "batch")

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
      wakeup == 2 ~ "Not very easy",
      wakeup == -1 ~ "Do not know",
      wakeup == 3 ~ "Fairly easy",
      wakeup == 4 ~ "Very easy"),
    wakeup = factor(wakeup, levels = c("Very easy", "Fairly easy",  "Not very easy", "Not at all easy")),
    chrono = factor(chrono, levels = c("Definitely morning", "Rather morning", "Don't know", "Rather evening", "Definitely evening")),
    ever_insomnia = case_when(ever_insomnia == 1 ~ "Never/rarely",
                              ever_insomnia == 2 ~ "Sometimes",
                              ever_insomnia == 3 ~ "Usually"),
    ever_insomnia = factor(ever_insomnia, levels = c("Never/rarely", "Sometimes", "Usually")),
    h_sleep = case_when(h_sleep > 0 & h_sleep < 7 ~ "Short (<7 h)",
                        h_sleep >= 7 & h_sleep <=9 ~ "Normal (7-9h)",
                        h_sleep > 9 ~ "Long (>9h)"),
    h_sleep = factor(h_sleep, levels = c("Normal (7-9h)", "Short (<7 h)", "Long (>9h)")))

vars_join <- sleep %>%
  left_join(job_vars) %>%
  mutate(
    c = case_when(
      chrono %in% c("Definitely morning", "Rather morning") ~ "Morning",
      chrono %in% c("Definitely evening", "Rather evening") ~ "Evening",
      TRUE ~ "Don't know"
    )
  ) %>%
  filter(!is.na(night_shift)) %>%
  filter(night_shift %in% c("Never", "Always", "Sometimes")) %>%
  unite("chrono_Nightshift", c, night_shift, sep = "_") %>%
  mutate(chrono_Nightshift = relevel(factor(chrono_Nightshift), ref = "Morning_Never")) %>% select(-shift_work)


#Given previously established U-shape relationships with health and cognition [20],
#we categorised sleep duration into
# short (<7 h), normal (7–9 h) and long (>9 h) based on recent guidelines


# 1. Identify DST transitions for each year
years <- unique(year(pred$date))

dst_transitions <- tibble(
  year = c(2006, 2007, 2008, 2009),
  spring_dst = as.Date(c("2006-03-26", "2007-03-25", "2008-03-30", "2009-03-29")),
  fall_dst   = as.Date(c("2006-10-29", "2007-10-28", "2008-10-26", "2009-10-25"))
)

# 2. Clean and engineer features from df_temp
df <- pred %>%
  mutate(
    date = as.POSIXct(date_bsampling, tz = "Europe/London"),
    is_dst = factor(dst(date), levels = c(FALSE, TRUE), labels = c("No", "Yes")),
    date_only = as.Date(date),
    y = year(date_bsampling),
    m = sprintf("%02d", month(date_bsampling)),
    d = sprintf("%02d", day(date_bsampling)),
    season = case_when(
      m %in% c("12", "01", "02") ~ "Winter",
      m %in% c("03", "04", "05") ~ "Spring",
      m %in% c("06", "07", "08") ~ "Summer",
      TRUE ~ "Fall"
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Fall")),
    day_of_week = factor(wday(date_bsampling, label = TRUE, week_start = 1), ordered = F),
    day_of_week = relevel(day_of_week, ref = "Mon"),
    is_weekend = wday(date_bsampling, week_start = 1) == 6,
    day_type = if_else(is_weekend, "Weekend", "Weekday"),
    year = y
  ) %>%
  left_join(dst_transitions, by = "year") %>%
  mutate(
    date_bsampling = as.Date(date_bsampling),

    # SPRING DST classification
    springDST = case_when(
      date_bsampling %in% c(spring_dst - 2, spring_dst - 1) ~ "before_spring_DST",
      date_bsampling %in% c(spring_dst + 1, spring_dst + 2) ~ "after_spring_DST",
      between(date_bsampling, spring_dst - 14, spring_dst - 3) ~ "baseline_spring",
      TRUE ~ NA_character_  # Set all other values to NA
    ),

    # AUTUMN DST classification
    autumnDST = case_when(
      date_bsampling %in% c(fall_dst - 2, fall_dst - 1) ~ "before_fall_DST",
      date_bsampling %in% c(fall_dst + 1, fall_dst + 2) ~ "after_fall_DST",
      between(date_bsampling, fall_dst - 14, fall_dst - 3) ~ "baseline_fall",
      TRUE ~ NA_character_
    ),

    # Optional: convert to factors (NA is preserved)
    springDST = factor(
      springDST,
      levels = c("baseline_spring", "before_spring_DST", "after_spring_DST")
    ),
    autumnDST = factor(
      autumnDST,
      levels = c("baseline_fall", "before_fall_DST", "after_fall_DST")
    )
  )

data <- df %>%
  left_join(covs) %>%
  left_join(a) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  left_join(gen_covs) %>%
  left_join(dep) %>%
  left_join(vars_join) %>%
  filter(age_recruitment > 40 & age_recruitment < 70)

data$res <- residuals(lm(pred_mean ~ time_day + as.factor(assessment_centre), data = data))


### trials

rs <- broom::tidy(lm(res ~ p30079 + sex + age_recruitment + assessment_centre + PC1 +
                       PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + season, data= data))

### chrono x night
r <- broom::tidy(lm(res ~ chrono*night_shift + sex + age_recruitment + assessment_centre + PC1 +
     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + year, data= data))

rs <- broom::tidy(lm(res ~ chrono*shift_work, data= data))

data %>%
  count(chrono, night_shift)

### PLOT DESC

# ps <- data %>%
#   filter(age_recruitment > 40 & age_recruitment < 70) %>%
#   ggplot(aes(x = age_recruitment, y = res, color = sex)) + geom_smooth() +
#   labs(x = "Age", y = "Circadian acceleration", color = "Sex") +
#   theme_classic(base_size = 16) +
#   theme(legend.position.inside = c(0.95, 0.95))
#
# ggsave("plots/FS_sexage_CA.png", ps, width = 8, height = 5)


### Plots CHRONO vs SEX

p_chrono <- sleep %>%
  left_join(covs) %>%
  filter(age_recruitment > 39 & age_recruitment <= 70) %>%
  filter(!is.na(chrono)) %>%
  group_by(age_recruitment,sex,  chrono) %>%
  count() %>% ungroup() %>%
  group_by(age_recruitment, sex) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = age_recruitment, y = p, color = chrono)) +
  geom_point() +
  labs(color = "Chronotype", y = "Proportion", x = "Age") +
  facet_grid(~sex) +
  theme_classic(base_size = 16)

ggsave("plots/FS_chronotypes_agesex.png", p_chrono, width = 10, height = 5)


ca_res <- data %>%
  filter(!is.na(chrono)) %>%
  ggplot(aes(x = age_recruitment, y = res, color = sex)) + geom_smooth() +
  labs(x = "Age", y = "Circadian acceleration", color = "Sex") +
  facet_grid(~chrono) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 8),
        legend.position.inside = c(0.95, 0.95))

ggsave("plots/FS_CA_chronotypes_agesex.png", ca_res, width = 10, height = 5)

## associations

m_0 <- tidy(lm(res ~ age_recruitment, data = data %>% filter(sex == "Female")))
# term            estimate std.error statistic p.value
# <chr>              <dbl>     <dbl>     <dbl>   <dbl>
#   1 (Intercept)     -0.0816   0.0504       -1.62   0.106
# 2 age_recruitment  0.00142  0.000881      1.61   0.107

m_1 <- tidy(lm(res ~ age_recruitment, data = data %>% filter(sex == "Male")))
# term            estimate std.error statistic p.value
# <chr>              <dbl>     <dbl>     <dbl>   <dbl>
#   1 (Intercept)      0.0704   0.0556        1.27   0.205
# 2 age_recruitment -0.00121  0.000964     -1.26   0.209

m_c <- tidy(lm(res ~ age_recruitment + sex*chrono, data = data))
m_1_c <- tidy(lm(res ~ age_recruitment*chrono, data = data %>% filter(sex == "Male")))

library(ggplot2)

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


tab_desc <- table1::table1(~ age_recruitment + sex + p30079 + TDI + bmi + smoking + season + is_weekend + autumnDST + springDST + chrono + h_sleep + wakeup + ever_insomnia + shift_work + night_shift + chrono_Nightshift,
                           data = data,
                           render.cont = my_render_cont, topclass="Rtable1-grid")


# --- 1. Predictor list ---
vars <- c("time_day", "age_recruitment", "sex", "chrono", "h_sleep", "ever_insomnia", "p30079",
          "season", "night_shift", "smoking", "bmi", "is_dst", "wakeup", "shift_work", "TDI", "autumnDST", "springDST", "chrono_Nightshift")

covars <- c("sex", "age_recruitment", "assessment_centre", paste0("PC", 1:20))

# --- 2. Fit models and extract results ---
results <- map_dfr(vars, function(v) {
  adj_vars <- if (v %in% c("sex", "age_recruitment")) paste0("PC", 1:10) else covars

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



saveRDS(results, "data_share/results_associations_phenotypes_CA.rds")


### PLOT PART



vars <- c("age_recruitment", "sex", "chrono", "h_sleep", "ever_insomnia", "p30079",
          "season", "night_shift", "smoking", "bmi", "is_dst", "wakeup", "shift_work", "TDI")

covars <- c("sex", "age_recruitment", paste0("PC", 1:10))

results <- bind_rows(readRDS("data_share/results_associations_phenotypes_CA.rds"),
                     readRDS("data_share/results_associations_medication_CA.rds") %>% filter(str_ends(term, "1")) %>%
                       mutate(term = str_remove(term, "1")))


# Define pretty labels and colors
pretty_predictor <- c(
  age_recruitment = "Age at Recruitment", sex = "Sex",
  chrono = "Chronotype", h_sleep = "Sleep Duration", ever_insomnia = "Insomnia",
  season = "Season", is_dst = "Daylight Savings", night_shift = "Night Shift",
  smoking = "Smoking", wakeup = "Waking Easiness", bmi = "BMI",
  shift_work = "Shift Work", TDI = "Townsend DI",
  has_prescription = "Any",
  antidepressants = "Antidepressants", antihypertensive = "Antihypertensives", sleep_medication = "Sedatives and hypnotics", mood_stabiliser = "Mood stabilisers", lithium = "Lithium")

domain_colors <- c(
  "Demographics" = "#d62728", "Job" = "#1f77b4",
  "Season" = "#9467bd", "Sleep \nquestionnaire" = "#ff7f0e", "Sleep medication" = "darkgreen")

# Define predictor order (standard names)
predictor_order <- c("sex", "age_recruitment", "bmi", "smoking", "TDI",
                     "season", "is_dst", "chrono", "wakeup", "h_sleep", "ever_insomnia",
                     "shift_work", "night_shift",
                     "has_prescription","sleep_medication", "antihypertensive", "antidepressants", "mood_stabiliser", "lithium")

# Define term order for each predictor
term_order <- list(
  sex = c("Female (ref)", "Male"),
  age_recruitment = c("Age at Recruitment"),
  bmi = c("BMI"),
  smoking = c("Never (ref)", "Previous", "Current"),
  TDI = c("Townsend DI"),
  season = c("Winter (ref)", "Spring", "Summer", "Fall"),
  is_dst = c("No (ref)", "Yes"),
  chrono = c("Definitely morning (ref)", "Rather morning", "Don't know",
             "Rather evening", "Definitely evening"),
  wakeup = c("Very easy (ref)", "Fairly easy", "Not very easy", "Not at all easy"),
  h_sleep = c("Normal (7-9h) (ref)", "Short (<7 h)", "Long (>9h)"),
  ever_insomnia = c("Never/rarely (ref)", "Sometimes", "Usually"),
  shift_work = c("Never/rarely (ref)", "Usually", "Always"),
  night_shift = c("Never. (ref)", "Sometimes.", "Usually.", "Always."),
  has_prescription = c("Any"),
  antihypertensive = c("Antihypertensives"),
  antidepressants= c("Antidepressants"),
  sleep_medication = c("Sedatives and hypnotics"),
  mood_stabiliser= c("Mood stabilisers"),
  lithium= c("Lithium")
)

# Create ordered dataframe
order_df <- imap_dfr(term_order, ~ tibble(
  predictor = .y,
  display_term = .x,
  term_rank = seq_along(.x))) %>%
  mutate(predictor = factor(predictor, levels = predictor_order)) %>%
  arrange(predictor, term_rank) %>%
  mutate(y_order = row_number())

# Process results - filter out time_day
res <- results %>%
  filter(predictor != "time_day") %>%
  mutate(term = case_when(
    str_detect(term, "night_shift") ~ paste0(term, "."),
    TRUE ~ term
  )) %>%
  mutate(
    lower = ifelse(reference, 0, estimate - 1.96 * std.error),
    upper = ifelse(reference, 0, estimate + 1.96 * std.error),
    Category = case_when(
      predictor %in% c("age_recruitment", "sex", "bmi", "smoking", "TDI") ~ "Demographics",
      predictor %in% c("chrono", "h_sleep", "ever_insomnia", "wakeup") ~ "Sleep \nquestionnaire",
      predictor %in% c("season", "is_dst") ~ "Season",
      predictor %in% c("shift_work", "night_shift") ~ "Job",
      TRUE ~ "Sleep medication"
    ),
    predictor_label = pretty_predictor[predictor],
    level = str_remove(term, paste0("^", predictor)),
    display_term = ifelse(reference,
                          paste0(level, " (ref)"),
                          ifelse(level %in% c(""), predictor_label, level))
  ) %>%
  mutate(predictor = factor(predictor, levels = predictor_order))

# Join with order_df and create plot data
plot_data <- res %>%
  right_join(order_df, by = c("predictor", "display_term")) %>%
  mutate(
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value),
    predictor_label = pretty_predictor[as.character(predictor)],
    Category = factor(Category, levels = c("Demographics", "Season", "Sleep \nquestionnaire", "Job", "Sleep medication"))
  ) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))

p1 <- plot_data %>%
  ggplot(aes(x = estimate, y = display_term, color = Category, alpha = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(aes(shape = reference), size = 3, fill = "white") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_color_manual(values = domain_colors) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  #scale_y_discrete(labels = term_labels) +
  #scale_x_continuous(limits = c(0.5, 1.2)) +
  facet_nested(
    rows = vars(Category, factor(predictor_label, levels=unique(plot_data$predictor_label))),
    scales = "free_y",
    space = "free_y"
  ) +
  labs(title = "Circadian acceleration", x = "Effect size (SE)", y = NULL, color = " ", alpha = "FDR < 5%") +
  theme_classic(base_size = 14) +
  theme(
    # Keep outer strip styling as fallback if needed
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
    strip.background = element_rect(fill = "antiquewhite2", color = "black", linewidth = 0.8),
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
p1


ggsave("plots/FX_phenotypes_CA.png", p1, width = 10, height = 10)



