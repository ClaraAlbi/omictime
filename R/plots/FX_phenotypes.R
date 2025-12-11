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
library(ggh4x)

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  filter(smoking != "-3") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
         assessment_centre = as.factor(assessment_centre)
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


#Given previously established U-shape relationships with health and cognition [20], we categorised sleep duration into short (<7 h), normal (7–9 h) and long (>9 h) based on recent guidelines

l <- c(list.files("/mnt/project/circadian/results/models",
                  pattern = "predictions", full.names = TRUE))

df_p <- tibble(f = l[str_detect(l, "tech_14")]) %>%
  mutate(d = map(f, readRDS)) %>%
  unnest(d) %>%
  rowwise() %>% mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  unnest()
df_p$res <- residuals(lm(pred_mean ~ time_day, data = df_p))

df_temp <- df_p %>% select(-f) %>%
  inner_join(readRDS("/mnt/project/biomarkers/time.rds")) %>%
  filter(!is.na(time_day)) %>%
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

data <- df %>%
  left_join(covs) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  left_join(gen_covs) %>%
  left_join(dep) %>%
  mutate(h = round(time_day, 0)) %>%
  filter(age_recruitment > 39)

data$res <- residuals(lm(pred_mean ~ time_day, data = data))
#data$res_q <- ntile(data$res, 5)


broom::tidy(lm(res ~ f_reasoning, data = data))

### chrono x night
r <- broom::tidy(lm(res ~ chrono*night_shift + sex + age_recruitment + assessment_centre + PC1 +
     PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data= data))

rs <- broom::tidy(lm(res ~ chrono*shift_work, data= data))

data %>%
  count(chrono, night_shift)

a <- data %>%
  mutate(c = case_when(chrono %in% c("Definitely morning", "Rather morning") ~ "Morning",
                       chrono %in% c("Definitely evening", "Rather evening") ~ "Evening",
                       TRUE ~ "Don't know")) %>%
  filter(!is.na(night_shift)) %>%
  filter(night_shift %in% c("Never", "Always", "Sometimes")) %>%
  #filter((chrono %in% c("Definitely morning", "Rather morning") & night_shift %in% c("Never", "Always")) | (chrono %in% c("Definitely evening", "Rather evening") & night_shift %in% c("Never", "Always")) |
  #         (chrono == "Don't know" & night_shift %in% c("Never", "Always")) ) %>% unite("comb", chrono, night_shift, sep = "_") %>%
  unite("comb", c, night_shift, sep = "_") %>%
  mutate(comb = relevel(as.factor(comb), ref = "Morning_Never")) %>%
  count(comb)

r <- broom::tidy(lm(res ~ comb + sex + age_recruitment + assessment_centre + PC1 +
                 PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data= a)) %>%
  filter(str_detect(term, "comb") | str_detect(term, "Intercept"))

### PLOT DESC

ps <- data %>%
  mutate(res_ab = abs(res)) %>%
  pivot_longer(c(res, res_ab)) %>%
  mutate(name = case_when(name == "res" ~ "Circadian acceleration",
                          name == "res_ab" ~ "Circadian misalignment")) %>%
  ggplot(aes(x = age_recruitment, y = value, color = sex)) + geom_smooth() +
  facet_wrap(~name, scales = "free_y") +
  labs(x = "Age", y = "Value", color = "Sex") +
  theme_classic(base_size = 16) +
  theme(legend.position.inside = c(0.95, 0.95))

ggsave("plots/FS_sexage_CA.png", ps, width = 10, height = 5)



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


tab_desc <- table1::table1(~ age_recruitment + sex + TDI + bmi + smoking + chrono + h_sleep + wakeup + ever_insomnia + season + is_dst + shift_work + night_shift ,
                           data = data,
                           render.cont = my_render_cont)


# --- 1. Predictor list ---
vars <- c("time_day", "age_recruitment", "sex", "chrono", "h_sleep", "ever_insomnia",
          "season", "night_shift", "smoking", "bmi", "is_dst", "wakeup", "shift_work", "TDI")

#vars <- colnames(MH)[-1]

covars <- c("sex", "age_recruitment", "assessment_centre", paste0("PC", 1:10))

#data <- data %>% filter(eth == 1001)

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



vars <- c("age_recruitment", "sex", "chrono", "h_sleep", "ever_insomnia",
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



