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

pred <- readRDS("/mnt/project/olink_int_panels14.rds") %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds")) %>%
  #filter(i == 0) %>%
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
                                `826-0.0` == 4 ~ "Always"),
         employment = case_when(`6142-0.0` == 1	~ "Employed",
                                 `6142-0.0` == 2	~ "Retired",
                                 `6142-0.0` == 3 ~ "Home/family manager",
                                 `6142-0.0` == 4 ~	"Disability",
                                 `6142-0.0` == 5 ~	"Unemployed",
                                 `6142-0.0` == 6 ~	"Voluntary work",
                                 `6142-0.0` == 7 ~	"Student")) %>%
  mutate(night_shift = factor(night_shift, levels = c("Never", "Sometimes", "Usually", "Always")),
         shift_work = factor(shift_work, levels = c("Never/rarely", "Sometimes", "Usually", "Always")),
         employment = fct_relevel(employment, "Employed", "Retired"))

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
         wakeup = `1170-0.0`,
         snoring = `1210-0.0`) %>%
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
                              ever_insomnia == 3 ~ "Usually",
                              TRUE ~ NA_character_),
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
    day_of_week = relevel(day_of_week, ref = "Wed"),
    is_weekend = wday(date_bsampling, week_start = 1) == 6,
    day_type = if_else(is_weekend, "Weekend", "Weekday"),
    day_type = factor(day_type, levels = c("Weekday", "Weekend")),
    fri_sun = if_else(day_of_week == "Fri", "Fri", if_else(is_weekend, "Sat", if_else(day_of_week == "Mon", "Mon",NA_character_))),
    fri_sun = factor(fri_sun, levels = c("Sat", "Fri", "Mon")),
    year = y
  ) %>%
  left_join(dst_transitions, by = "year") %>%
  mutate(
    date_bsampling = as.Date(date_bsampling),

    # SPRING DST classification
    springDST = case_when(
      date_bsampling %in% c(spring_dst - 2, spring_dst - 1) ~ "Before spring DST",
      date_bsampling %in% c(spring_dst + 1, spring_dst + 2) ~ "After spring DST",
      between(date_bsampling, spring_dst - 14, spring_dst - 3) ~ "Baseline spring",
      TRUE ~ NA_character_  # Set all other values to NA
    ),

    # AUTUMN DST classification
    autumnDST = case_when(
      date_bsampling %in% c(fall_dst - 2, fall_dst - 1) ~ "Before autumn DST",
      date_bsampling %in% c(fall_dst + 1, fall_dst + 2) ~ "After autumn DST",
      between(date_bsampling, fall_dst - 14, fall_dst - 3) ~ "Baseline autumn",
      TRUE ~ NA_character_
    ),

    # Optional: convert to factors (NA is preserved)
    springDST = factor(
      springDST,
      levels = c("Baseline spring", "Before spring DST", "After spring DST")
    ),
    autumnDST = factor(
      autumnDST,
      levels = c("Baseline autumn", "Before autumn DST", "After autumn DST")
    )
  )

data <- df %>%
  left_join(covs) %>%
  left_join(a) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  left_join(gen_covs) %>%
  left_join(dep) %>%
  #left_join(phy) %>%
  #left_join(mh)
  #left_join(sleep_new %>% select(eid, rmeq_score, rmeq_chronotype)) %>%
  left_join(vars_join %>% select(eid, chrono_Nightshift)) %>%
  left_join(cohort_f) %>%
  mutate(has_prescription = ifelse(eid %in% cohort$eid, 1, 0),
         antihypertensive = ifelse(!is.na(C0), as.integer(C0 == TRUE), 0),
         sleep_medication = ifelse(!is.na(N05C), as.integer(N05C == TRUE), 0),
         antidepressants = ifelse(!is.na(N06), as.integer(N06 == TRUE), 0),
         mood_stabiliser = ifelse(!is.na(N03), as.integer(N03 == TRUE), 0),
         lithium = ifelse(!is.na(N05A), as.integer(N05A == TRUE), 0),
         across(c(has_prescription, antihypertensive, sleep_medication, antidepressants, mood_stabiliser, lithium), as.factor))


data$res <- residuals(lm(pred_mean ~ time_day, data = data))
#data$res_all <- residuals(lm(`14panels` ~ time_day, data = data))
#data$female <- residuals(lm(female ~ time_day, data = data, na.action = "na.exclude"))
#data$male <- residuals(lm(male ~ time_day, data = data, na.action = "na.exclude"))
#data$raw <- residuals(lm(raw ~ time_day, data = data))
#data$tech <- residuals(lm(tech ~ time_day, data = data, na.action = "na.exclude"))

library(table1)
my_render_cont <- function(x){
  with(
    stats.apply.rounding(stats.default(x)),
    c(
      "",
      `Mean (SD)` = sprintf("%s (%s)", MEAN, SD),
      `Median [Q1, Q3]` = sprintf("%s [%s, %s]", MEDIAN, Q1, Q3),
      `Min, Max` = sprintf("%s, %s", MIN, MAX)
    )
  )
}


tab_desc <- table1::table1(~ age_recruitment + sex + p30079 + TDI + bmi + smoking + season + day_type + fri_sun + autumnDST + springDST + chrono + h_sleep + wakeup + ever_insomnia + shift_work + night_shift + employment + chrono_Nightshift +has_prescription + antihypertensive + sleep_medication + antidepressants + mood_stabiliser + lithium,
                           data = data,
                           render.cont = my_render_cont, topclass="Rtable1-grid")


vars <- c("time_day", "age_recruitment", "sex", "chrono", "h_sleep", "ever_insomnia", "p30079", #"rmeq_chronotype", "rmeq_score",
          "is_weekend", "day_of_week",
          #"sbp", "dbp", "rec_rest", "rec_dep", "rec_tired", "rec_enth",
          "season", "night_shift", "smoking", "bmi", "is_dst", "wakeup", "shift_work", "employment", "TDI", "autumnDST", "springDST",
          "day_type", "fri_sun",
          "chrono_Nightshift",
          "has_prescription",
          "antihypertensive",
          "sleep_medication",
          "antidepressants",
          "mood_stabiliser",
          "lithium")

covars <- c("sex", "age_recruitment", "assessment_centre", paste0("PC", 1:20))


results <- map_dfr(vars, function(v) {
  adj_vars <- if (v %in% c("sex", "age_recruitment")) paste0("PC", 1:20) else covars

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

saveRDS(results, "data_share/results_associations_phenotypes_CA_mean.rds")

