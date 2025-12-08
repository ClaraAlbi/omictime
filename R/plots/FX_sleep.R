library(ggplot2)

sleep_new <- data.table::fread("/mnt/project/sleep_online.csv") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  mutate(
    score_wake_up = case_when(
      p30425 == "5:00am - 6:30am" ~ 5,
      p30425 == "6:30am - 7:45am" ~ 4,
      p30425 == "7:45am - 9:45am" ~ 3,
      p30425 == "9:45am - 11:00am" ~ 2,
      p30425 == "11:00am - 12 noon" ~ 1,
      TRUE ~ NA_real_
    ),

    score_tired_after_waking = case_when(
      p30426 == "Very tired" ~ 1,
      p30426 == "Fairly tired" ~ 2,
      p30426 == "Fairly refreshed" ~ 3,
      p30426 == "Very refreshed" ~ 4,
      TRUE ~ NA_real_
    ),

    score_evening_tiredness = case_when(
      p30427 == "8:00pm - 9:00pm" ~ 5,
      p30427 == "9:00pm - 10:15pm" ~ 4,
      p30427 == "10:15pm - 12:45am" ~ 3,
      p30427 == "12:45am - 2:00am" ~ 2,
      p30427 == "2:00am - 3:00am" ~ 1,
      TRUE ~ NA_real_
    ),

    score_best_time = case_when(
      p30428 == "5:00am - 8:00am" ~ 5,
      p30428 == "8:00am - 10:00am" ~ 4,
      p30428 == "10:00am - 5:00pm" ~ 3,
      p30428 == "5:00pm - 10:00pm" ~ 2,
      p30428 == "10:00pm - 5:00am" ~ 1,
      TRUE ~ NA_real_
    ),

    score_chronotype = case_when(
      p30429 == "Definitely a morning-type" ~ 6,
      p30429 == "Rather more a morning-type than an evening-type" ~ 4,
      p30429 == "Rather more an evening-type than a morning-type" ~ 2,
      p30429 == "Definitely an evening-type" ~ 0,
      TRUE ~ NA_real_
    ),
    rmeq_score = if_else(
      !is.na(score_wake_up) &
        !is.na(score_tired_after_waking) &
        !is.na(score_evening_tiredness) &
        !is.na(score_best_time) &
        !is.na(score_chronotype),

      score_wake_up +
        score_tired_after_waking +
        score_evening_tiredness +
        score_best_time +
        score_chronotype,

      NA_real_
    ),
    rmeq_chronotype = case_when(
      rmeq_score >= 22 ~ "Definitely morning",
      rmeq_score >= 18 & rmeq_score <= 21 ~ "Rather morning",
      rmeq_score >= 12 & rmeq_score <= 17 ~ "Neither",
      rmeq_score >= 8  & rmeq_score <= 11 ~ "Rather evening",
      rmeq_score <= 7  ~ "Definitely evening",
      TRUE ~ NA_character_
    ),
    rmeq_chronotype = factor(rmeq_chronotype, levels = c("Definitely morning", "Rather morning", "Neither", "Rather evening", "Definitely evening")),
    rmeq_chronotype = relevel(rmeq_chronotype, ref = "Definitely morning"),

  )

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  ungroup() %>%
  left_join(sleep_new)
df$res <- residuals(lm(pred_mean ~ time_day, data = df))


df %>%
  ggplot(aes(x = res, y = rmeq_chronotype)) +
  geom_boxplot()

summary(lm(res ~ rmeq_chronotype, data = df))
summary(lm(res ~ p30429, data = df))


#### CHECK IF PEOPLE THAT CHANGE GROUP HAD HIGHER DYSREGULATION

