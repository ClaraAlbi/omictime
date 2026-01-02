

phase <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  filter(i == 0) %>%
  left_join(data.table::fread("sleep_phase.csv")) %>%
  left_join(data.table::fread("/mnt/project/chronotype2.tsv") %>%
              select(eid, chrono = `1180-0.0`) %>%
              mutate(chrono = case_when(
                chrono == 1 ~ "Definitely morning",
                chrono == 2 ~ "Rather morning",
                chrono == -1~ "Don't know",
                chrono == 3 ~ "Rather evening",
                chrono == 4 ~ "Definitely evening",
                TRUE ~ NA_character_)))

phase$res <- residuals(lm(pred_mean ~ time_day, data = phase))

phase <- phase %>%
  mutate(
    condition_wrong_time = case_when(
      p30544 == "Yes" ~ TRUE,
      p30544 == "No"  ~ FALSE,
      TRUE ~ NA
    ),

    condition_stay_up_late = case_when(
      p30545 == "Yes" ~ TRUE,
      p30545 == "No"  ~ FALSE,
      is.na(p30545) & !condition_wrong_time ~ FALSE,
      TRUE ~ NA
    ),

    condition_sleep_normal_waking_hours = case_when(
      p30546 == "Yes" ~ TRUE,
      p30546 == "No"  ~ FALSE,
      is.na(p30546) & !condition_wrong_time ~ FALSE,
      TRUE ~ NA
    ),

    condition_eveningness = chrono  %in% c(
      "Definitely evening",
      "Rather evening"
    ),

    delayed_sleep_wake_phase_disorder =
      condition_wrong_time &
      condition_stay_up_late &
      condition_sleep_normal_waking_hours &
      condition_eveningness
  )

table(phase$delayed_sleep_wake_phase_disorder)

tidy(lm(res ~ delayed_sleep_wake_phase_disorder, data = phase %>% filter(chrono %in% c("Definitely evening", "Rather evening"))))


dat <- phase %>%
  mutate(
    condition_shift_schedule = case_when(
      is.na(shift_work_schedule) ~ NA,
      grepl("shifts", shift_work_schedule) ~ TRUE,
      TRUE ~ FALSE
    ),

    condition_sleep_problem = case_when(
      p30432 %in% c(
        "Yes, a minor problem",
        "Yes, a considerable problem",
        "Yes, a serious problem"
      ) ~ TRUE,
      p30432 == "No" ~ FALSE,
      !condition_shift_schedule ~ FALSE,
      TRUE ~ NA
    ),

    condition_doze_off = case_when(
      p30434 %in% c(
        "Slightly likely",
        "Moderately likely",
        "Highly likely"
      ) ~ TRUE,
      p30434 == "Not likely at all" ~ FALSE,
      !condition_shift_schedule ~ FALSE,
      TRUE ~ NA
    ),

    condition_poor_wellbeing = case_when(
      p30433 %in% c("fairly bad", "very bad") ~ TRUE,
      p30433 %in% c("very good", "fairly good") ~ FALSE,
      !condition_shift_schedule ~ FALSE,
      TRUE ~ NA
    ),

    shift_work_disorder =
      condition_shift_schedule &
      (condition_sleep_problem | condition_doze_off) &
      condition_poor_wellbeing
  )



