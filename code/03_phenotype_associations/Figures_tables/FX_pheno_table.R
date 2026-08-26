


library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(forcats)
library(scales)

# groups
demo   <- c("age_recruitment", "sex", "p30079", "bmi", "smoking", "TDI")
sleep  <- c("chrono", "wakeup", "h_sleep", "ever_insomnia")
season <- c("season", "day_type", "fri_sun", "springDST", "autumnDST")
meds   <- c("has_prescription","sleep_medication", "antihypertensive", "antidepressants",
            "mood_stabiliser", "lithium")
shift  <- c("employment","shift_work", "night_shift")

# group name lookup (predictor -> group)
predictor_group <- c(
  setNames(rep("Demographics", length(demo)), demo),
  setNames(rep("Sleep questionnaire", length(sleep)), sleep),
  setNames(rep("Seasonality", length(season)), season),
  setNames(rep("Medication", length(meds)), meds),
  setNames(rep("Employment", length(shift)), shift)
)

# pretty predictor labels (one named vector for everything)
pretty_predictor <- c(
  age_recruitment = "Age at Recruitment",
  sex = "Sex",
  p30079 = "Genetic ancestry",
  bmi = "BMI",
  smoking = "Smoking",
  TDI = "Townsend DI",

  chrono = "Chronotype",
  wakeup = "Easiness to wake up",
  h_sleep = "Sleep Duration",
  ever_insomnia = "Insomnia",

  season = "Season",
  day_type = "Social jetlag",
  fri_sun = "Fri/Sun indicator",   # keep if you end up using it
  autumnDST = "Fall DST",
  springDST = "Spring DST",

  has_prescription = "Any prescription",
  sleep_medication = "Sedatives and hypnotics",
  antihypertensive = "Antihypertensives",
  antidepressants = "Antidepressants",
  mood_stabiliser = "Mood stabilisers",
  lithium = "Lithium",

  employment = "Employment status",
  shift_work = "Shift Work",
  night_shift = "Night Shift"
)

# your desired display order per predictor (THIS drives the final table order)
term_order <- list(
  # demographics
  sex = c("Female (ref)", "Male"),
  age_recruitment = c("Age at Recruitment"),
  p30079 = c("European ancestry (EUR) (ref)","African ancestry (AFR)", "Central/South Asian ancestry (CSA)",
             "Middle Eastern ancestry (MID)", "East Asian ancestry (EAS)", "Admixed American ancestry (AMR)"),
  bmi = c("BMI"),
  smoking = c("Never (ref)", "Previous", "Current"),
  TDI = c("Townsend DI"),

  # sleep
  chrono = c("Definitely morning (ref)", "Rather morning", "Don't know",
             "Rather evening", "Definitely evening"),
  wakeup = c("Very easy (ref)", "Fairly easy", "Not very easy", "Not at all easy"),
  h_sleep = c("Normal (7-9h) (ref)", "Short (<7 h)", "Long (>9h)"),
  ever_insomnia = c("Never/rarely (ref)", "Sometimes", "Usually"),

  # seasonality
  season = c("Winter (ref)", "Spring", "Summer", "Fall"),
  day_type = c("Weekday (ref)", "Weekend"),
  autumnDST = c("Baseline autumn (ref)", "Before autumn DST", "After autumn DST"),
  springDST = c("Baseline spring (ref)", "Before spring DST", "After spring DST"),

  # meds (all continuous/binary-ish shown as single line in your plot)
  has_prescription = c("Any prescription"),
  antihypertensive = c("Antihypertensives"),
  antidepressants = c("Antidepressants"),
  sleep_medication = c("Sedatives and hypnotics"),
  mood_stabiliser = c("Mood stabilisers"),
  lithium = c("Lithium"),

  # employment
  employment = c("Employed (ref)", "Retired", "Home/family manager","Disability","Unemployed","Voluntary work","Student"),
  shift_work = c("Never/rarely (ref)", "Sometimes", "Usually", "Always"),
  night_shift = c("Never (ref)", "Sometimes", "Usually", "Always")
)

# master ordering dataframe (predictor + display_term -> y_order)
order_df <- imap_dfr(term_order, ~ tibble(
  predictor = .y,
  display_term = .x,
  term_rank = seq_along(.x)
)) %>%
  mutate(
    group = unname(predictor_group[predictor]),
    group = factor(group, levels = c("Demographics","Seasonality","Sleep questionnaire","Employment","Medication")),
    predictor = factor(predictor, levels = c(demo, season, sleep, shift, meds))
  ) %>%
  arrange(group, predictor, term_rank) %>%
  mutate(y_order = row_number())


format_assoc_table <- function(results, order_df, pretty_predictor, predictor_group) {

  all_predictors <- unique(order_df$predictor) |> as.character()

  res <- results %>%
    filter(predictor %in% all_predictors) %>%
    mutate(
      predictor = case_when(str_detect(term, "p30079") ~ "p30079", TRUE ~ predictor),
      term = str_remove(term, "1"),
      term = str_remove(term, "0$"),
      in_mins = estimate * 60,
      FDR = p.adjust(p.value),

      reference = coalesce(as.logical(reference), FALSE),

      predictor_label = unname(pretty_predictor[predictor]),
      level_raw = str_remove(term, paste0("^", predictor)),
      display_term = case_when(
        reference ~ paste0(level_raw, " (ref)"),
        level_raw == "" ~ predictor_label,
        TRUE ~ level_raw
      ),

      lower = if_else(reference, 0, estimate - 1.96 * std.error),
      upper = if_else(reference, 0, estimate + 1.96 * std.error)
    )

  order_df %>%
    right_join(res, by = c("predictor", "display_term")) %>%
    arrange(y_order) %>%
    transmute(
      group = group,                 # comes from order_df
      predictor,
      predictor_label,
      display_term,
      reference,
      n,
      n_total,
      beta = estimate,
      se = std.error,
      ci_low = lower,
      ci_high = upper,
      p = p.value,
      FDR,
      in_mins,
      y_order
    )
}
nice_table <- format_assoc_table(results, order_df, pretty_predictor, predictor_group) %>%
  select(Data = group, Predictor = predictor_label, Category = display_term, n, n_total, beta, beta_mins = in_mins, se, ci_low, ci_high, p, p_adj = FDR)


writexl::write_xlsx(nice_table, "tables/phenotype_associations.xlsx")
saveRDS(nice_table, "tables/phenotype_associations.rds")

