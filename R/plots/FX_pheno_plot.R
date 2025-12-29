library(tidyverse)
library(ggh4x)
library(patchwork)

### PLOT PART

results <- readRDS("data_share/results_associations_phenotypes_CA.rds") %>%
  filter(!str_ends(term, "0")) %>%
  mutate(predictor = case_when(str_detect(term, "p30079") ~ "p30079",
                               TRUE ~ predictor),
         term = str_remove(term, "1"),
         in_mins = estimate * 60)


# 1. demographics

demo <- c("age_recruitment", "sex", "p30079", "bmi", "smoking", "TDI")

pretty_predictor <- c(
  age_recruitment = "Age at Recruitment", sex = "Sex",
  p30079 = "Genetic ancestry",
  bmi = "BMI", smoking = "Smoking",
  TDI = "Townsend DI")

term_order <- list(
  sex = c("Female (ref)", "Male"),
  age_recruitment = c("Age at Recruitment"),
  p30079 = c("European ancestry (EUR) (ref)","African ancestry (AFR)", "Central/South Asian ancestry (CSA)", "Middle Eastern ancestry (MID)", "East Asian ancestry (EAS)", "Admixed American ancestry (AMR)"),
  bmi = c("BMI"),
  smoking = c("Never (ref)", "Previous", "Current"),
  TDI = c("Townsend DI"))

# Create ordered dataframe
order_df <- imap_dfr(term_order, ~ tibble(
  predictor = .y,
  display_term = .x,
  term_rank = seq_along(.x))) %>%
  mutate(predictor = factor(predictor, levels = demo)) %>%
  arrange(predictor, term_rank) %>%
  mutate(y_order = row_number())

# Process results - filter out time_day
res <- results %>%
  filter(predictor %in% demo) %>%
   mutate(
     lower = ifelse(reference, 0, estimate - 1.96 * std.error),
     upper = ifelse(reference, 0, estimate + 1.96 * std.error),
    predictor_label = pretty_predictor[predictor],
    level = str_remove(term, paste0("^", predictor)),
    display_term = ifelse(reference,
                          paste0(level, " (ref)"),
                          ifelse(level %in% c(""), predictor_label, level))
  ) %>%
  mutate(predictor = factor(predictor, levels = demo))

# Join with order_df and create plot data
plot_data_demo <- res %>%
  right_join(order_df, by = c("predictor", "display_term")) %>%
  mutate(
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value),
    predictor_label = pretty_predictor[as.character(predictor)]) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))


p1 <- plot_data_demo %>%
  ggplot(aes(x = estimate, y = display_term, alpha = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1, color = "#d62728") +
  geom_point(aes(shape = reference), size = 2, fill = "white", color = "#d62728") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  facet_nested(
    rows = vars(factor(predictor_label, levels=unique(plot_data_demo$predictor_label))),
    scales = "free_y",
    space = "free_y"
  ) +
  labs(title = "Socio-demographics", x = "Effect size", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 10) +
  theme(
    # Keep outer strip styling as fallback if needed
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold", size = 8),
    strip.background = element_rect(fill = alpha("#A9CBB7", 0.4), color = "black", linewidth = 0.8),
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
p1

# 2. SLEEP

sleep <- c("chrono", "wakeup", "h_sleep", "ever_insomnia")

pretty_predictor <- c(
  chrono = "Chronotype", h_sleep = "Sleep Duration", ever_insomnia = "Insomnia", wakeup = "Easiness to wake up")

term_order <- list(
  chrono = c("Definitely morning (ref)", "Rather morning", "Don't know",
             "Rather evening", "Definitely evening"),
  wakeup = c("Very easy (ref)", "Fairly easy", "Not very easy", "Not at all easy"),
  h_sleep = c("Normal (7-9h) (ref)", "Short (<7 h)", "Long (>9h)"),
  ever_insomnia = c("Never/rarely (ref)", "Sometimes", "Usually"))

# Create ordered dataframe
order_df <- imap_dfr(term_order, ~ tibble(
  predictor = .y,
  display_term = .x,
  term_rank = seq_along(.x))) %>%
  mutate(predictor = factor(predictor, levels = sleep)) %>%
  arrange(predictor, term_rank) %>%
  mutate(y_order = row_number())

# Process results - filter out time_day
res <- results %>%
  filter(predictor %in% sleep) %>%
  mutate(
    lower = ifelse(reference, 0, estimate - 1.96 * std.error),
    upper = ifelse(reference, 0, estimate + 1.96 * std.error),
    predictor_label = pretty_predictor[predictor],
    level = str_remove(term, paste0("^", predictor)),
    display_term = ifelse(reference,
                          paste0(level, " (ref)"),
                          ifelse(level %in% c(""), predictor_label, level))
  ) %>%
  mutate(predictor = factor(predictor, levels = sleep))

# Join with order_df and create plot data
plot_data_sleep <- res %>%
  right_join(order_df, by = c("predictor", "display_term")) %>%
  mutate(
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value),
    predictor_label = pretty_predictor[as.character(predictor)]) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))


p2 <- plot_data_sleep %>%
  ggplot(aes(x = estimate, y = display_term, alpha = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1, color = "#ff7f0e") +
  geom_point(aes(shape = reference), size = 2, fill = "white", color = "#ff7f0e") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  facet_nested(
    rows = vars(factor(predictor_label, levels=unique(plot_data_sleep$predictor_label))),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title = "Sleep questionnaire", x = "Effect size", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 10) +
  theme(
    # Keep outer strip styling as fallback if needed
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold", size = 8),
    strip.background = element_rect(fill = alpha("#A9CBB7", 0.4), color = "black", linewidth = 0.8),
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
p2

### Seasonality

season <- c("season", "is_weekend","springDST", "autumnDST")

pretty_predictor <- c(
  season = "Season", is_weekend = "Social jetlag",
  autumnDST = "Fall DST", springDST = "Spring DST")

term_order <- list(
  season = c("Winter (ref)", "Spring", "Summer", "Fall"),
  is_weekend = c("TRUE"),
  autumnDST = c("baseline_fall (ref)", "before_fall_DST", "after_fall_DST"),
  springDST = c("baseline_spring (ref)", "before_spring_DST", "after_spring_DST"))



# Create ordered dataframe
order_df <- imap_dfr(term_order, ~ tibble(
  predictor = .y,
  display_term = .x,
  term_rank = seq_along(.x))) %>%
  mutate(predictor = factor(predictor, levels = season)) %>%
  arrange(predictor, term_rank) %>%
  mutate(y_order = row_number())

# Process results - filter out time_day
res <- results %>%
  filter(predictor %in% season) %>%
  mutate(
    lower = ifelse(reference, 0, estimate - 1.96 * std.error),
    upper = ifelse(reference, 0, estimate + 1.96 * std.error),
    predictor_label = pretty_predictor[predictor],
    level = str_remove(term, paste0("^", predictor)),
    display_term = ifelse(reference,
                          paste0(level, " (ref)"),
                          ifelse(level %in% c(""), predictor_label, level))
  ) %>%
  mutate(predictor = factor(predictor, levels = season))

# Join with order_df and create plot data
plot_data_season <- res %>%
  right_join(order_df, by = c("predictor", "display_term")) %>%
  mutate(
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value),
    predictor_label = pretty_predictor[as.character(predictor)]) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))


p3 <- plot_data_season %>%
  ggplot(aes(x = estimate, y = display_term, alpha = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1, color = "darkgreen") +
  geom_point(aes(shape = reference), size = 2, fill = "white", color = "darkgreen") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  facet_nested(
    rows = vars(factor(predictor_label, levels=unique(plot_data_season$predictor_label))),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title = "Seasonal effects", x = "Effect size", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 10) +
  theme(
    # Keep outer strip styling as fallback if needed
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold", size = 8),
    strip.background = element_rect(fill = alpha("#A9CBB7", 0.4), color = "black", linewidth = 0.8),
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
p3


### MEDICATION

meds <- c("has_prescription","sleep_medication", "antihypertensive", "antidepressants", "mood_stabiliser", "lithium")

pretty_predictor <- c(
  has_prescription = "Any",
  antidepressants = "Antidepressants", antihypertensive = "Antihypertensives", sleep_medication = "Sedatives and hypnotics", mood_stabiliser = "Mood stabilisers", lithium = "Lithium")

term_order <- list(
  has_prescription = c("Any"),
  antihypertensive = c("Antihypertensives"),
  antidepressants= c("Antidepressants"),
  sleep_medication = c("Sedatives and hypnotics"),
  mood_stabiliser= c("Mood stabilisers"),
  lithium= c("Lithium"))


# Create ordered dataframe
order_df <- imap_dfr(term_order, ~ tibble(
  predictor = .y,
  display_term = .x,
  term_rank = seq_along(.x))) %>%
  mutate(predictor = factor(predictor, levels = meds)) %>%
  arrange(predictor, term_rank) %>%
  mutate(y_order = row_number())

# Process results - filter out time_day
res <- results %>%
  filter(predictor %in% meds) %>%
  mutate(
    lower = ifelse(reference, 0, estimate - 1.96 * std.error),
    upper = ifelse(reference, 0, estimate + 1.96 * std.error),
    predictor_label = pretty_predictor[predictor],
    level = str_remove(term, paste0("^", predictor)),
    display_term = ifelse(reference,
                          paste0(level, " (ref)"),
                          ifelse(level %in% c(""), predictor_label, level))
  ) %>%
  mutate(predictor = factor(predictor, levels = meds), p = "Medication type")

# Join with order_df and create plot data
plot_data_med <- res %>%
  right_join(order_df, by = c("predictor", "display_term")) %>%
  mutate(
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value),
    predictor_label = pretty_predictor[as.character(predictor)]) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))


p4 <- plot_data_med %>%
  ggplot(aes(x = estimate, y = display_term, alpha = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1, color = "#1f77b4") +
  geom_point(aes(shape = reference), size = 2, fill = "white", color = "#1f77b4") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  scale_x_continuous(limits = c(-1, 1)) +
  facet_nested(
    rows = vars(p),
    scales = "free_y",
    space = "free_y"
  ) +
  labs(title = "Sleep presciptions", x = "Effect size", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 10) +
  theme(
    # Keep outer strip styling as fallback if needed
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold", size = 8),
    strip.background = element_rect(fill = alpha("#A9CBB7", 0.4), color = "black", linewidth = 0.8),
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
p4


#### JOB

shift <- c("shift_work", "night_shift")

pretty_predictor <- c(shift_work = "Shift Work", night_shift = "Night Shift")

term_order <- list(
  shift_work = c("Never/rarely (ref)", "Usually", "Always"),
  night_shift = c("Never (ref)", "Sometimes", "Usually", "Always"))


# Create ordered dataframe
order_df <- imap_dfr(term_order, ~ tibble(
  predictor = .y,
  display_term = .x,
  term_rank = seq_along(.x))) %>%
  mutate(predictor = factor(predictor, levels = shift)) %>%
  arrange(predictor, term_rank) %>%
  mutate(y_order = row_number())

# Process results - filter out time_day
res <- results %>%
  filter(predictor %in% shift) %>%
  mutate(
    lower = ifelse(reference, 0, estimate - 1.96 * std.error),
    upper = ifelse(reference, 0, estimate + 1.96 * std.error),
    predictor_label = pretty_predictor[predictor],
    level = str_remove(term, paste0("^", predictor)),
    display_term = ifelse(reference,
                          paste0(level, " (ref)"),
                          ifelse(level %in% c(""), predictor_label, level))
  ) %>%
  mutate(predictor = factor(predictor, levels = shift))

# Join with order_df and create plot data
plot_data_shift <- res %>%
  right_join(order_df, by = c("predictor", "display_term")) %>%
  # ensure rows are ordered by your y_order
  arrange(y_order) %>%
  # create a unique internal id so factor levels are unique across predictors
  mutate(
    y_id = paste0(predictor, "___", display_term),   # unique key (predictor + label)
    # make y_id a factor with levels in the arranged order (preserves y_order)
    y_id = factor(y_id, levels = unique(y_id)),
    # compute any plotting flags
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value),
    predictor_label = pretty_predictor[as.character(predictor)]
  ) %>% arrange(y_order) %>%
  mutate(
    y_order_rev = -y_order,
    y_id = factor(y_id, levels = y_id[order(y_order_rev)])
  )

p5 <- ggplot(plot_data_shift, aes(x = estimate, y = y_id, alpha = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1, color = "#9467bd") +
  geom_point(aes(shape = reference), size = 2, fill = "white", color = "#9467bd") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  facet_nested(
    rows = vars(factor(predictor_label, levels = unique(plot_data_shift$predictor_label))),
    scales = "free_y",
    space = "free_y",
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title = "Employment - shift work", x = "Effect size", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 10) +
  theme(
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold", size = 8),
    strip.background = element_rect(fill = alpha("#A9CBB7", 0.4), color = "black", linewidth = 0.8),
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  # replace y-axis labels with just the display_term part (after "___")
  scale_y_discrete(labels = function(x) sub(".*___", "", x))
p5



### COMBINE

all_ps <-
(p1 | p2) /
  (p3 | (p5 / p4)) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    plot.margin = margin(5, 5, 5, 5)
  )

ggsave("plots/FX_phenotypes.png", all_ps, width = 10, height = 8)
