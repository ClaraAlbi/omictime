
install.packages("ggh4x")
library(ggh4x)

demo <- c("age_recruitment", "sex", "p30079", "bmi", "smoking", "TDI")
sleep <- c("chrono", "wakeup", "h_sleep", "ever_insomnia")
season <- c("season", "day_type", "fri_sun", "springDST", "autumnDST")
meds <- c("has_prescription","sleep_medication", "antihypertensive", "antidepressants", "mood_stabiliser", "lithium")
shift <- c("employment","shift_work", "night_shift")

results <- readRDS("data_share/results_associations_phenotypes_CA_mean.rds") %>%
  filter(!str_ends(term, "0")) %>%
  mutate(predictor = case_when(str_detect(term, "p30079") ~ "p30079",
                               TRUE ~ predictor),
         term = str_remove(term, "1"),
         in_mins = estimate * 60) %>%
  filter(predictor %in% c(demo, sleep, season, meds, shift)) %>%
  mutate(FDR = p.adjust(p.value))



# 2. SLEEP


pretty_predictor <- c(
  chrono = "Chronotype", h_sleep = "Sleep Duration", ever_insomnia = "Insomnia", wakeup = "Sleep inertia")

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
    predictor_label = pretty_predictor[as.character(predictor)]) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))


p2 <- plot_data_sleep %>%
  ggplot(aes(x = estimate, y = display_term, shape = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(size = 3, fill = "white") +
  scale_shape_manual(
    values = c(`TRUE` = 15, `FALSE` = 22),
    na.value = 20
  ) +
  #scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  facet_nested(
    rows = vars(factor(predictor_label, levels=unique(plot_data_sleep$predictor_label))),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title = "Sleep questionnaire", x = "CA (β, 95% CI)", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 12) +
  guides(alpha = "none") +
  theme(
    panel.grid.major = element_line(color = "gray"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title.position = "plot"
  )
p2



pretty_predictor <- c(
  season = "Season", day_type = "Social jetlag",
  #fri_sun = "Week day",
  autumnDST = "Fall DST", springDST = "Spring DST")

term_order <- list(
  season = c("Winter (ref)", "Spring", "Summer", "Fall"),
  day_type = c("Weekday (ref)", "Weekend"),
  #fri_sun = c("Sat (ref)", "Fri", "Mon"),
  autumnDST = c("Baseline autumn (ref)", "Before autumn DST", "After autumn DST"),
  springDST = c("Baseline spring (ref)", "Before spring DST", "After spring DST"))


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
    predictor_label = pretty_predictor[as.character(predictor)]) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))


p3 <- plot_data_season %>%
  ggplot(aes(x = estimate, y = display_term, shape = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(size = 3, fill = "white") +
  scale_shape_manual(
    values = c(`TRUE` = 15, `FALSE` = 22),
    na.value = 20
  ) +
  facet_nested(
    rows = vars(factor(predictor_label, levels=unique(plot_data_season$predictor_label))),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(title = "Seasonal effects", x = "CA (β, 95% CI)", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 12) +
  guides(alpha = "none") +
  theme(
    panel.grid.major = element_line(color = "gray"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title.position = "plot"
  )
p3



## employment

pretty_predictor <- c(employment = "Employment status", shift_work = "Shift Work", night_shift = "Night Shift")

term_order <- list(
  employment = c("Employed (ref)", "Retired", "Homemaker","Disability","Unemployed","Voluntary work","Student"),
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
  mutate(term = case_when(term == "employmentHome/family manager" ~ "employmentHomemaker",
                          TRUE ~ term)) %>%
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
    predictor_label = pretty_predictor[as.character(predictor)]
  ) %>% arrange(y_order) %>%
  mutate(
    y_order_rev = -y_order,
    y_id = factor(y_id, levels = y_id[order(y_order_rev)])
  )

p5 <- plot_data_shift %>%
  ggplot(aes(x = estimate, y = y_id, shape = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
  geom_point(size = 3, fill = "white") +
  scale_shape_manual(
    values = c(`TRUE` = 15, `FALSE` = 22),
    na.value = 20
  ) +
  facet_nested(
    rows = vars(factor(predictor_label, levels = unique(plot_data_shift$predictor_label))),
    scales = "free_y",
    space = "free_y",
  ) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_discrete(labels = function(x) sub(".*___", "", x)) +
  labs(title = "Employment", x = "CA (β, 95% CI)", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 12) +
  guides(alpha = "none") +
  theme(
    panel.grid.major = element_line(color = "gray"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.title.position = "plot"
  )

p5

