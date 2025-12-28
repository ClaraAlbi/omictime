library(tidyverse)
library(ggh4x)

### PLOT PART

results <- bind_rows(readRDS("data_share/results_associations_phenotypes_CA.rds"),
                     readRDS("data_share/results_associations_medication_CA.rds") %>% filter(str_ends(term, "1")) %>%
                       mutate(term = str_remove(term, "1"))) %>%
  mutate(predictor = case_when(str_detect(term, "p30079") ~ "p30079",
                               TRUE ~ predictor))


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
plot_data <- res %>%
  right_join(order_df, by = c("predictor", "display_term")) %>%
  mutate(
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value),
    predictor_label = pretty_predictor[as.character(predictor)]) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))


p1 <- plot_data %>%
  ggplot(aes(x = estimate, y = display_term, alpha = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1, color = "#d62728") +
  geom_point(aes(shape = reference), size = 2, fill = "white", color = "#d62728") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  facet_nested(
    rows = vars(factor(predictor_label, levels=unique(plot_data$predictor_label))),
    scales = "free_y",
    space = "free_y"
  ) +
  labs(title = "Demographics", x = "Effect size", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 10) +
  theme(
    # Keep outer strip styling as fallback if needed
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold", size = 10),
    strip.background = element_rect(fill = "antiquewhite2", color = "black", linewidth = 0.8),
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
p1

# 2. SLEEP

sleep <- c("chrono", "h_sleep", "ever_insomnia", "wakeup")

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
plot_data <- res %>%
  right_join(order_df, by = c("predictor", "display_term")) %>%
  mutate(
    is_ref = if_else(estimate == 0 & is.na(p.value), TRUE, FALSE),
    FDR = p.adjust(p.value),
    predictor_label = pretty_predictor[as.character(predictor)]) %>%
  arrange(y_order) %>%
  mutate(
    display_term = fct_reorder(display_term, -y_order))


p1 <- plot_data %>%
  ggplot(aes(x = estimate, y = display_term, alpha = FDR < 0.05)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1, color = "#d62728") +
  geom_point(aes(shape = reference), size = 2, fill = "white", color = "#d62728") +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 21), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
  facet_nested(
    rows = vars(factor(predictor_label, levels=unique(plot_data$predictor_label))),
    scales = "free_y",
    space = "free_y"
  ) +
  labs(title = "Demographics", x = "Effect size", y = NULL,  alpha = "FDR < 5%") +
  theme_classic(base_size = 10) +
  theme(
    # Keep outer strip styling as fallback if needed
    strip.text.y.right = element_text(angle = 0, hjust = 0.5, face = "bold", size = 10),
    strip.background = element_rect(fill = "antiquewhite2", color = "black", linewidth = 0.8),
    panel.grid.major.x = element_line(linewidth = 0.1),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )
p1



