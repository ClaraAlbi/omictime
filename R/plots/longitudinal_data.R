library(tidyverse)

data <- readRDS("data/predictions_internal_time_updated.rds") %>%
  mutate(time_extended = time_day + 24 * (as.numeric(visitday) - 1))

data <- data %>%
  ungroup() %>%
  mutate(
    pred_scaled = (pred_mean - min(pred_mean, na.rm = TRUE)) /
      (max(pred_mean, na.rm = TRUE) - min(pred_mean, na.rm = TRUE)) *
      24
  )

data <- data %>%
  mutate(pred_scaled = pred_scaled %% 24)

grid_df <- expand_grid(
  participantid = unique(data$participantid),
  time_day = seq(0, 23, by = 1),
  visitday = as.factor(c(1,2))) %>%
  left_join(data %>% select(participantid, time_day, visitday) %>% mutate(has_sample = 1,
                                                                          time_day = round(time_day, 0)),
            by = c("participantid", "time_day", "visitday")) %>%
  mutate(has_sample = ifelse(is.na(has_sample), 0, has_sample))

ggplot(grid_df, aes(x = time_day, y = participantid, fill = factor(has_sample))) +
  geom_tile(color = "white") +
  scale_fill_manual(
    values = c("0" = "grey95", "1" = "#2b8cbe"),
    labels = c("No sample", "Has sample"),
    name = NULL
  ) +
  facet_grid(~visitday) +
  labs(
    x = "Time of day (hours)",
    y = "Participant"
    ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "top"
  )

data %>%
  ggplot(aes(x = time_day, y = pred_scaled, color = participantid)) +
  geom_point() +
  geom_smooth(se = F, ) +
  facet_grid(~visitday) +
  theme_minimal()

data %>%
  ggplot(aes(x = time_day, y = pred_scaled, color = participantid)) +
  geom_point() +
  geom_smooth(se = F, ) +
  theme_minimal()


harmonic_method <- function(formula, data, weights = NULL, ...) {
  # extract y and x names from the formula that ggplot supplies
  yvar <- all.vars(formula)[1]
  xvar <- all.vars(formula)[2]

  # construct harmonic regression dynamically
  f <- as.formula(paste0(
    yvar, " ~ sin(2*pi*", xvar, "/24) + cos(2*pi*", xvar, "/24)"
  ))

  lm(f, data = data)
}

ggplot(data, aes(x = time_day, y = pred_scaled)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = harmonic_method, se = FALSE, color = "red") +
  theme_minimal()

fit <- lm(pred_mean ~ sin(2*pi*time_extended/24) + cos(2*pi*time_extended/24),
          data = data)

# Add fitted values and residuals to data
data_plot <- data %>%
  ungroup() %>%
  mutate(
    fitted = fitted(fit),
    resid  = resid(fit),
    pred_extended = fitted + 24 * (as.numeric(visitday) - 1)
  )

# Plot observed vs fitted, with residual lines
ggplot(data_plot, aes(x = time_day, y = pred_mean)) +
  geom_point(alpha = 0.8, aes(color = participantid)) +
  # residual lines (vertical from fitted to observed)
  geom_segment(aes(xend = time_day, y = fitted, yend = pred_mean),
               color = "grey60", alpha = 0.6) +
  # fitted curve
  geom_line(aes(y = fitted), color = "red", linewidth = 0.5) +
  labs(
    title = "Longitudinal data",
    x = "Recorded time",
    y = "Mean proteomic predicted time",
    subtitle = sprintf("R² = %.2f", summary(fit)$r.squared)
  ) +
  theme_minimal()

ggplot(data_plot, aes(x = resid)) + geom_histogram() +
  theme_minimal()


data_plot %>%
  filter(resid > -10) %>%
  group_by(participantid) %>%
  mutate(mean_res = mean(resid, na.rm = TRUE)) %>%
  ggplot(aes(x = time_day, y = resid, color = visitday)) +
  geom_point() +
  geom_hline(aes(yintercept = mean_res), color = "black", linetype = "dashed") +
  facet_grid(~ participantid, scales = "free_y") +
  theme_minimal(base_size = 14) +
  scale_x_continuous(n.breaks = 3)

data_plot %>%
  mutate(time_day = round(time_day, 0)) %>%
  pivot_wider(names_from = visitday, values_from = resid, names_prefix = "visit") %>%
  rowwise() %>%
  mutate(mean_r = mean(c(visit1, visit2))) -> a%>%
  ggplot(aes(x = time_day, y = mean_r, color = participantid)) +
  geom_point()



