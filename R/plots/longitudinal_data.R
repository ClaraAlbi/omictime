library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
install.packages("paletteer")
install.packages("ggmisc")
install.packages("forcats")

data <- readRDS("data_share/predictions_internal_time_updated.rds") %>%
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

summary(lm(time_day ~ pred_mean, data %>% filter(time_day >= 9 & time_day <= 20)))

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





## USEFUL



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


n <- 12
r <- sprintf("%.2f", summary(fit)$r.squared)

# Plot observed vs fitted, with residual lines
p_long <- ggplot(data_plot, aes(x = time_day, y = pred_mean)) +
  geom_point(alpha = 0.8, aes(color = participantid), size = 2) +
  # residual lines (vertical from fitted to observed)
  #geom_segment(
  #  aes(xend = time_day, y = fitted, yend = pred_mean),
  #  color = "grey60", alpha = 0.6
  #) +
  ggpmisc::stat_poly_eq(
    aes(label = paste("italic(R^2) ==", after_stat(r))),
    formula = formula,
    parse = TRUE,
    label.x = 0.02,
    label.y = 1,
    size = 4,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    aes(label = paste("italic(n) ==", after_stat(n))),
    formula = formula,
    parse = TRUE,
    label.x = 0.02,
    label.y = 0.90,
    size = 4,
    color = "black"
  ) +
  geom_line(aes(y = fitted), color = "red", linewidth = 0.7) +
  labs(
    title = "Longitudinal data",
    x = "Recorded time of day",
    color = "Participant",
    y = "Predicted Proteomic time"
  ) +
  paletteer::scale_color_paletteer_d("dichromat::DarkRedtoBlue_12", direction = -1) +
  guides(
    color = guide_legend(
      byrow = TRUE,
      keyheight = unit(10, "pt"),
      keywidth  = unit(10, "pt"),
      default.unit = "pt",
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.background = element_blank(),
    #legend.spacing.y = unit(0, "pt"),
    strip.text = element_text(size = 10, face = "bold", hjust = 0),
    #legend.text = element_text(size = 8),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -10, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0),
    #legend.box.spacing = unit(0, "pt"),
    axis.title = element_text(face = "bold"), legend.position = "right"
  )


ggsave("plots/F3_long.png", p_long, width = 7, height = 3)


p1 <- ggplot(data_plot, aes(x = time_day, y = pred_mean)) +
  geom_point(aes(color = participantid)) +
  theme_classic()

# blank white plot
blank <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "white", color = NA))

# combine: main plot on left, white filler on right
p_ext <- plot_grid(p_long, blank, ncol = 2, rel_widths = c(0.7, 0.3))


library(forcats)

data_plot %>%
  group_by(participantid) %>%
  count()

d_grouped <-data_plot %>%
  group_by(participantid) %>% mutate(m_res = mean(resid))

broom::tidy(lm(resid ~ participantid, data = d_grouped)) %>%
  mutate(p_adj = p.adjust(p.value))

p_c <- d_grouped %>%
  ggplot(aes(x = fct_reorder(as.factor(participantid), m_res), y = resid, fill = participantid)) +
  geom_hline(yintercept = 0, linetype = 2)  +
  geom_boxplot() +
  labs(y = "Acceleration", x = "Participant ID", title = "Longitudinal") +
  paletteer::scale_fill_paletteer_d("dichromat::DarkRedtoBlue_12", direction = -1) +
  theme_classic(base_size = 14) +
  theme(legend.position = 'none',
        plot.title   = element_text(face = "bold", size = 16),
        axis.title   = element_text(face = "bold"))

ggsave("plots/F3_long_acc.png", p_c, width = 7, height = 3)



