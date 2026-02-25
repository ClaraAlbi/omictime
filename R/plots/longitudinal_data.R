library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

install.packages("broom")

install.packages("readxl")
library(readxl)
library(lubridate)
install.packages("paletteer")

install.packages("ggmisc")
library(ggmisc)
install.packages("forcats")
library(forcats)

### changes
data_x <- readRDS("data_share/predictions_internal_time_updated.rds") %>%
  mutate(visitnumber = as.numeric(visitnumber),
         participantid = as.numeric(participantid))

data_2 <- data_x %>% filter(participantid == 2) %>%
  select(eid, starts_with("pred"), visitnumber, visitday)

data <- read_xlsx("data_share/modified_dates_long.xlsx") %>%
  mutate(visitnumber = as.numeric(str_remove(visit_num, "V")),
         time_day = hour(ymd_hms(time_h)) + minute(ymd_hms(time_h))/60,
         participantid = 2) %>% select(-date, -x, -time_h, -visit_num) %>%
  left_join(data_2) %>%
  bind_rows(., data_x %>% filter(participantid != 2)) %>%
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


### FIGURE GRID

summary(lm(time_day ~ pred_mean, data %>% filter(time_day > 9 & time_day < 20)))

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


### HISTOGRAM
p_hist <- data %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

ggsave("plots/time_histogram_TREASURE.png", p_hist, width = 8, height = 8)



## USEFUL



fit <- lm(pred_mean ~ sin(2*pi*time_extended/24) + cos(2*pi*time_extended/24),
          data = data)

summary(fit)

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
  geom_point(alpha = 0.8, aes(color = as.factor(participantid)), size = 1.5) +
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
    title = "TREASURE",
    x = "Recorded time",
    color = "PID",
    y = "Predicted proteomic time"
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
    legend.position = "right"
  )

p_long
ggsave("plots/F3_long.png", p_long, width = 7, height = 3)


library(forcats)

data_plot %>%
  group_by(participantid) %>%
  count()

d_grouped <-data_plot %>%
  group_by(participantid) %>% mutate(m_res = mean(resid))

tests <- broom::tidy(lm(resid ~ 0 + participantid, data = d_grouped)) %>%
  mutate(p_adj = p.adjust(p.value)) %>%
  filter(p_adj < 0.05)

stars <- d_grouped %>%
  mutate(star = case_when(participantid %in% c(2, 5, 6, 10, 12) ~ "*",
         TRUE ~ "")) %>%
  group_by(participantid, star) %>%
  mutate(y = max(resid) + 0.3)

p_c <- d_grouped %>%
  ggplot(aes(x = fct_reorder(as.factor(participantid), m_res), y = resid, fill = m_res)) +
  geom_hline(yintercept = 0, linetype = 2)  +
  geom_boxplot() +
  geom_text(
    data = stars,
    aes(x = as.factor(participantid), y = y, label = star),
    size = 6
  ) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange-Blue Diverging",
                                     direction = -1,
                                     limits = c(-5, 5)) +
  labs(y = "CA", x = "Participant ID", title = "Samples within 40h (TREASURE)") +
  theme_classic(base_size = 12) +
  guides(
    fill = guide_colourbar(
      title = "Mean CA",
      title.position = "top",
      #direction = "horizontal",
      title.hjust    = 0.5,
      barwidth       = 1.5,
      barheight      = 6,
      reverse = TRUE
    )
  ) +
  theme(legend.position = "left",
        #legend.position = c(0.5, 0.95),   # x, y (0–1 scale)
        #legend.justification = c(1, 1),    # anchor legend box
        #legend.background = element_rect(fill = "white", colour = "black"),
        plot.title = element_text(face = "bold", size = 13),
        axis.title = element_text(size = 10))

p_c

ggsave("plots/F3_long_acc.png", p_c, width = 7, height = 3)

### ICC


m2 <- lmer(
  resid ~ time_day + (1 | participantid),
  data = data_plot
)

p2 <- performance::icc(m2)


