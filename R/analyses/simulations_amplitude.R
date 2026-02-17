# Circadian protein simulation (R)
# Creates:
# 1) Reference (Ref)
# 2) Advanced phase individual (same amplitude)
# 3) Lower amplitude individual (same phase)
# Then plots all three.

# Required packages: install.packages(c('tidyverse'))

library(tidyverse)

set.seed(42)

# Simulation parameters
duration_h <- 24          # total hours
dt <- 0.1                 # time resolution (hours)
time <- seq(0, duration_h, by = dt)

# Circadian cosine model
circadian <- function(t, amplitude = 1, phase_hours = 0, baseline = 0.5, period = 24){
  omega <- 2 * pi / period
  baseline + amplitude * cos(omega * (t - phase_hours))
}

# ----------------------------
# 1) Reference
# ----------------------------
ref_amp <- 1.0
ref_phase <- 0

ref_signal <- circadian(time, amplitude = ref_amp, phase_hours = ref_phase)

# ----------------------------
# 2) Advanced phase ONLY
# (same amplitude, earlier peak)
# ----------------------------
adv_amp <- 1.0
adv_phase <- 2   # 4-hour advance

adv_signal <- circadian(time, amplitude = adv_amp, phase_hours = adv_phase)

# ----------------------------
# 3) Lower amplitude ONLY
# (same phase, reduced amplitude)
# ----------------------------
low_amp <- 0.5
low_phase <- 0

low_signal <- circadian(time, amplitude = low_amp, phase_hours = low_phase)

# Optional measurement noise
noise_sd <- 0.03

ref_noisy <- ref_signal + rnorm(length(time), sd = noise_sd)
adv_noisy <- adv_signal + rnorm(length(time), sd = noise_sd)
low_noisy <- low_signal + rnorm(length(time), sd = noise_sd)

# Build tidy dataframe safely

df_ref <- tibble(
  time = time,
  value = ref_noisy,
  true = ref_signal,
  subject = "Ref"
)

df_adv <- tibble(
  time = time,
  value = adv_noisy,
  true = adv_signal,
  subject = "Advanced phase"
)

df_low <- tibble(
  time = time,
  value = low_noisy,
  true = low_signal,
  subject = "Lower amplitude"
)

df <- bind_rows(df_ref, df_adv, df_low)

# Sample every 2 hours to mimic discrete sampling
sample_times <- seq(0, duration_h, by = 2)
sampled <- df %>% filter(time %in% sample_times)

# LOESS smoothing
smoothed <- df %>%
  group_by(subject) %>%
  arrange(time) %>%
  mutate(loess = predict(loess(value ~ time, span = 0.1), time)) %>%
  ungroup()


# Plot
p <- smoothed %>%
  ggplot() +
  geom_line(data = smoothed,
            aes(x = time, y = loess, color = subject),
            linewidth = 1.2) +
  geom_point(data = sampled,
             aes(x = time, y = value, color = subject),
             size = 1.8, alpha = 0.85) +
  geom_vline(xintercept = c(6, 18)) +
  scale_color_manual(values = c(
    "Ref" = "#1B3A8A",              # deep blue
    "Advanced phase" = "#C1121F",   # strong red
    "Lower amplitude" = "#2A9D8F"   # teal
  ), breaks = c("Ref", "Advanced phase", "Lower amplitude")) +
  labs(x = "Time of day",
       y = "Protein levels",
       color = NULL) +
  scale_x_continuous(breaks = seq(0, duration_h, by = 4)) +
  scale_y_continuous(limits = c(-0.6, 2)) +
  guides(color = guide_legend(
    nrow = 3
  )) +
    theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.5, 0.95),   # inside plot (center top)
    legend.justification = c(0.5, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(color = "black", linewidth = 0.3),
    legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

p


ggsave("plots/example_adv_low_phase.png", p)
