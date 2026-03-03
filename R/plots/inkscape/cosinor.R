library(tidyverse)
library(ggplot2)

# ---- Parameters ----
mesor <- 0
amplitude <- 1
acrophase <- 20     # peak at 6h
period <- 24

# ---- Generate time ----
time <- seq(0, 24, length.out = 500)

# ---- Cosinor equation ----
y <- mesor + amplitude * cos(2 * pi * (time - acrophase) / period)
y2 <- mesor + (amplitude + 0.01) * cos(2 * pi * (time - acrophase -0.5) / period)

is <- sample(which(time > 9 & time < 20), 20)

y3 <- y[is] + rnorm(length(is), sd = 0.3)

data <- data.frame(time, y, y2) %>%
  pivot_longer(y:y2)

data_p <- data.frame(time = time[is], name = "y2", value = y3)

# ---- Plot ----
p <- ggplot(data, aes(x = time, y = value, color = name)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = data_p, size = 0.2)  +
  #geom_vline(xintercept = 6, linetype = "dashed", color = "darkgray") +

  #geom_vline(xintercept = 18, linetype = "dashed", color = "darkgray") +

  scale_x_continuous(
    limits = c(0, 24),
    breaks = seq(0, 24, by = 6),
    expand = c(0, 0)
  ) +
  annotate("rect",
           xmin = 9, xmax = 20,
           ymin = -Inf, ymax = Inf,
           alpha = 0.15) +
  scale_color_manual(
    values = c(
      "y"  = "darkgray",
      "y2"  = "#2374AB"),
    labels = c("CR", "UKB")) +

  labs(
    x = "Time of day",
    y = "Amplitude", color = "Data"
  ) +
  theme_classic(base_size = 8) +
  theme(panel.background = element_blank(),
        legend.text = element_text(size = 4),
    axis.title =element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.5)
  )

p

ggsave("plots/cosinor_plot.pdf",
       plot = p,
       width = 1.7,
       height = 0.7,
       units = "in")



p_ckb <- ggplot(data, aes(x = time, y = y)) +
  geom_line(linewidth = 0.5) +

  scale_x_continuous(
    limits = c(0, 48),
    breaks = seq(0, 24, by = 12),
    expand = c(0, 0)
  ) +
  annotate("rect",
           xmin = 8, xmax = 21,
           ymin = -Inf, ymax = Inf,
           alpha = 0.15) +

  labs(title = "CKB") +
  theme_classic(base_size = 6) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 8),
          axis.title =element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.5)
  )

p_ckb



p_fg <- ggplot(data, aes(x = time, y = y)) +
  geom_line(linewidth = 0.5) +

  scale_x_continuous(
    limits = c(0, 48),
    breaks = seq(0, 24, by = 12),
    expand = c(0, 0)
  ) +
  annotate("rect",
           xmin = 11, xmax = 13,
           ymin = -Inf, ymax = Inf,
           alpha = 0.15) +

  labs(title = "FinnGen") +
  theme_classic(base_size =6) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 8),
          axis.title =element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.5)
  )
p_fg



time <- seq(0, 48, length.out = 500)

# ---- Cosinor equation ----
y <- mesor + amplitude * cos(2 * pi * (time - acrophase) / period)

data2 <- data.frame(time, y)

p_tre <- ggplot(data2, aes(x = time, y = y)) +
  geom_line(linewidth = 0.5) +

  scale_x_continuous(
    limits = c(0, 48),
    breaks = seq(0, 48, by = 12),
    expand = c(0, 0)
  ) +
  annotate("rect",
           xmin = 0, xmax = 40,
           ymin = -Inf, ymax = Inf,
           alpha = 0.15) +

  labs(title = "TREASURE") +
  theme_classic(base_size = 6) +
  theme(panel.background = element_blank(),
        plot.title = element_text(size = 8),
        axis.title =element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.5)
  )
p_tre


library(patchwork)

p_pout <- (p_ckb / p_fg / p_tre) +
  plot_layout(heights = c(1, 1, 1)) &
  theme(plot.margin = margin(1, 2, 1, 1))

ggsave("plots/cosinor_plot_ext.pdf",
       plot = p_pout,
       width = 1.1,
       height = 1.8,
       units = "in")
