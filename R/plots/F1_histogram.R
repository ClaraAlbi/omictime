
library(dplyr)
library(tidyr)
library(glue)
library(ggplot2)

time <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  filter(time_day > 9 & time_day < 20)

light_band <- data.frame(
  xmin = 5.4,
  xmax = 20.5,
  ymin = -Inf,
  ymax = Inf
)

night_band <- data.frame(
  xmin = c(0, 20.5),
  xmax = c(5.4, 24),
  ymin = -Inf,
  ymax = Inf
)

p_hist <- time %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60, color = "black", fill = "#355F71") +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  scale_x_continuous(limits = c(0, 24), n.breaks = 23) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

ggsave("plots/plot_histogram_i0_1.png", p_hist, width = 8, height = 8)
