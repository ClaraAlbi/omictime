library(ggplot2)

set.seed(123)

# ---- Parameters ----
n <- 100
rho <- 0.7

# ---- Generate correlated standard normals ----
x0 <- rnorm(n)
y0 <- rho * x0 + sqrt(1 - rho^2) * rnorm(n)

# ---- Rescale both to 9–20 ----
rescale_9_20 <- function(z) {
  9 + (z - min(z)) / (max(z) - min(z)) * (20 - 9)
}

x <- rescale_9_20(x0)
y <- rescale_9_20(y0)

df <- data.frame(x, y)

# Check correlation (optional)
cor(df$x, df$y)

# ---- Plot ----
p <- ggplot(df, aes(x, y)) +

  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "gray") +

  scale_x_continuous(
    limits = c(8, 21),
    breaks = seq(8, 20, by = 4),
    expand = c(0, 0)
  ) +

  scale_y_continuous(
    limits = c(8, 21),
    breaks = seq(8, 20, by = 4),
    expand = c(0, 0)
  ) +

  coord_fixed() +   # keeps square aspect ratio

  theme_classic(base_size = 8) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.title = element_text(face = "bold")
  ) +

  labs(
    x = "Clock time",
    y = "Predicted time"
  )

ggsave("plots/pred_cartoon.pdf",
       plot = p,
       width = 1.1,
       height = 1.8,
       units = "in")

