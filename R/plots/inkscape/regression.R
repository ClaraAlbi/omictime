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
p

ggsave("plots/pred_cartoon.pdf",
       plot = p,
       width = 1.1,
       height = 1.8,
       units = "in")




#### CA PLOT


# ---- Parameters ----
n <- 500
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


mod <- lm(y ~ x, data = df)
df$pred <- predict(mod)
df$col <- residuals(mod)



top2 <- df[c(35, 173),]

# ---- Plot ----
p <- ggplot(df, aes(x, y, color = col)) +

  geom_point(size = 1, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "red") +

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
  geom_point(
    data  = top2,
    aes(x = x, y = y),
    color = "black",
    size  = 1.5, shape = 21,
  ) +
  geom_segment(
    data = top2,
    aes(
      x = x,
      y = y,
      xend = x,
      yend = pred
    ),
    arrow = arrow(length = unit(0.15, "cm")),
    color = "black",
    linewidth = 0.5
  ) +

  coord_fixed() +   # keeps square aspect ratio

  paletteer::scale_color_paletteer_c("ggthemes::Orange-Blue Diverging",
                                     direction = -1,
                                     limits = c(-5, 5)) +
  guides(
    colour = guide_colourbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = 2,
      barheight      = 1.2,
      reverse = TRUE
    )
  ) +

  labs(
    x = "Recorded time",
    y = "Circadian phase",
    color = "Circadian Acceleration"
  ) +
  theme_classic(base_size = 8) +
  theme(axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        legend.text = element_text(size = 4),
        legend.position   = "top",
        legend.margin = margin(0,0,0,0))
p

ggsave("plots/pred_cartoon_CA.pdf",
       plot = p,
       width = 1.1,
       height = 1.8,
       units = "in")



