
library(ggplot2)

df <- readRDS("olink_int_replication.rds") %>%
  mutate(gap = time_day - pred_lasso) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  mutate(across(c(y, m), as.numeric)) %>%
  filter(n == 3)

df$res <- residuals(lm(time_day ~ pred_lasso, data = df))
df$r2  <- residuals(lm(res ~ pred_lasso, data = df))

metrics <- df %>%
  group_by(eid) %>%
  summarize(
    mean_value = mean(res, na.rm = TRUE),
    mean_time = mean(time_day, na.rm = TRUE),
    sd_value   = sd(res,   na.rm = TRUE),
    range_value= max(res, na.rm = TRUE) - min(res, na.rm = TRUE)
  ) %>%
  mutate(
    cv_value = sd_value / mean_value
  )

top_30 <- metrics %>%
  arrange(desc(range_value)) %>%
  slice(1:50)

metrics %>%
  ggplot(aes(x = mean_time, y = mean_value)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value))

# 4) Spaghetti plot for a random subset of 50 individuals
set.seed(123)  # for reproducibility
sample_ids <- sample(unique(df$eid), 50)
df %>%
  left_join(metrics) %>%
  #filter(eid %in% top_30$eid) %>%
  ggplot(aes(x = y, y = gap, group = eid, color = range_value)) +
  geom_line(alpha = 0.5) +
  labs(
    x = "Timepoint",
    y = "Biomarker Value",
    title = "Individual Trajectories (n = 50)"
  ) +
  theme_minimal()
