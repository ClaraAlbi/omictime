library(stringr)
library(tidyr)
library(dplyr)
library(purrr)
install.packages("cowplot")
install.packages("ggpmisc")
library(cowplot)
library(ggplot2)

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  mutate(i = case_when(i == 0  ~ "i0: Initial assessment \n(2006-2010)",
                       i == 2  ~ "i2: Imaging \n(2014+)",
                       i == 3  ~ "i3: First repeat imaging \n(2019+)"),
         i = factor(i, levels = c("i0: Initial assessment \n(2006-2010)", "i2: Imaging \n(2014+)", "i3: First repeat imaging \n(2019+)"))) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)))

df$res <- residuals(lm(pred_mean ~ time_day, data = df))

formula <- y ~ x

pr <- df %>%
  ggplot(aes(x = time_day, y = pred_mean)) +
  geom_point(alpha = 0.7, size = 1.5, color = "#76B041") +
  geom_smooth(method = "lm", color = "red", size = 1.2, se = FALSE, formula = formula) +
  facet_grid(~i, ) +
  ggpmisc::stat_poly_eq(
    aes(label = paste(after_stat(rr.label), sep = "*\", \"*")),
    formula = formula,
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 4,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    aes(label = paste("italic(n) ==", after_stat(n))),
    formula = formula,
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 4,
    color = "black"
  ) +
  scale_y_continuous(breaks = c(10, 15, 20), limits = c(9, 20)) +
  scale_x_continuous(breaks = c(10, 15, 20), limits = c(9, 20)) +
  labs(
    x = "Recorded time of day",
    y = "Predicted Proteomic Time"
  ) +
  theme_classic(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold", hjust = 0),
    axis.title = element_text(face = "bold")
  )

ggsave("plots/F3_internal_olink.png", pr, width = 6, height = 3)

##### NMR

df_nmr <- readRDS("/mnt/project/nmr_int_replication.rds") %>%
  filter(i == 1) %>%
  mutate(i = case_when(i == 0  ~ "i0: Initial assessment \n(2006-2010)",
                       i == 1 ~ "i1: First repeat assessment \n(2012-13)"),
         i = factor(i, levels = c("i0: Initial assessment \n(2006-2010)", "i1: First repeat assessment \n(2012-13)"))) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)))

df_nmr$res <- residuals(lm(pred_mean ~ time_day, data = df_nmr))

formula <- y ~ x

pr_nmr <- df_nmr %>%
  ggplot(aes(x = time_day, y = pred_mean)) +
  geom_point(alpha = 0.7, size = 1.5, color = "#2374AB") +
  geom_smooth(method = "lm", color = "red", size = 1.2, se = FALSE, formula = formula) +
  facet_grid(~i) +
  ggpmisc::stat_poly_eq(
    aes(label = paste(after_stat(rr.label), sep = "*\", \"*")),
    formula = formula,
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 4,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    aes(label = paste("italic(n) ==", after_stat(n))),
    formula = formula,
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 4,
    color = "black"
  ) +
  scale_y_continuous(breaks = c(10, 15, 20), limits = c(9, 20)) +
  scale_x_continuous(breaks = c(10, 15, 20), limits = c(9, 20)) +
  labs(
    x = "Recorded time of day",
    y = "Predicted Metabolic Time"
  ) +
  theme_classic(base_size =9) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold", hjust = 0),
    axis.title = element_text(face = "bold")
  )

ggsave("plots/F3_internal_nmr.png", pr_nmr, width = 3, height = 3)

p_comb <- cowplot::plot_grid(pr, pr_nmr, rel_widths = c(0.7, 0.3))

ggsave("plots/F3_internal.png", p_comb, width = 10, height = 3)

p_f <- plot_grid(pl, p_comb, nrow = 2)

ggsave("plots/F3.png", p_f, width = 10, height = 7)

