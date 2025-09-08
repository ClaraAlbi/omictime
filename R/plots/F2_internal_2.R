library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("cowplot")
install.packages("ggpmisc")

df <- readRDS("olink_int_replication.rds") %>%
  mutate(i = case_when(i == 0  ~ "Initial assessment \n(2006-2010)",
                       i == 2  ~ "Imaging \n(2014+)",
                       i == 3  ~ "First repeat imaging \n(2019+)"),
         i = factor(i, levels = c("Initial assessment \n(2006-2010)", "Imaging \n(2014+)", "First repeat imaging \n(2019+)"))) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T)

df$res <- residuals(lm(pred_lasso ~ time_day, data = df))

formula <- y ~ x

pr <- df %>%
  ggplot(aes(x = time_day, y = pred_lasso)) +
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
  labs(
    x = "Recorded Time of day",
    y = "Predicted Proteomic Time"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold", hjust = 0),
    axis.title = element_text(face = "bold")
  )

ggsave("plots/F3_internal_olink.png", pr, width = 6, height = 3)

##### NMR

df <- df_nmr %>%
  mutate(i = case_when(i == 0  ~ "Initial assessment \n(2006-2010)",
                       i == 1 ~ "First repeat assessment \n(2012-13)"),
         i = factor(i, levels = c("Initial assessment \n(2006-2010)", "First repeat assessment \n(2012-13)"))) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T)

df$res <- residuals(lm(pred_lasso ~ time_day, data = df))

formula <- y ~ x

pr <- df %>%
  ggplot(aes(x = time_day, y = pred_lgb)) +
  geom_point(alpha = 0.7, size = 1.5, color = "#2374AB") +
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
  labs(
    x = "Recorded Time of day",
    y = "Predicted Proteomic Time"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold", hjust = 0),
    axis.title = element_text(face = "bold")
  )

ggsave("plots/F3_internal_nmr.png", pr, width = 6, height = 3)

