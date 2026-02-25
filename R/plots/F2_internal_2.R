install.packages("patchwork")
library(patchwork)
install.packages("ggpmisc")
library(ggpmisc)


library(stringr)
library(tidyr)
library(dplyr)
library(purrr)


df <- readRDS("/mnt/project/biomarkers_res/olink_int_replication_v2.rds") %>%
  filter(i != 0) %>%
  bind_rows(readRDS("/mnt/project/olink_int_panels14.rds") %>% mutate(i = 0)) %>%
  filter(!is.na(time_day)) %>%
  mutate(i = case_when(i == 0  ~ "i0: Initial visit \n(2006-2010)",
                       i == 2  ~ "i2: Imaging visit \n(2014+)",
                       i == 3  ~ "i3: Repeat imaging visit \n(2019+)"),
         i = factor(i, levels = c("i0: Initial visit \n(2006-2010)", "i2: Imaging visit \n(2014+)", "i3: Repeat imaging visit \n(2019+)"))) %>%
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
  scale_y_continuous(breaks = c(10, 15, 20)) +
  scale_x_continuous(breaks = c(10, 15, 20)) +
  coord_cartesian(xlim = c(9, 20), ylim = c(9, 20)) +
  labs(
    x = "Recorded time",
    y = "Predicted proteomic time"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold", hjust = 0)
  )

pr
ggsave("plots/F3_internal_olink.png", pr, width = 6, height = 3)

##### NMR

df_nmr <- readRDS("/mnt/project/nmr_int_replication.rds") %>%
  filter(i == 1) %>%
  mutate(i = case_when(i == 0  ~ "i0: Initial visit \n(2006-2010)",
                       i == 1 ~ "i1: Repeat visit \n(2012-13)"),
         i = factor(i, levels = c("i0: Initial visit \n(2006-2010)", "i1: Repeat visit \n(2012-13)"))) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)))

df_nmr$res <- residuals(lm(pred_mean ~ time_day, data = df_nmr))


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
  scale_y_continuous(breaks = c(10, 15, 20)) +
  scale_x_continuous(breaks = c(10, 15, 20)) +
  coord_cartesian(xlim = c(9, 20), ylim = c(9, 20)) +
  labs(
    x = "Recorded time",
    y = "Predicted metabolomic time"
  ) +
  theme_classic(base_size =10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold", hjust = 0)
  )

ggsave("plots/F3_internal_nmr.png", pr_nmr, width = 3, height = 3)



# blank white plot
blank <- ggplot() + labs(title = "   FinnGen",
                         y = "Predicted proteomic time",
                         x = "Recorded time") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 10),
        strip.background = element_blank())

blank2 <- ggplot() + labs(title = "   CKB",
                          y = "Predicted proteomic time",
                          x = "Recorded time") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.title = element_text(size = 10),
        strip.background = element_blank())


p_ext <- pl / (pr + pr_nmr) / (blank2 + blank + p_long)

p_ext <-
  pl /
  (pr + pr_nmr + plot_layout(widths = c(0.8, 0.2))) /
  (blank2 + blank + p_long + plot_layout(widths = c(0.7, 0.4, 1.3))) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold")
  )

ggsave("plots/Figure_3.png", p_ext, width = 10, height = 10, dpi = 320)
ggsave("plots/Figure_3.pdf", p_ext, width = 10, height = 10, dpi = 320)


