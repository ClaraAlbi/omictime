

library(tidyverse)
library(forcats)
library(ggtext)
library(paletteer)

fields <- data.table::fread("data/field.tsv")

df_r2 <- readRDS("data/combined_variance.rds")

df_top <- df_r2 %>%
  filter(term == "time_day") %>%
  mutate(pval = p.adjust(p.value)) %>%
  filter(pval < 0.05) %>%
  distinct(phen, .keep_all = TRUE) %>%
  arrange(desc(pr2)) %>%
  slice_head(n = 30) %>%
  select(phen) %>%
  inner_join(df_r2, by = "phen") %>%
  filter(term != "Residuals") %>%
  mutate(model = case_when(!term %in% c("bmi", "sex", "age_recruitment", "fasting",
                                        "smoking", "time_day") ~ "technical",
                           TRUE ~ term)) %>%
  group_by(phen, model, color_var, type, title) %>%
  summarise(t_r2 = sum(pr2))

# trait_order <- df_top %>%
#   filter(model == "time_day") %>%
#   arrange(desc(t_r2)) %>%
#   pull(title)


facet_levels <- df_top %>%
  filter(model == "time_day") %>%
  arrange(desc(t_r2)) %>%
  # recreate the exact HTML string you’ll use below
  mutate(f_html = sprintf("<span style='color:%s'>%s</span>", color_var, title)) %>%
  pull(f_html)


time_lookup <- df_top %>%
  filter(model == "time_day") %>%
  group_by(phen) %>%
  slice_max(t_r2, n = 1, with_ties = FALSE) %>%   # or slice_head(n = 1)
  ungroup() %>%
  transmute(phen, time_var = round(t_r2, 3)*100)

#plot_bars_v <-
df_plot <- df_top %>%
  mutate(
    f_html = sprintf("<span style='color:%s'>%s</span>", color_var, title),
    facet_html = factor(f_html, levels = facet_levels),
    model = case_when(model == "age_recruitment" ~ "Age",
                      model == "time_day" ~ "Time of day",
                      model == "technical" ~ "Technical",
                      model == "smoking" ~ "Smoking",
                      model == "sex" ~ "Sex",
                      model == "bmi" ~ "BMI",
                      model == "fasting" ~ "Fasting"),
    model = factor(
      model,
      levels = c("Fasting", "BMI", "Smoking", "Sex", "Age", "Technical",  "Time of day")
    )
  ) %>%
  left_join(time_lookup)

df_lab <- df_plot %>%
  group_by(facet_html, title) %>%
  summarise(
    x_end    = sum(t_r2, na.rm = TRUE),
    time_var = dplyr::first(time_var),    # common per bar
    .groups  = "drop"
  ) %>%
  filter(!is.na(time_var)) %>%
  mutate(
    label = paste0(scales::number(time_var, accuracy = 1), "%"),
    x_lab =  0.01                   # nudge to the right (2%)
  )
plot_bars_v <- ggplot(df_plot, aes(y = title, x = t_r2, fill = model)) +
  geom_col(width = 1) +
  geom_text(
    data = df_lab,
    aes(x = x_lab, y = title, label = label),
    inherit.aes = FALSE,
    hjust = 0, size = 3.3
  ) +
  facet_wrap(~facet_html, scales = "free_y", nrow = 20, ncol = 2) +
  scale_fill_paletteer_d("rcartocolor::Temps") +
  labs(fill = "Covariate", y = "", x = "R2") +
  guides(fill = guide_legend(reverse = TRUE, ncol = 3)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.background = element_blank(),
        strip.text  = element_markdown(size = 14, hjust = 0),
        strip.placement = "inside",
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
  )


# library(cowplot)
#
# right <- plot_grid(plot_bars_h, plot_round, nrow = 2, rel_heights = c(0.55, 0.45),
#                    labels = c("B", "C"))
#
# plot_final <- plot_grid(plot_bars_v, right, ncol = 2, labels = c("A", ""))
#
# ggsave("plot_1.png", plot_final, width = 5.2, height = 5.2)

ggsave("plots/plot_vars_h.png", plot_bars_v, width = 6, height = 10)

