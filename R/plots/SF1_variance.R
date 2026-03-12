install.packages("ggtext")
library(ggtext)
install.packages("paletteer")
library(paletteer)

df_r2 <- readRDS("tables/variance_covariates.rds")

df_top <- df_r2 %>%
  mutate(pfdr = p.adjust(p.value_time_day)) %>%
  filter(pfdr < 0.05) %>%
  distinct(Label, .keep_all = TRUE) %>%
  filter(t_r2_time_day > 0.01) %>%
  mutate(color = case_when(Type == "Proteins" ~ "#76B041",
                           Type == "Metabolites"  ~ "#2374AB",
                           Type == "Cell counts" ~ "#8F3985",
                           Type == "Biochemistry" ~ "#E85F5C")) %>%
  pivot_longer(
    cols = starts_with("t_r2_"),
    names_to = "term",
    names_prefix = "t_r2_",
    values_to = "value"
  ) %>%
  filter(term != "Residuals") %>%
  mutate(model = case_when(!term %in% c("bmi", "sex", "age_recruitment", "fasting",
                                        "smoking","PCs","month_attending",  "time_day") ~ "technical",
                           TRUE ~ term)) %>%
  group_by(Biomarker, Label, Type, model, color) %>%
  summarise(t_r2 = sum(value, na.rm = T)) %>% ungroup() %>%
  mutate(
    f_html = sprintf("<span style='color:%s'>%s</span>", color, Label),
    #facet_html = factor(f_html, levels = facet_levels),
    model = case_when(model == "age_recruitment" ~ "Age",
                      model == "time_day" ~ "Time of day",
                      model == "technical" ~ "Technical",
                      model == "smoking" ~ "Smoking",
                      model == "sex" ~ "Sex",
                      model == "bmi" ~ "BMI",
                      model == "month_attending" ~ "Month",
                      model == "PCs"  ~ "20 genetic PCs",
                      model == "fasting" ~ "Fasting"),
    model = factor(
      model,
      levels = c("Technical", "Fasting", "BMI", "Smoking", "Sex", "Age", "20 genetic PCs", "Month",  "Time of day")
    )
  )

df_lab <-
  df_top %>% filter(model == "Time of day") %>%
  mutate(
    f_html = fct_reorder(f_html, t_r2 ,  .desc = TRUE),
    l = paste0(100*round(t_r2, 2), "%"),
    x_lab =  0.01
  )

df_top <- df_top %>%
  group_by(f_html) %>%
  mutate(order_r2 = t_r2[model == "Time of day"][1]) %>%
  ungroup() %>%
  mutate(f_html = fct_reorder(f_html, order_r2, .desc = TRUE))

plot_bars_v <- ggplot(df_top, aes(y = Label, x = t_r2, fill = model)) +
  geom_col(width = 1) +
  geom_text(
    data = df_lab,
    aes(label = l, x = -0.05, y = Label),
    inherit.aes = FALSE,
    hjust = 0, size = 3.3
  ) +
  facet_wrap(~f_html, scales = "free_y", ncol = 4, nrow = 35) +
  scale_fill_paletteer_d("trekcolors::tholian") +
  labs(fill = "Covariate", y = "", x = "R2") +
  scale_x_continuous(breaks = seq(0, 0.4, 0.1)) +
  guides(fill = guide_legend(reverse = TRUE, ncol = 3)) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 10),
        strip.background = element_blank(),
        strip.text  = element_markdown(size = 8, hjust = 0),
        strip.placement = "inside",
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)
  )


#plot_bars_v
ggsave("plots/plot_vars_h.png", plot_bars_v, width = 10, height = 12.5, dpi = 320)
