

library(tidyverse)
library(forcats)
library(ggtext)
library(paletteer)

fields <- data.table::fread("data/field.tsv")

df_r2 <- bind_rows(readRDS("data/aov_labs.rds") %>% mutate(type = "Biochemistry") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("data/aov_counts.rds") %>% mutate(type = "Cell_counts") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("data/aov_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("data/aov_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen)
  ) %>%
  mutate(color_var = case_when(type == "Proteomics-Olink" ~ "#76B041",
                               type == "Metabolomics-NMR" ~ "#2374AB",
                               type == "Cell_counts" ~ "#8F3985",
                               type == "Biochemistry" ~ "#E85F5C")) %>%
  mutate(title = case_when(title == "White blood cell (leukocyte) count" ~ "Leukocyte count",
                           title == "Phospholipids to Total Lipids in Small HDL percentage" ~ "Phosphlipid ratio SHDL",
                           title == "Cholesterol to Total Lipids in Very Large HDL percentage" ~ "Cholesterol ratio VLHDL",
                           title == "Phospholipids to Total Lipids in Very Large HDL percentage" ~ "Phospholipids ratio VLHDL",
                           title == "Cholesterol to Total Lipids in Small HDL percentage" ~ "Cholesterol ratio SHLD",
                           title == "Spectrometer-corrected alanine" ~ "Alanine",
                           TRUE ~ title))

df_top <- df_r2 %>%
  group_by(type, phen) %>%
  filter(any(term == "time_day" & p.value < 0.05)) %>%
  ungroup() %>%
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

plot_bars_v <- df_top %>%
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

  ggplot(aes(y = title, x = t_r2, fill = model)) +
  geom_col(width = 1) +
  facet_wrap(~facet_html, scales = "free_y", nrow = 20, ncol = 2) +
  scale_fill_paletteer_d("rcartocolor::Temps") +
  labs(fill = "Covariate", y = "", x = "R2") +
  guides(fill = guide_legend(reverse = TRUE, ncol = 3)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        strip.background = element_blank(),
        strip.text  = element_markdown(size = 10, hjust = 0),
        strip.placement = "inside",
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0.5, "lines"),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12)
  )


# library(cowplot)
#
# right <- plot_grid(plot_bars_h, plot_round, nrow = 2, rel_heights = c(0.55, 0.45),
#                    labels = c("B", "C"))
#
# plot_final <- plot_grid(plot_bars_v, right, ncol = 2, labels = c("A", ""))
#
# ggsave("plot_1.png", plot_final, width = 5.2, height = 5.2)

ggsave("plots/plot_vars_h.png", plot_bars_v, width = 4, height = 7)

