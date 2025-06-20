library(tidyverse)
library(cowplot)
library(ggpubr)

fields <- data.table::fread("data/field.tsv")

df_effects <- bind_rows(readRDS("data/effects_labs.rds") %>% mutate(type = "Biochemistry") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                        readRDS("data/effects_counts.rds") %>% mutate(type = "Cell_counts") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                        readRDS("data/effects_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("data/effects_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen)
  ) %>%
  mutate(color_var = case_when(type == "Proteomics-Olink" ~ "#76B041",
                               type == "Metabolomics-NMR" ~ "#2374AB",
                               type == "Cell_counts" ~ "#8F3985",
                               type == "Biochemistry" ~ "#E85F5C"),
         term = case_when(term == "(Intercept)" ~ "mesor",
                          term == "cos(2 * pi * time_day/24)" ~ "beta_cos1",
                          term == "sin(2 * pi * time_day/24)" ~ "beta_sin1")) %>%
  pivot_wider(id_cols = c(phen, color_var, type, title), values_from = c(estimate, p.value, std.error), names_from = term) %>%
  mutate(amplitude_24hfreq = sqrt(estimate_beta_cos1^2 + estimate_beta_sin1^2),
         acrophase_24hfreq = (atan2(estimate_beta_sin1, estimate_beta_cos1) / (2 * pi) * 24 + 24) %% 24,
         q = as.integer(round(acrophase_24hfreq, 0)))


df_r2 <- bind_rows(readRDS("data/aov_labs.rds") %>% mutate(type = "Biochemistry") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("data/aov_counts.rds") %>% mutate(type = "Cell_counts") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("data/aov_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("data/aov_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen))


light_band <- data.frame(
  xmin = 6,
  xmax = 20,
  ymin = -Inf,
  ymax = Inf
)

night_band <- data.frame(
  xmin = c(0, 20),
  xmax = c(6, 24),
  ymin = -Inf,
  ymax = Inf
)

df_top <- df_r2 %>%
  group_by(type, phen) %>%
  filter(any(term == "time_day" & p.value < 0.05)) %>%
  ungroup() %>%
  distinct(phen, .keep_all = TRUE) %>%
  arrange(desc(pr2)) %>%
  slice_head(n = 20) %>%
  select(phen)

labes <- full_join(df_effects[df_effects$amplitude_24hfreq > 0.4,], df_effects[df_effects$phen %in% df_top$phen,]) %>%
  mutate(title = case_when(title == "White blood cell (leukocyte) count" ~ "Leukocyte count",
                           title == "Phospholipids to Total Lipids in Small HDL percentage" ~ "Phosphlipid ratio",
                           title == "Cholesterol to Total Lipids in Very Large HDL percentage" ~ "Cholesterol ratio",
                           TRUE ~ title))

plot_round <- df_effects %>%
  filter(phen %in% df_r2$phen[df_r2$term == "time_day" & df_r2$pr2 > 0.01]) %>%
  ggplot(aes(x = acrophase_24hfreq, y = amplitude_24hfreq, color = type, label = title)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  #geom_col(aes(fill = type.x), alpha = 0.2) +
  geom_point(size = 0.5) +
  ggrepel::geom_label_repel(min.segment.length = 0,
    data        = labes,
    aes(label   = title, color = type),
    fill        = alpha("white", 0.5),  # only the box is 50% transparent
    size        = 3,
    label.size  = 0.2,
    show.legend = FALSE,
    max.overlaps = 20
  ) +
  coord_polar(start = 0) +
  labs(x = "Acrophase", y = "Amplitude") +
  scale_x_continuous(limits = c(0, 24), breaks = 0:23) +
  scale_color_manual(
    name   = "Omic type",
    values = c(
      "Proteomics-Olink"  = "#76B041",
      "Metabolomics-NMR"  = "#2374AB",
      "Cell_counts"       = "#8F3985",
      "Biochemistry"      = "#E85F5C"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 15,
        size  = 6
      ), nrow = 2, byrow = TRUE
    )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        text = element_text(size = 10))

ggsave("plots/plot_harmonic.png", plot_round  + theme(legend.position = "none"), width = 4, height = 4)

legend_grob <- get_legend(plot_round)


legend_plot <- as_ggplot(legend_grob)

ggsave(
  filename = "plots/legend_only.png",
  plot     = legend_plot,
  width    = 6,
  height   = 1.5
  )
