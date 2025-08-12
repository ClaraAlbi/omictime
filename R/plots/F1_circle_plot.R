library(tidyverse)
library(cowplot)
library(ggpubr)

df_effects <- readRDS("data/combined_effects.rds") %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05)

df_r2 <- readRDS("data/combined_variance.rds") %>%
  filter(term == "time_day") %>%
  mutate(pval = p.adjust(p.value)) %>%
  filter(pval < 0.05)

df_r2 %>%
  distinct(phen, .keep_all = TRUE) %>%
  arrange(desc(pr2)) %>%
  slice_head(n = 20) %>%
  pull(title) %>% paste0(collapse = ", ")

labes <- df_effects %>%
  full_join(df_r2, by = c("phen", "color_var", "type_clean", "title")) %>%
  filter(pr2 > 0.01)

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


plot_round <- labes %>%
  ggplot(aes(x = acrophase_24hfreq, y = amplitude_24hfreq, color = type_clean, label = title)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue", alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.3, linetype = 2, color = "darkgray") +
  #geom_col(aes(fill = type.x), alpha = 0.2) +
  geom_point(size = 1.5, alpha = 0.8) +
  ggrepel::geom_label_repel(min.segment.length = 0,
    data        = subset(labes, amplitude_24hfreq > 0.3),
    aes(label   = title, color = type_clean),
    fill        = alpha("white", 0.3),  # only the box is 50% transparent
    size        = 3,
    label.size  = 0.1,
    show.legend = FALSE,
    max.overlaps = 20
  ) +
  coord_polar(start = 0) +
  labs(x = "Acrophase", y = "Amplitude") +
  scale_x_continuous(limits = c(0, 24), breaks = 0:23) +
  scale_color_manual(
    name   = "Data type",
    values = c(
      "Proteins"  = "#76B041",
      "Metabolites"  = "#2374AB",
      "Cell counts"       = "#8F3985",
      "Biochemistry"      = "#E85F5C"
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 15,
        size  = 8
      ), nrow = 1, byrow = TRUE,
    )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank(),
        text = element_text(size = 12))

ggsave("plots/plot_harmonic.png", plot_round, width = 6, height = 7)


legend_grob <- get_legend(plot_round)


legend_plot <- as_ggplot(legend_grob)

ggsave(
  filename = "plots/legend_only.png",
  plot     = legend_plot,
  width    = 6,
  height   = 1.5
  )
