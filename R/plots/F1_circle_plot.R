library(tidyverse)
library(cowplot)
library(ggpubr)

df_effects <- readRDS("data/combined_effects.rds") %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05)

df_effects %>% group_by(type_clean) %>% count()
df_effects %>% group_by(type_clean) %>% summarise(mean_amplitude = mean(amplitude_24hfreq),
                                                  sd_amplitude = sd(amplitude_24hfreq),
                                                  n = n())

df_r2 <- readRDS("data/combined_variance.rds") %>%
  filter(term == "time_day") %>%
  mutate(pval = p.adjust(p.value)) %>%
  filter(pval < 0.05)

df_r2 %>%
  distinct(phen, .keep_all = TRUE) %>%
  arrange(desc(t_r2)) %>%
  slice_head(n = 20) %>%
  pull(title) %>% paste0(collapse = ", ")

labes <- df_effects %>%
  inner_join(df_r2, by = c("phen", "color_var", "type_clean", "title"))

labes %>% filter(amplitude_24hfreq > 0.3) %>% count()


data.table::fwrite(labes, "data/resulting_phases.txt")

light_band <- data.frame(
  xmin = 5.4,
  xmax = 20.5,
  ymin = -Inf,
  ymax = Inf
)

night_band <- data.frame(
  xmin = c(0, 20.5),
  xmax = c(5.4, 24),
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
  geom_point(size = 0.2, alpha = 0.2) +
  ggrepel::geom_label_repel(min.segment.length = 0,
    data        = subset(labes, amplitude_24hfreq > 0.3 | t_r2 > 0.03),
    aes(label   = title, color = type_clean),
    fill        = alpha("white", 0.7),
    size        = 3,
    label.padding = 0.1,
    box.padding = 0.25,
    label.size  = 0.1,
    show.legend = FALSE,
    max.overlaps = 50
  ) +
  coord_polar(start = 0, clip = "off") +
  labs(x = "Acrophase", y = "Amplitude") +
  scale_x_continuous(limits = c(0, 24), breaks = 0:23,
                     expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(
    name   = "Data type",
    values = c(
      "Proteins"  = "#76B041",
      "Metabolites"  = "#2374AB",
      "Cell counts" = "#8F3985",
      "Biochemistry" = "#E85F5C"
    )
  ) +
  theme_classic() +
  theme(text = element_text(size = 12),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = margin(t = 0, 0, 0, 0), vjust = 5),
        legend.position = "none",
        legend.title = element_blank(),
        plot.margin = margin(0,0,0,0))

ggsave("plots/plot_harmonic.png", plot_round, width = 6, height = 7)

pf <- cowplot::plot_grid(plot_round, p_olink, ncol = 2,
                   rel_widths = c(1, 1), rel_heights = c(1, 1),
                   align = "hv", labels = c("A", "B"))


