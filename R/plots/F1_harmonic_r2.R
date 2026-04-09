library(ggrepel)
library(patchwork)

df_r2 <- readRDS("tables/variance_covariates.rds")%>%
  mutate(pfdr = p.adjust(p.value_time_day)) %>%
  filter(pfdr < 0.05)

df_eff <- readRDS("tables/harmonic_models.rds") %>%
  mutate(pfdr = p.adjust(pvalue_h)) %>%
  filter(pfdr < 0.05)

d1 <- inner_join(df_r2, df_eff, by = c("Biomarker", "Label", "Type", "FID"))

p1 <- ggplot(d1, aes(x = amplitude_24hfreq)) +
  geom_hline(yintercept = 0.01, linetype = 2, alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.7, color = "black") +
  #geom_density(aes(color = Type)) +
  geom_point(aes(color = Type, y = t_r2_time_day), size = 1, alpha = 0.7) +
  geom_text_repel(
    data = d1 %>% filter((amplitude_24hfreq > 0.1 & t_r2_time_day > 0.01) | t_r2_time_day > 0.02 | amplitude_24hfreq > 0.37),
    size = 3,
    aes(
      x = amplitude_24hfreq,
      y = t_r2_time_day,
      label = Label,
      color = Type),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  scale_color_manual(
    values = c(
      "Proteins"  = "#76B041",
      "Metabolites"  = "#2374AB",
      "Cell counts"       = "#8F3985",
      "Biochemistry"      = "#E85F5C"
    )) +
  guides(
    color = guide_legend(override.aes = list(size = 5))
  ) +
  labs(y = "R2 time day", x = "Amplitude", color = "Data type") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray"),
        legend.position = c(0.99, 0.02),
        legend.justification = c("right", "bottom"))

p1
ggsave("plots/r2_Vs_amplitude.png", p1, width = 10, height = 10)




light_band <- data.frame(
  xmin = 6.5,
  xmax = 18.5,
  ymin = -Inf,
  ymax = Inf
)

night_band <- data.frame(
  xmin = c(0, 18.5),
  xmax = c(6.5, 24),
  ymin = -Inf,
  ymax = Inf
)

px <- d1 %>%
  filter(amplitude_24hfreq > 0.1) %>%
  ggplot(aes(x = acrophase_24hfreq, y = amplitude_24hfreq)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue", alpha = 0.2, inherit.aes = FALSE) +
  geom_vline(xintercept = 9, linetype = 2, alpha = 0.7, color = "black") +
  geom_vline(xintercept = 20, linetype = 2, alpha = 0.7, color = "black") +
  #geom_density(aes(color = Type)) +
  #geom_hline(yintercept = 0.1, linetype = 2, alpha = 0.7, color = "black") +
  #geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.7, color = "black") +
  geom_point(aes(color = Type, size = t_r2_time_day), alpha = 0.7) +
  geom_text_repel(
    data = d1 %>% filter(amplitude_24hfreq > 0.1 & t_r2_time_day > 0.01),
    size = 2.5,
    aes(label = Label),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  #coord_polar() +
  scale_x_continuous(limits = c(0, 24), breaks = 0:23,
                     expand = c(0,0)) +
  scale_color_manual(
    values = c(
      "Proteins"  = "#76B041",
      "Metabolites"  = "#2374AB",
      "Cell counts"       = "#8F3985",
      "Biochemistry"      = "#E85F5C"
    )) +
  guides(
    color = guide_legend(override.aes = list(size = 5), ncol = 1),
    size = guide_legend(ncol = 1)
  ) +
  labs(x = "Acrophase", y = "Amplitude", color = "Data type", size = "R2 time day") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "right",
        legend.title.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.direction = "vertical",
        #legend.justification = c("center", "top"),
        legend.background = element_rect(
          color = "black", fill = "white", linewidth = 0.2
        ),
        #axis.text.x = element_text(size = 14),
        #axis.title = element_text(size = 16),
        axis.line = element_blank())  +
  ggtitle("A") +
  theme(
    plot.title = element_text(face = "bold")
  )
px

ggsave("plots/plot_harmonic_all.png", px, width = 10, height = 4)

#FULL PLOT FIGURE 2



P_F <- px /
  free(p_phase + p_olink + p_pha) /
  free(p_amp + a_olink + p_sp_amp) /
  free(p_tissue) +
  plot_layout(heights = c(1, 0.8, 0.8, 1))

ggsave("plots/Figure_2.png", P_F, width = 10, height = 15, dpi = 320)
