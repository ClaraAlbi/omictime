library(tidyverse)
library(paletteer)
#install.packages('ggpmisc')
library(ggpmisc)

fields <- data.table::fread("data/field.tsv")

df_effects <- readRDS("data/combined_effects.rds") %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05) #%>%
  filter(type_clean == "Proteins")

# Count how many overnight
df_effects %>%
  filter(acrophase_24hfreq < 9 | acrophase_24hfreq > 20) %>%
  count()

df_effects %>%
  filter(acrophase_24hfreq < 8 & acrophase_24hfreq > 0) %>%
  group_by(type_clean) %>%
  summarise(m = mean(acrophase_24hfreq),
            msd = sd(acrophase_24hfreq),
            n = n())

rats <- df_effects %>%
  filter(amplitude_24hfreq > 0.1) %>%
  #mutate(t = round(acrophase_24hfreq, 0)) %>%
  #group_by(t, type_clean) %>% count() %>% ungroup() %>%
  #group_by(type_clean) %>% mutate(total = sum(n)) %>%
  #ggplot(aes(x = t, y = n/total, color = type_clean)) + geom_smooth() +
  ggplot(aes(x = acrophase_24hfreq, color = type_clean)) + geom_density() +
  #facet_grid(~type_clean) +
  labs(x = "Time of day", y = "Density") +
  coord_polar() +
  scale_color_manual(
    name   = "Data type",
    values = c(
      "Proteins"  = "#76B041",
      "Metabolites"  = "#2374AB",
      "Cell counts" = "#8F3985",
      "Biochemistry" = "#E85F5C"
    )
  ) +
  scale_x_continuous(limits = c(0, 24)) +
  scale_y_continuous(limits = c(0.00, 0.2)) +
  theme_classic()

ggsave("plots/F2S_density_datatype.png", rats, width = 8, height = 6)


df_r2 <- readRDS("data/combined_variance.rds") %>%
  filter(term == "time_day") %>%
  mutate(pval = p.adjust(p.value)) %>%
  filter(pval < 0.05) #%>%
  filter(type_clean == "Proteins")

circ_diff24 <- function(a, b) {
  # returns difference in hours in [-12,12]
  delta <- (a - b + 12) %% 24 - 12
  abs(delta)
}

# 1) Read & join
df_pha <- readxl::read_xlsx("data/1-s2.0-S2352721823002401-mmc1.xlsx", skip = 1) %>%
  rename(
    acro_fund   = `acrophase of fundamental harmonic in 2-harmonic fit  (time to DLMO in hours)`,
    acro_1st    = `acrophase of first harmonic in 2-harmonic fit  (time to DLMO in hours)`,
    acro_1h = `acrophase of 1-harmonic fit (time to DLMO in hours)`
  ) %>%
  inner_join(
    df_effects %>%
      mutate(Gene = toupper(phen)) %>%
      select(Gene, acrophase_24hfreq, phen, type_clean, color_var, title, amplitude_24hfreq),
    by = "Gene"
  )
length(unique(df_pha$Gene))


p_pha <- df_pha %>%
  ggplot(aes(x = acrophase_24hfreq, y = acro_1h, color = `Rhythmic Category`)) + geom_point() +
  scale_y_continuous(
    "Acrophase SomaScan",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = paletteer_dynamic("cartography::green.pal", 2),
    guide = guide_legend(nrow = 2,
                         override.aes = list(
                           shape = 15,    # solid square
                           size  = 6      # enlarge the square
                         )
    )
  ) +
  scale_x_continuous(
    "Acrophase Olink/UKB",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

#
# df_pha %>%
#     ggplot(aes(x = acro_1st, y = acro_fund, color = `Rhythmic Category`)) + geom_point()

df_long <- df_pha %>%
  pivot_longer(
    cols      = c(acro_fund, acro_1st),
    names_to  = "harmonic",
    values_to = "acro2h"
  ) %>%
  mutate(
    # circular distance to 1‐harmonic fit
    dist24 = circ_diff24(acro2h, acrophase_24hfreq),
    # harmonic_label = case_when(harmonic == "acro_fund" ~ "2har-fundamental",
    #                            harmonic == "acro_1st" ~ "2har-1st",
    #                            harmonic == "acro_1h" ~ "1har")
  ) %>%
  group_by(`Sequence ID (Somalogic reference)`, Gene, `Rhythmic Category`, acrophase_24hfreq) %>%
  summarise(
    acro2h  = acro2h[which.min(dist24)],
    best_harmonic = harmonic[which.min(dist24)]
  )

df_long %>%
  ungroup() %>%
  count(best_harmonic)

p_best <- df_long %>%
  ggplot(aes(x = acrophase_24hfreq, y = acro2h, color = best_harmonic)) + geom_point() +
  scale_y_continuous(
    "Acrophase SomaScan",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = paletteer_dynamic("cartography::orange.pal", 2),
    guide = guide_legend(nrow = 2,
                         override.aes = list(
                           shape = 15,    # solid square
                           size  = 6      # enlarge the square
                         )
    )
  ) +
  scale_x_continuous(
    "Acrophase Olink/UKB",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

#
# p_phase <- ggplot(df_best, aes(
#   x     = acrophase_24hfreq,
#   y     = acro2h,
#   color = `Rhythmic Category`,
#   label = Gene
# )) +
#   geom_abline(linetype = "dashed", color = "grey50") +
#   geom_point(size = 2) +
#   ggrepel::geom_text_repel(data = df_best[df_best$amplitude_24hfreq > 0.3,],
#                            size = 5) +
#   scale_y_continuous(
#     "Acrophase SomaScan",
#     breaks = seq(0, 24, by = 4),
#     limits = c(0, 24),
#     expand = c(0, 0)
#   ) +
#   scale_color_manual(
#     values = paletteer_dynamic("cartography::green.pal", 2),
#     guide = guide_legend(nrow = 2,
#                          override.aes = list(
#                            shape = 15,    # solid square
#                            size  = 6      # enlarge the square
#                          )
#     )
#   ) +
#   scale_x_continuous(
#     "Acrophase Olink/UKB",
#     breaks = seq(0, 24, by = 4),
#     limits = c(0, 24),
#     expand = c(0, 0)
#   ) +
#   coord_fixed() +
#   theme_minimal(base_size = 18) +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank()
#   )


# AMPLITUDE
df_cmp <- readxl::read_xlsx("data/1-s2.0-S2352721823002401-mmc1.xlsx", skip = 1) %>%
  rename(
    amp_fund   = `amplitude of fundamental harmonic in 2-harmonic fit`,
    amp_1st    = `amplitude of first (12h) harmonic in 2-harmonic fit`,
    amp_1h = `amplitude of 1 harmonic fit`
  ) %>%
  inner_join(
    df_effects %>%
      mutate(Gene = toupper(phen)) %>%
      select(Gene, acrophase_24hfreq, phen, type_clean, title, amplitude_24hfreq),
    by = "Gene"
  ) %>%
  left_join(df_r2)

cor.test(df_cmp$amplitude_24hfreq, df_cmp$amp_1h, method = "pearson")


#
# p_amp <- ggplot(df_cmp, aes(
#   x     = amplitude_24hfreq,
#   y     = amp_1h,
#   color = `Rhythmic Category`,
#   label = Gene
# )) +
#   geom_abline(linetype = "dashed", color = "grey50") +
#   geom_point(size = 2) +
#   ggrepel::geom_text_repel(data = df_cmp[df_cmp$amplitude_24hfreq > 0.3 | df_cmp$amp_1h >0.3,],
#                            size = 5) +
#   scale_y_continuous(
#     "Amplitude SomaScan"
#   ) +
#   scale_x_continuous(
#     "Amplitude Olink/UKB"
#   ) +
#   scale_color_manual(
#     values = paletteer_dynamic("cartography::green.pal", 2),
#     guide = guide_legend(nrow = 2,
#       override.aes = list(
#         shape = 15,    # solid square
#         size  = 6      # enlarge the square
#       )
#     )
#   ) +
#   coord_fixed() +
#   theme_minimal(base_size = 18) +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank()
#   )
#
#
# p1 <- cowplot::plot_grid(p_pha, p_best, ncol = 2, labels = "AUTO")

#ggsave("plots/validation_harmonic_Specht.png", p1, width = 12, height = 6)




### Rebecca

ref <- data.table::fread("~/OneDrive - Nexus365/projects/circadian/clara_results_top19.csv") %>%
  inner_join(df_effects, by = c("Assay" = "title")) %>%
  mutate(acrophase_hr = case_when(acrophase_24hfreq > 20 & acrophase_hr < 4 ~ acrophase_hr + 24,
                                  TRUE ~ acrophase_hr))

t <- cor.test(ref$acrophase_24hfreq, ref$acrophase_hr, method = "pearson")

cor.test(ref$amplitude_24hfreq, ref$amplitude, method = "pearson")

t1 <- to_rad(ref$acrophase_24hfreq)
t2 <- to_rad(ref$acrophase_hr)

x <- cbind(t1, t2)
f <- BAMBI::circ_cor(x, type = "js")



p_olink <- ref %>%
  ggplot(aes(
    x     = acrophase_24hfreq,
    y     = acrophase_hr,
    label = Assay
  )) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 0.2) +
  stat_correlation(
    method = "pearson",
    aes(label = paste(after_stat(r.label), after_stat(p.value.label), sep = "*\", \"*")),
    parse = TRUE
  ) +
  #ggrepel::geom_text_repel(size = 2, max.overlaps = 20) +
  ggrepel::geom_label_repel(min.segment.length = 0,
                            fill        = alpha("white", 0.3),
                            size        = 3,
                            label.padding = 0.1,
                            box.padding = 0.1,
                            label.size  = 0.1,
                            show.legend = FALSE,
                            max.overlaps = 40
  ) +
  scale_y_continuous(
    "TREASURE acrophase",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 28),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    "UKB acrophase",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  ggtitle("Correlation with external 24h dataset") +
  coord_fixed() +
  theme_classic(base_size = 12) +
  theme(axis.title = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
ggsave("plots/acrophase_validation.png", p_olink, width = 4, height = 4.5)

a_olink <- ref %>%
  ggplot(aes(
    x     = amplitude_24hfreq,
    y     = amplitude,
    label = Assay
  )) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 0.2) +
  #ggrepel::geom_text_repel(size = 2, max.overlaps = 20) +
  ggrepel::geom_label_repel(min.segment.length = 0,
                            fill        = alpha("white", 0.3),
                            size        = 3,
                            label.padding = 0.1,
                            box.padding = 0.1,
                            label.size  = 0.1,
                            show.legend = FALSE,
                            max.overlaps = 40
  ) +
  scale_y_continuous(
    "TREASURE amplitude",
    #breaks = seq(0, 24, by = 4),
    limits = c(0, 0.7),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    "UKB amplitude",
    limits = c(0, 0.7),
    expand = c(0, 0)
  ) +
  coord_fixed() +
  theme_classic(base_size = 12) +
  theme(axis.title = element_text(size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank()
  )

ggsave("plots/amplitude_validation.png", a_olink, width = 4, height = 4)


### Calculate correlation
install.packages("BAMBI")

# OLINK
#
# olink_c <- df_effects %>%
#   inner_join(ref, by = c("title" ="Assay"))
#
# harmonic == "acro_fund" ~ "2har-fundamental",
# harmonic == "acro_1st" ~ "2har-1st",
# harmonic == "acro_1h" ~ "1har")
to_rad <- function(x, period = 24) (x %% period) / period * 2*pi

d_c <- df_long %>%
  filter(`Rhythmic Category` == "diurnal")

t1 <- to_rad(d_c$acrophase_24hfreq)
t2 <- to_rad(d_c$acro2h)

x <- cbind(t1, t2)
f <- BAMBI::circ_cor(x, type = "js")





# How many repeated
rep <- df_pha %>%
  group_by(`Sequence ID (Somalogic reference)`,
           Gene,
           `Rhythmic Category`) %>%
  count() %>% arrange(desc(n)) %>%
  ungroup()

gene_flags <- rep %>%
  count(Gene, `Rhythmic Category`, name = "n_rows") %>%
  pivot_wider(names_from = `Rhythmic Category`,
              values_from = n_rows,
              values_fill = 0) %>%           # make indicators
  mutate(
    circadian = as.integer(circadian > 0),
    diurnal   = as.integer(diurnal   > 0),
    both      = circadian & diurnal
  )

n_genes_both <- sum(gene_flags$both)
genes_both   <- gene_flags %>% filter(both == 1) %>% pull(Gene)
