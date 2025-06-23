library(tidyverse)

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
              mutate(title = phen)) %>%
  filter(term == "time_day")

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
      select(Gene, acrophase_24hfreq, phen, type, color_var, title, amplitude_24hfreq),
    by = "Gene"
  ) %>%
  pivot_longer(
    cols      = c(acro_fund, acro_1st, acro_1h),
    names_to  = "harmonic",
    values_to = "acro2h"
  ) %>%
  mutate(
    # circular distance to 1‐harmonic fit
    dist24 = circ_diff24(acro2h, acrophase_24hfreq),
    best_harmonic = case_when(harmonic == "acro_fund" ~ "2h-fundamental",
                              harmonic == "acro_1st" ~ "2h-1st",
                              harmonic == "acro_1h" ~ "1h")
  )

# 4) pick the single “closest” per Gene
df_best <- df_pha %>%
  group_by(Gene) %>%
  slice_min(dist24, with_ties = FALSE) %>%
  ungroup() %>%
  # build a combined legend variable
  mutate(
    combo = paste0(`Rhythmic Category`, " / ", best_harmonic)
  ) %>%
  left_join(df_r2)


p_phase <- ggplot(df_best, aes(
  x     = acrophase_24hfreq,
  y     = acro2h,
  color = combo,
  label = Gene
)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = df_best[df_best$amplitude_24hfreq > 0.2 & df_best$p.value < 0.05*3000 & df_best$pr2 > 0.01,],
                           size = 5) +
  scale_y_continuous(
    "Acrophase from 24h data (h)",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    "Acrophase from 12h data (h)",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  coord_fixed() +
  scale_color_manual(
    "Category / Closest harmonic",
    values = c(
      "circadian / 1h"   = "#D73027",
      "circadian / 2h-1st"= "#FC8D59",
      "circadian / 2h-fundamental"= "darkred",
      "diurnal / 1h"   = "#4575B4",
      "diurnal / 2h-1st"= "#91BFDB",
      "diurnal / 2h-fundamental"= "darkblue"
    ),
    guide = guide_legend(
      override.aes = list(
        shape = 15,    # solid square
        size  = 6      # enlarge the square
      )
    )
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("plots/validation_harmonic_FS2.png",p_phase, width = 10, height = 7)
