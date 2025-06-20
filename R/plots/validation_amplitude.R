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


circ_diff24 <- function(a, b) {
  delta <- (a - b + 12) %% 24 - 12
  abs(delta)
}

# 1) Read & join
df_cmp <- readxl::read_xlsx("data/1-s2.0-S2352721823002401-mmc1.xlsx", skip = 1) %>%
  rename(
    amp_fund   = `amplitude of fundamental harmonic in 2-harmonic fit`,
    amp_1st    = `amplitude of first (12h) harmonic in 2-harmonic fit`,
    amp_1h = `amplitude of 1 harmonic fit`
  ) %>%
  inner_join(
    df_effects %>%
      mutate(Gene = toupper(phen)) %>%
      select(Gene, acrophase_24hfreq, phen, type, color_var, title, amplitude_24hfreq),
    by = "Gene"
  )

p_amp <- ggplot(df_cmp, aes(
  x     = amplitude_24hfreq,
  y     = amp_1h,
  color = `Rhythmic Category`,
  label = Gene
)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = df_cmp[df_cmp$amplitude_24hfreq > 0.2 | df_cmp$amp_1h > 0.2,],
                           size = 5) +
  scale_y_continuous(
    "Amplitude 24h data"
  ) +
  scale_x_continuous(
    "Amplitude 12h data"
  ) +
  scale_color_manual(values = paletteer_dynamic("cartography::green.pal", 2)) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave("plots/validation_harmonic_FS3.png", p_amp, width = 10, height = 7)
