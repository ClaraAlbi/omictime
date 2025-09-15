library(tidyverse)

fields <- data.table::fread("data/field.tsv")

df_effects <- readRDS("data/combined_effects.rds") %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05)

df_r2 <- readRDS("data/combined_variance.rds") %>%
  filter(term == "time_day") %>%
  mutate(pval = p.adjust(p.value)) %>%
  filter(pval < 0.05)

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
  ) %>%
  pivot_longer(
    cols      = c(acro_fund, acro_1st, acro_1h),
    names_to  = "harmonic",
    values_to = "acro2h"
  ) %>%
  mutate(
    # circular distance to 1‐harmonic fit
    dist24 = circ_diff24(acro2h, acrophase_24hfreq),
    best_harmonic = case_when(harmonic == "acro_fund" ~ "2har-fundamental",
                              harmonic == "acro_1st" ~ "2har-1st",
                              harmonic == "acro_1h" ~ "1har")
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
  left_join(df_r2) %>%
  left_join(df_effects)


p_phase <- ggplot(df_best, aes(
  x     = acrophase_24hfreq,
  y     = acro2h,
  color = `Rhythmic Category`,
  label = Gene
)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = df_best[df_best$amplitude_24hfreq > 0.3,],
                           size = 5) +
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

p_amp <- ggplot(df_cmp, aes(
  x     = amplitude_24hfreq,
  y     = amp_1h,
  color = `Rhythmic Category`,
  label = Gene
)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = df_cmp[df_cmp$amplitude_24hfreq > 0.3,],
                           size = 5) +
  scale_y_continuous(
    "Amplitude SomaScan"
  ) +
  scale_x_continuous(
    "Amplitude Olink/UKB"
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
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


p1 <- cowplot::plot_grid(p_phase, p_amp, ncol = 2, labels = "AUTO")


ggsave("plots/validation_harmonic_Specht.png", p1, width = 12, height = 6)




### Rebecca

ref <- data.table::fread("~/OneDrive - Nexus365/projects/circadian/clara_results_top19.csv")


p_olink <- df_best %>%
  inner_join(ref, by = c("Gene" ="Assay")) %>%
  ggplot(aes(
    x     = acrophase_24hfreq,
    y     = acrophase_hr,
    label = Gene
  )) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 5) +
  scale_y_continuous(
    "Acrophase from >24h data. Olink",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    "Acrophase from 12h data. Olink",
    breaks = seq(0, 24, by = 4),
    limits = c(0, 24),
    expand = c(0, 0)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )


a_olink <- df_best %>%
  inner_join(ref, by = c("Gene" ="Assay")) %>%
  ggplot(aes(
    x     = amplitude_24hfreq,
    y     = amplitude,
    label = Gene
  )) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(size = 5) +
  scale_y_continuous(
    "Amplitude from >24h data. Olink", limits = c(0, 0.6)
  ) +
  scale_x_continuous(
    "Amplitude from 12h data. Olink", limits = c(0, 0.6)
  ) +
  coord_fixed() +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p2 <- cowplot::plot_grid(p_olink, a_olink, ncol = 2, labels = "AUTO")

ggsave("plots/validation_harmonic_olink.png", p2, width = 12, height = 6)


### Calculate correlation
to_rad <- function(x, period = 24) (x %% period) / period * 2*pi

install.packages("BAMBI")

t1 <- to_rad(d_wide$acrophase_24hfreq)
t2 <- to_rad(d_wide$acro2h)

x <- cbind(t1, t2)
f <- BAMBI::circ_cor(x, type = "js")

olink_c <- df_best %>%
  inner_join(ref, by = c("Gene" ="Assay"))

t1 <- to_rad(olink_c$acrophase_24hfreq)
t2 <- to_rad(olink_c$acrophase_hr)

x <- cbind(t1, t2)
BAMBI::circ_cor(x, type = "js")


# Permutation p-value for any statistic (JS or FL)
perm_pval <- function(th1, th2, stat_fun = rho_js, B = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  obs <- stat_fun(cbind(th1, th2))
  cnt <- 0L
  for (b in 1:B) {
    th2p <- sample(th2, replace = FALSE)
    if (abs(stat_fun(cbind(th1, th2p))) >= abs(obs)) cnt <- cnt + 1L
  }
  list(stat = obs, p = (cnt + 1) / (B + 1))
}


d_wide <- df_best %>%
  filter(`Rhythmic Category` == "circadian")

# Example: hours → radians
perm_pval(t1, t2, stat_fun = BAMBI::circ_cor, B = 5000)

cor.test(t2, t1)

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
