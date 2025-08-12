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
              mutate(title = phen)) %>%
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
         q = as.integer(round(acrophase_24hfreq, 0)),
         p.val_joined = min(p.value_beta_sin1, p.value_beta_cos1),
         p.val_fdr = p.adjust(p.val_joined))

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
  left_join(df_r2)


p_phase <- ggplot(df_best, aes(
  x     = acrophase_24hfreq,
  y     = acro2h,
  color = combo,
  label = Gene
)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(data = df_best[df_best$amplitude_24hfreq > 0.2 & df_best$p.value < 0.05/3000 & df_best$pr2 > 0.01,],
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
      "circadian / 1har"   = "#D73027",
      "circadian / 2har-1st"= "#FC8D59",
      "circadian / 2har-fundamental"= "darkred",
      "diurnal / 1har"   = "#4575B4",
      "diurnal / 2har-1st"= "#91BFDB",
      "diurnal / 2har-fundamental"= "darkblue"
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

ggsave("plots/validation_harmonic_FS2.png", p_phase, width = 10, height = 7)



### Calculate correlation

to_rad <- function(x, period = 24) (x %% period) / period * 2*pi

circ_mean <- function(theta) atan2(mean(sin(theta)), mean(cos(theta)))

rho_js <- function(th1, th2) {
  stopifnot(length(th1) == length(th2))
  mu1 <- circ_mean(th1); mu2 <- circ_mean(th2)
  s1 <- sin(th1 - mu1); s2 <- sin(th2 - mu2)
  num <- sum(s1 * s2)
  den <- sqrt(sum(s1^2) * sum(s2^2))
  if (den == 0) return(NA_real_)
  num / den
}


# Permutation p-value for any statistic (JS or FL)
perm_pval <- function(th1, th2, stat_fun = rho_js, B = 2000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  obs <- stat_fun(th1, th2)
  cnt <- 0L
  for (b in 1:B) {
    th2p <- sample(th2, replace = FALSE)
    if (abs(stat_fun(th1, th2p)) >= abs(obs)) cnt <- cnt + 1L
  }
  list(stat = obs, p = (cnt + 1) / (B + 1))
}


d_wide <- df_best %>%
  filter(`Rhythmic Category` == "circadian")

# Example: hours → radians
t1 <- to_rad(d_wide$acrophase_24hfreq); t2 <- to_rad(d_wide$acro2h)
t1 <- to_rad(df_best$acrophase_24hfreq); t2 <- to_rad(df_best$acro2h)
perm_pval(t1, t2, stat_fun = rho_js, B = 5000)


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
