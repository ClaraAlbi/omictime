library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
install.packages("broom")
install.packages("ggplot2")

labels <- data.table::fread("data_share/supplementary_data1_explanatory_labels.csv") |>
  mutate(
    FID   = str_squish(FID),
    LABEL = str_squish(Label)
  ) %>%
  mutate(LABEL = case_when(Type == "Proteins"~ Name,
         TRUE ~ LABEL))

labels <- labels |>
  mutate(
    LABEL = if_else(
      LABEL == tolower(LABEL) & str_detect(LABEL, "^[a-z0-9]+$"),
      toupper(LABEL),
      LABEL
    )
  )

rs <- data.table::fread("data_share/supplementary_data1.csv") %>%
  mutate(pfdr = p.adjust(p.value_time_day)) %>%
  filter(pfdr < 0.05) #%>%
  filter(r2_time_day > 0.01)

amp <- data.table::fread("data_share/supplementary_data2.csv") %>%
  mutate(pfdr = p.adjust(pvalue_h)) %>%
  filter(pfdr < 0.05) #%>%
  filter(amplitude_24hfreq > 0.1)

d1 <- amp %>%
  inner_join(rs, by = c("FID", "Name", "Type"))

d1 <- d1 |>
  left_join(labels |> select(FID, LABEL), by = "FID") |>
  mutate(FID = coalesce(LABEL, FID)) |>
  select(-LABEL)



p1 <- ggplot(d1, aes(x = amplitude_24hfreq, y = r2_time_day)) +
  geom_hline(yintercept = 0.01, linetype = 2, alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.7, color = "black") +
  geom_point(aes(color = Type), size = 1, alpha = 0.7) +
  geom_text_repel(
    data = d1 %>% filter((amplitude_24hfreq > 0.1 & r2_time_day > 0.01) | r2_time_day > 0.03 | amplitude_24hfreq > 0.37),
    size = 2,
    aes(label = FID),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  labs(y = "R2 time day", x = "Amplitude", color = "") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(),
        legend.position = c(0.99, 0.02),
        legend.justification = c("right", "bottom"))

ggsave("see.png", p1)

set_prots <- intersect(rs$FID, amp$FID)set_prots <- intersect(rs$FID, amp$FID)set_prots <- intersect(rs$FID, amp$FID)

olink <- readRDS("/mnt/project/biomarkers_res/res_olink_tech.rds") %>%
  select(eid, any_of(rs$FID))

nmr <- readRDS("/mnt/project/biomarkers_res/res_nmr_tech.rds") %>%
  select(eid, any_of(rs$FID))

bio <- readRDS("/mnt/project/biomarkers_res/res_labs_tech.rds") %>%
  select(eid, any_of(rs$FID))

cc <- readRDS("/mnt/project/biomarkers_res/res_counts_tech.rds") %>%
  select(eid, any_of(rs$FID))

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight / (height/100)^2,
         age_group = case_when(
           age_recruitment >= 40 & age_recruitment < 50 ~ "40-50",
           age_recruitment >= 50 & age_recruitment < 60 ~ "50-60",
           age_recruitment >= 60 & age_recruitment < 70 ~ "60-70",
      TRUE ~ NA_character_
    ),
    age_group = factor(age_group, levels = c("40-50", "50-60", "60-70")),
    sex = factor(sex, labels = c("Female", "Male")),
    smoking = factor(smoking),
    month_attending = factor(month_attending)) %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds")) %>%
  mutate(cos_time = cos(2 * pi * time_day / 24),
         sin_time = sin(2 * pi * time_day / 24))


d <- bio %>%
  left_join(nmr) %>%
  left_join(cc) %>%
  left_join(olink) %>%
  left_join(covs)


###

run_time_interaction <- function(data, outcome, interaction_var, alpha = 0.05) {

  f_base <- paste0("`",
                   outcome, "`",
                   " ~ cos_time + sin_time + age_recruitment + sex + bmi + smoking + fasting + month_attending"
  )
  print(outcome)

  f_int <- paste0("`",
    outcome, "`",
    " ~ (cos_time + sin_time) * ", interaction_var,
    " + age_recruitment + sex + bmi + smoking + fasting + month_attending"
  )

  m0 <- lm(f_base, data = data)
  m1 <- lm(f_int, data = data)

  anova_res <- anova(m0, m1) %>%
    as.data.frame() %>%
    slice(2) %>%
    transmute(
      outcome = outcome,
      interaction = interaction_var,
      df = Df,
      f_stat = F,
      p_value = `Pr(>F)`
    )

  use_interaction <- anova_res$p_value < alpha

  est_res <- if (use_interaction) {
    broom::tidy(m1) %>%
      mutate(
        outcome = outcome,
        model = paste0("time_x_", interaction_var)
      )
  } else {
    broom::tidy(m0) %>%
      mutate(
        outcome = outcome,
        model = "baseline"
      )
  }

  list(
    anova = anova_res,
    estimates = est_res,
    interaction_kept = use_interaction
  )
}


biomarkers <- colnames(d)[2:135]

base_formula <- function(outcome) {
  as.formula(paste0(
    outcome,
    " ~ cos_time + sin_time + age_recruitment + sex + ",
    "bmi + smoking + fasting + month_attending"
  ))
}

interaction_vars <- c("sex")

results <- map(
  biomarkers,
  function(biomarker) {
    map(
      interaction_vars,
      function(int_var) {
        run_time_interaction(
          data = d,
          outcome = biomarker,
          interaction_var = int_var)})})

anova_table <- map_dfr(
  results,
  ~ map_dfr(.x, "anova")
) %>% mutate(pfdr = p.adjust(p_value)) %>%
  filter(pfdr < 0.05)

estimate_table <- map_dfr(
  results,
  ~ map_dfr(.x, "estimates")
) %>% filter(outcome %in% anova_table$outcome)




ts <- olink %>% left_join(covs %>% select(eid, sex, age_recruitment)) %>%
  pivot_longer(-c(eid, sex, age_recruitment)) %>%
  filter(name %in% anova_table$outcome) %>%
  filter(age_recruitment > 39) %>%
  group_by(name, sex, age_recruitment) %>% summarise(m_d = mean(value, na.rm = T), sd_d  = sd(value, na.rm = T))

ggplot(ts, aes(x = age_recruitment, y = m_d, color = sex)) +
  geom_point() +
  facet_wrap(~name, scales = "free")



ts <- d %>%
  select(any_of(biomarkers), age_recruitment, sex, time_day) %>%
  pivot_longer(-c(age_recruitment, sex, time_day)) %>%
  filter(name %in% anova_table$outcome) %>%
  filter(age_recruitment > 39) %>%
  mutate(t = round(time_day, 0)) %>%
  group_by(name, sex, t) %>% summarise(m_d = mean(value, na.rm = T), sd_d  = sd(value, na.rm = T))


ts2 <- ts %>% ungroup() %>%
  mutate(name = recode(name, !!!setNames(fields$title, fields$field_id)))


name_lookup <- tibble::tribble(
  ~original, ~short,
  "Spectrometer-corrected alanine", "Ala_corr",
  "Alanine", "Alanine",
  "Isoleucine", "Isoleucine",
  "Leucine", "Leucine",
  "Valine", "Valine",
  "Phenylalanine", "Phenylala",
  "Total Concentration of Branched-Chain Amino Acids (Leucine + Isoleucine + Valine)", "BCAA_total",
  "Glucose", "Glucose",
  "Lactate", "Lactate",
  "3-Hydroxybutyrate", "Hydroxybut",
  "Acetate", "Acetate",
  "Acetoacetate", "Acetoacet",
  "Acetone", "Acetone",
  "Phosphate", "Phosphate",
  "Total bilirubin", "Bilirubin",
  "Triglycerides", "Triglyc",
  "Monounsaturated Fatty Acids", "MUFA",
  "Monounsaturated Fatty Acids to Total Fatty Acids percentage", "MUFA_pct",
  "Phospholipids in HDL", "HDL_PL",
  "Average Diameter for HDL Particles", "HDL_Diam",
  "Concentration of Very Large HDL Particles", "HDL_VL_cnt",
  "Total Lipids in Very Large HDL", "HDL_VL_TL",
  "Phospholipids in Very Large HDL", "HDL_VL_PL",
  "Free Cholesterol in Very Large HDL", "HDL_VL_FC",
  "Triglycerides in Very Large HDL", "HDL_VL_TG",
  "Total Lipids in Large HDL", "HDL_L_TL",
  "Phospholipids in Large HDL", "HDL_L_PL",
  "Free Cholesterol in Large HDL", "HDL_L_FC",
  "Triglycerides in Large HDL", "HDL_L_TG",
  "Triglycerides in Medium HDL", "HDL_M_TG",
  "Free Cholesterol in Small LDL", "LDL_S_FC",
  "Free Cholesterol to Total Lipids in Very Small VLDL percentage", "VLDL_VS_FCp",
  "Free Cholesterol to Total Lipids in Small LDL percentage", "LDL_S_FCp",
  "Phospholipids to Total Lipids in Very Large HDL percentage", "HDL_VL_PLp",
  "Cholesterol to Total Lipids in Very Large HDL percentage", "HDL_VL_Chp",
  "Cholesteryl Esters to Total Lipids in Very Large HDL percentage", "HDL_VL_CEp",
  "Free Cholesterol to Total Lipids in Very Large HDL percentage", "HDL_VL_FCp",
  "Free Cholesterol to Total Lipids in Large HDL percentage", "HDL_L_FCp",
  "Phospholipids to Total Lipids in Small HDL percentage", "HDL_S_PLp",
  "White blood cell (leukocyte) count", "Leukocytes",
  "Lymphocyte count", "Lymphocytes",
  "Monocyte count", "Monocytes",
  "Neutrophill count", "Neutrophils",
  "Basophill count", "Basophills",
  "actn2", "ACTN2",
  "adamts15", "ADAMTS15",
  "adm", "ADM",
  "agr3", "AGR3",
  "angptl1", "ANGPTL1",
  "angptl4", "ANGPTL4",
  "c1qtnf5", "C1QTNF5",
  "ccn3", "CCN3",
  "fas", "FAS",
  "fst", "FST",
  "gip", "GIP",
  "glb1", "GLB1",
  "hspb6", "HSPB6",
  "il6", "IL6",
  "inhbb", "INHBB",
  "lgals1", "LGALS1",
  "mep1a", "MEP1A",
  "muc13", "MUC13",
  "mybpc1", "MYBPC1",
  "myl3", "MYL3",
  "pgf", "PGF",
  "pla2g10", "PLA2G10",
  "pomc", "POMC",
  "pth", "PTH",
  "sdc1", "SDC1",
  "smoc1", "SMOC1",
  "spon2", "SPON2",
  "timp4", "TIMP4",
  "tmprss15", "TMPRSS15",
  "tnr", "TNR"
)


ts2 <- ts2 |>
  left_join(name_lookup, by = c("name" = "original")) |>
  mutate(name = dplyr::coalesce(short, name)) |>
  select(-short)

p2 <- ggplot(ts2, aes(x = t, y = m_d, color = sex)) +
  geom_point() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Biomarker residuals", x = "Time of day") +
  theme_minimal(base_size = 8) +
  theme(strip.text = element_text(size = 6),
        axis.title = element_text(size = 10),
        legend.position = "bottom")

ggsave("plots/sex_time_int.png", p2, width = 10, height = 10)

ts2 %>%
  pivot_wider(id_cols = c(name, t), names_from = sex, values_from = m_d) %>%
  mutate(f_top = Female > Male) %>%
  group_by(name) %>%
  summarise(f_t = sum(f_top)) %>%
  arrange(desc(f_t)) %>% pull(f_t) %>%
  hist(.)
