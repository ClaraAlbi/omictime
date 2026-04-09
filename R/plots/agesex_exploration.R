library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
install.packages("ggplot2")
library(ggplot2)
install.packages("broom")
install.packages("forcats")
install.packages("patchwork")
library(patchwork)
install.packages("ggrepel")
library(ggrepel)
install.packages("ggpmisc")
library(ggpmisc)



df_r2 <- readRDS("tables/variance_covariates.rds")%>%
  mutate(pfdr = p.adjust(p.value_time_day)) %>%
  filter(pfdr < 0.05)

df_eff <- readRDS("tables/harmonic_models.rds") %>%
  mutate(pfdr = p.adjust(pvalue_h)) %>%
  filter(pfdr < 0.05)

d1 <- inner_join(df_r2, df_eff, by = c("Biomarker", "Label", "Type", "FID"))

set_prots <- d1 %>%
  filter(r2_time_day > 0.01) %>%
  filter(amplitude_24hfreq > 0.1)


olink <- readRDS("/mnt/project/biomarkers_res/res_olink_tech.rds") %>%
  select(eid, any_of(set_prots$FID))

nmr <- readRDS("/mnt/project/biomarkers_res/res_nmr_tech.rds") %>%
  select(eid, any_of(set_prots$FID))

bio <- readRDS("/mnt/project/biomarkers_res/res_labs_tech.rds") %>%
  select(eid, any_of(set_prots$FID))

cc <- readRDS("/mnt/project/biomarkers_res/res_counts_tech.rds") %>%
  select(eid, any_of(set_prots$FID))

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


biomarkers <- colnames(d)[2:132]

base_formula <- function(outcome) {
  as.formula(paste0(
    outcome,
    " ~ cos_time + sin_time + age_recruitment + sex + ",
    "bmi + smoking + fasting + month_attending"
  ))
}

interaction_vars <- c("sex")

# results <- map(
#   biomarkers,
#   function(biomarker) {
#     map(
#       interaction_vars,
#       function(int_var) {
#         run_time_interaction(
#           data = d,
#           outcome = biomarker,
#           interaction_var = int_var)})})
#

#saveRDS(results, "data_share/associations_sexbytime.rds")

library(tidyverse)

results <- readRDS("data_share/associations_sexbytime.rds")

anova_table <- map_dfr(
  results,
  ~ map_dfr(.x, "anova")
) %>% mutate(pfdr = p.adjust(p_value)) %>%
  filter(pfdr < 0.05)

estimate_table <- map_dfr(
  results,
  ~ map_dfr(.x, "estimates")
) %>% filter(outcome %in% anova_table$outcome)

compute_amp_phase <- function(beta_cos, beta_sin) {
  amp <- sqrt(beta_cos^2 + beta_sin^2)
  phase_rad <- atan2(beta_sin, beta_cos)
  phase_hr <- (24 / (2 * pi)) * phase_rad
  phase_hr <- (phase_hr + 24) %% 24  # wrap to [0,24)

  tibble(
    amplitude = amp,
    acrophase = phase_hr
  )
}

sex_amp_phase <- estimate_table %>%
  filter(model == "time_x_sex") %>%
  filter(term %in% c(
    "cos_time", "sin_time",
    "cos_time:sexMale", "sin_time:sexMale"
  )) %>%
  select(outcome, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  mutate(
    # Female (reference)
    female = pmap(
      list(cos_time, sin_time),
      compute_amp_phase
    ),
    # Male (reference + interaction)
    male = pmap(
      list(cos_time + `cos_time:sexMale`,
           sin_time + `sin_time:sexMale`),
      compute_amp_phase
    )
  ) %>%
  transmute(
    outcome,
    female = map(female, ~ mutate(.x, sex = "Female")),
    male   = map(male,   ~ mutate(.x, sex = "Male"))
  ) %>%
  pivot_longer(c(female, male), values_to = "res") %>%
  unnest(res)

sex_amp_phase_wide <- sex_amp_phase %>%
  left_join(d1, by = c("outcome" = "FID")) %>%
  select(Label, sex, amplitude, acrophase) %>%
  pivot_wider(
    id_cols   = Label,
    names_from  = sex,
    values_from = c(amplitude, acrophase)
  ) %>% mutate(dif_amp = amplitude_Female - amplitude_Male,
               dif_phase = acrophase_Female - acrophase_Male)

t2_amp <- cor.test(sex_amp_phase_wide$amplitude_Female, sex_amp_phase_wide$amplitude_Male)
# data:  sex_amp_phase_wide$amplitude_Female and sex_amp_phase_wide$amplitude_Male
# t = 10.597, df = 70, p-value = 3.398e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.6761208 0.8601282
# sample estimates:
#   cor
# 0.7848532

t2 <- cor.test(sex_amp_phase_wide$acrophase_Female, sex_amp_phase_wide$acrophase_Male)
# Pearson's product-moment correlation
# data:  sex_amp_phase_wide$acrophase_Female and sex_amp_phase_wide$acrophase_Male
# t = 35.276, df = 70, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9570797 0.9830760
# sample estimates:
#       cor
# 0.9730078

p_amp <- ggplot(sex_amp_phase_wide,
       aes(x = amplitude_Female, y = amplitude_Male)) +
  geom_point() +
  geom_label_repel(
    data = sex_amp_phase_wide %>% filter(abs(dif_amp) > 0.3),
    size = 2.5,
    aes(label = Label),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  ggpmisc::stat_poly_eq(
    mapping    = aes(label = paste("italic(R) ==", round(t2_amp$estimate, 2))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("italic(p) ==", formatC(t2_amp$p.value, format = "e", digits = 0))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.90,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("n == 72")),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 3.5,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "Female amplitude", y = "Male amplitude") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray"))+
  ggtitle("E") +
  theme(
    plot.title = element_text(face = "bold")
  )

p_phase <- ggplot(sex_amp_phase_wide,
       aes(x = acrophase_Female, y = acrophase_Male)) +
  geom_point() +
  geom_label_repel(
    data = sex_amp_phase_wide %>% filter(abs(dif_phase) > 4),
    size = 2.5,
    aes(label = Label),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  ggpmisc::stat_poly_eq(
    mapping    = aes(label = paste("italic(R) ==", round(t2$estimate, 2))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.95,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("italic(p) ==", formatC(t2$p.value, format = "e", digits = 0))),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.90,
    size = 3.5,
    color = "black"
  ) +
  ggpmisc::stat_poly_eq(
    mapping = aes(label = paste("n == 72")),
    parse = TRUE,
    label.x = 0.05,
    label.y = 0.85,
    size = 3.5,
    color = "black"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  labs(x = "Female acrophase (h)", y = "Male acrophase (h)") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray")) +
  ggtitle("B") +
  theme(
    plot.title = element_text(face = "bold")
  )








####

ts1 <- d %>%
  select(eid, any_of(biomarkers), time_day) %>%
  pivot_longer(-c(eid, time_day)) %>%
  mutate(t = round(time_day, 0)) %>%
  group_by(name, t) %>% summarise(n        = n(),
                                  mean_val = mean(value, na.rm = T),
                                  sd_val   = sd(value, na.rm = T))


fid_levels <- ts1 %>%
  left_join(d1, by = c("name" = "FID")) %>%
  distinct(FID_clean, r2_time_day) %>%
  ungroup() %>%
  arrange(desc(r2_time_day)) %>%
  pull(FID_clean)

p_raw <- ts1 %>%
  left_join(d1, by = c("name" = "FID")) %>%
  mutate(
    FID_clean = factor(FID_clean, levels = fid_levels),
    se_val   = sd_val / sqrt(n),
    ci_lower = mean_val - 1.96 * se_val,
    ci_upper = mean_val + 1.96 * se_val
  ) %>%
  ggplot(aes(x = t, y = mean_val)) +
  geom_point(size = 0.2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  facet_wrap(~FID_clean, ncol = 10) +
  labs(y = "Biomarker residuals", x = "Time of day") +
  theme_minimal(base_size = 8)

ggsave("plots/biom_time_int.png", p_raw, width = 10, height = 14)




ts <- d %>%
  select(eid, any_of(biomarkers), time_day, sex) %>%
  pivot_longer(-c(eid, time_day, sex)) %>%
  mutate(t = round(time_day, 0)) %>%
  filter(name %in% anova_table$outcome) %>%
  group_by(name, sex, t)  %>% summarise(n        = n(),
                                        mean_val = mean(value, na.rm = T),
                                        sd_val   = sd(value, na.rm = T))



fid_levels_sex <- ts %>%
  ungroup() %>%
  left_join(d1, by = c("name" = "FID")) %>%
  distinct(FID_clean, r2_time_day) %>%
  arrange(desc(r2_time_day)) %>%
  pull(FID_clean)

p2 <- ts %>%
  ungroup() %>%
  left_join(d1, by = c("name" = "FID")) %>%
  mutate(
    FID_clean = factor(FID_clean, levels = fid_levels_sex),
    se_val   = sd_val / sqrt(n),
    ci_lower = mean_val - 1.96 * se_val,
    ci_upper = mean_val + 1.96 * se_val
  ) %>%
  ggplot(aes(x = t, y = mean_val, color = sex)) +
  geom_point(size = 0.2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  facet_wrap(~FID_clean, ncol = 10) +
  labs(y = "Biomarker residuals", x = "Time of day") +
  theme_minimal(base_size = 8)

ggsave("plots/sex_time_int.png", p2, width = 10, height = 12)

ts %>%
  pivot_wider(id_cols = c(name, t), names_from = sex, values_from = mean_val) %>%
  mutate(f_top = Female > Male) %>%
  group_by(name) %>%
  summarise(f_t = sum(f_top)) %>%
  arrange(desc(f_t)) %>% count(f_t)
  hist(.)
