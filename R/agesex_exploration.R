library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
install.packages("ggrepel")
library(ggrepel)
library(ggplot2)
install.packages("broom")
install.packages("forcats")

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
  mutate(FID_clean = coalesce(LABEL, FID)) |>
  select(-LABEL)



p1 <- ggplot(d1, aes(x = amplitude_24hfreq, y = r2_time_day)) +
  geom_hline(yintercept = 0.01, linetype = 2, alpha = 0.7, color = "black") +
  geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.7, color = "black") +
  geom_point(aes(color = Type), size = 1, alpha = 0.7) +
  geom_text_repel(
    data = d1 %>% filter((amplitude_24hfreq > 0.1 & r2_time_day > 0.01) | r2_time_day > 0.02 | amplitude_24hfreq > 0.37),
    size = 2,
    aes(label = FID_clean),
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
  labs(y = "R2 time day", x = "Amplitude", color = "Data type") +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray"),
        legend.position = c(0.99, 0.02),
        legend.justification = c("right", "bottom"))

ggsave("plots/r2_Vs_amplitude.png", p1)



px <- ggplot(d1, aes(x = acrophase_24hfreq, y = amplitude_24hfreq)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue", alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 0.1, linetype = 2, alpha = 0.7, color = "black") +
  #geom_vline(xintercept = 0.1, linetype = 2, alpha = 0.7, color = "black") +
  geom_point(aes(color = Type, size = r2_time_day), alpha = 0.7) +
  geom_text_repel(
    data = d1 %>% filter((amplitude_24hfreq > 0.1 & r2_time_day > 0.01) | r2_time_day > 0.02 | amplitude_24hfreq > 0.37),
    size = 2.5,
    aes(label = FID_clean),
    max.overlaps = 30,
    box.padding = 0.2,
    point.padding = 0.3,
    segment.alpha = 0.4
  ) +
  coord_polar() +
  scale_x_continuous(limits = c(0, 24), breaks = 0:23,
                     expand = c(0,0)) +
  scale_color_manual(
    values = c(
      "Proteins"  = "#76B041",
      "Metabolites"  = "#2374AB",
      "Cell counts"       = "#8F3985",
      "Biochemistry"      = "#E85F5C"
    )) +
  labs(x = "Acrophase", y = "Amplitude", color = "Data type", size = "R2 time day") +
  theme_classic(base_size = 12) +
  theme(panel.grid.major = element_line(color = "gray"),
        legend.position = c(0.99, 0.02),
        axis.text.x = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 16),
        axis.line = element_blank(),
        legend.justification = c("right", "bottom"))

ggsave("plots/plot_circle.png", px, width = 10, height = 10)

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
