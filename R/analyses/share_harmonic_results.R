

# Code used with the individual-level data to get beta_cos and beta_sin and intercept
m <- lm(rint_b ~ . + sex * age_recruitment, data = d %>% select(-eid, -raw, -time_day, -bmi, -smoking), na.action = "na.exclude")

# The . values are:
# time_day,
# sex,
# age_recruitment,
# fasting,
# assessment_centre,
# month_attending,
# any_of(paste0("PC", 1:20)),
# Batch,
# ppp_sel,
# bmi,
# smoking,
# all_of(p)

res_rint <- residuals(m)
m_time_only_cos <- lm(res_rint ~ cos(2 * pi * time_day / 24) + sin(2 * pi * time_day / 24), data = d)

# Transformation of regression estimates into amplitude and acrophase
df_effects <- readRDS("data/effects_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
  mutate(term = case_when(term == "(Intercept)" ~ "mesor",
                          term == "cos(2 * pi * time_day/24)" ~ "beta_cos1",
                          term == "sin(2 * pi * time_day/24)" ~ "beta_sin1")) %>%
  pivot_wider(id_cols = c(phen, type, phen), values_from = c(estimate, p.value, std.error), names_from = term) %>%
  mutate(amplitude_24hfreq = sqrt(estimate_beta_cos1^2 + estimate_beta_sin1^2),
         acrophase_24hfreq = (atan2(estimate_beta_sin1, estimate_beta_cos1) / (2 * pi) * 24 + 24) %% 24,
         q = as.integer(round(acrophase_24hfreq, 0)))

df_r2 <- readRDS("data/aov_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
  filter(term == "time_day") %>%
  select(phen, type, term, pr2, p.value)

all <- inner_join(df_effects, df_r2)
saveRDS(all, "data/olink_UKB_timeday_effects_albinana.rds")
