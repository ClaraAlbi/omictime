library(dplyr)
library(tidyr)
install.packages("broom")
library(broom)
library(stringr)

sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid,
         h_sleep = `1160-0.0`,
         chrono = `1180-0.0`,
         ever_insomnia = `1200-0.0`,
         wakeup = `1170-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ "Definitely morning",
    chrono == 2 ~ "Rather morning",
    chrono == -1~ "Neither",
    chrono == 3 ~ "Rather evening",
    chrono == 4 ~ "Definitely evening",
    TRUE ~ NA_character_)) %>% rename(chronotype_assess = chrono) %>%
  mutate(chronotype_assess = factor(chronotype_assess, levels = c("Definitely morning", "Rather morning", "Neither", "Rather evening", "Definitely evening")),
         chronotype_assess = relevel(chronotype_assess, ref = "Definitely morning"))

sleep_new <- data.table::fread("/mnt/project/sleep_online.csv") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  mutate(
    score_wake_up = case_when(
      p30425 == "5:00am - 6:30am" ~ 5,
      p30425 == "6:30am - 7:45am" ~ 4,
      p30425 == "7:45am - 9:45am" ~ 3,
      p30425 == "9:45am - 11:00am" ~ 2,
      p30425 == "11:00am - 12 noon" ~ 1,
      TRUE ~ NA_real_
    ),

    score_tired_after_waking = case_when(
      p30426 == "Very tired" ~ 1,
      p30426 == "Fairly tired" ~ 2,
      p30426 == "Fairly refreshed" ~ 3,
      p30426 == "Very refreshed" ~ 4,
      TRUE ~ NA_real_
    ),

    score_evening_tiredness = case_when(
      p30427 == "8:00pm - 9:00pm" ~ 5,
      p30427 == "9:00pm - 10:15pm" ~ 4,
      p30427 == "10:15pm - 12:45am" ~ 3,
      p30427 == "12:45am - 2:00am" ~ 2,
      p30427 == "2:00am - 3:00am" ~ 1,
      TRUE ~ NA_real_
    ),

    score_best_time = case_when(
      p30428 == "5:00am - 8:00am" ~ 5,
      p30428 == "8:00am - 10:00am" ~ 4,
      p30428 == "10:00am - 5:00pm" ~ 3,
      p30428 == "5:00pm - 10:00pm" ~ 2,
      p30428 == "10:00pm - 5:00am" ~ 1,
      TRUE ~ NA_real_
    ),

    score_chronotype = case_when(
      p30429 == "Definitely a morning-type" ~ 6,
      p30429 == "Rather more a morning-type than an evening-type" ~ 4,
      p30429 == "Rather more an evening-type than a morning-type" ~ 2,
      p30429 == "Definitely an evening-type" ~ 0,
      TRUE ~ NA_real_
    ),
    rmeq_score = if_else(
      !is.na(score_wake_up) &
        !is.na(score_tired_after_waking) &
        !is.na(score_evening_tiredness) &
        !is.na(score_best_time) &
        !is.na(score_chronotype),

      score_wake_up +
        score_tired_after_waking +
        score_evening_tiredness +
        score_best_time +
        score_chronotype,

      NA_real_
    ),
    rmeq_chronotype = case_when(
      rmeq_score >= 22 ~ "Definitely morning",
      rmeq_score >= 18 & rmeq_score <= 21 ~ "Rather morning",
      rmeq_score >= 12 & rmeq_score <= 17 ~ "Neither",
      rmeq_score >= 8  & rmeq_score <= 11 ~ "Rather evening",
      rmeq_score <= 7  ~ "Definitely evening",
      TRUE ~ NA_character_
    ),
    rmeq_chronotype = factor(rmeq_chronotype, levels = c("Definitely morning", "Rather morning", "Neither", "Rather evening", "Definitely evening")),
    rmeq_chronotype = relevel(rmeq_chronotype, ref = "Definitely morning"),
    p30429 = factor(p30429, levels = c("Definitely a morning-type", "Rather more a morning-type than an evening-type","Do not know", "Rather more an evening-type than a morning-type", "Definitely an evening-type"), labels = c("Definitely morning", "Rather morning", "Neither", "Rather evening", "Definitely evening")),
    p30429 = relevel(p30429, ref = "Definitely morning")
  )



sleep_new %>%
  mutate(y = lubridate::year(p30489)) %>%
  count(y)

df <- readRDS("/mnt/project/olink_int_panels14.rds") %>%
  filter(!is.na(time_day)) %>%
  #filter(i == 0) %>%
  left_join(sleep_new) %>%
  left_join(covs) %>%
  left_join(a) %>%
  left_join(gen_covs) %>%
  left_join(dep)
df$res <- residuals(lm(pred_mean ~ time_day, data = df))


df %>%
  ggplot(aes(x = res, y = rmeq_chronotype)) +
  geom_boxplot()

vars <- c("time_day", "age_recruitment", "sex", "chrono", "h_sleep", "ever_insomnia", "p30079", #"rmeq_chronotype", "rmeq_score",
          "is_weekend", "day_of_week",
          "season", "night_shift", "smoking", "bmi", "is_dst", "wakeup", "shift_work", "employment", "TDI", "autumnDST", "springDST")

covars <- c("sex", "age_recruitment", "assessment_centre", paste0("PC", 1:20))

v <-

f_CNT <- as.formula(paste("res ~", paste(c("p30429", covars), collapse = " + ")))
f_rMEQ <- as.formula(paste("res ~", paste(c("rmeq_chronotype", covars), collapse = " + ")))


fit1 <- lm(f_CNT, data = df)
fit2 <- lm(f_rMEQ, data = df)

data <- bind_rows(tidy(fit1, conf.int = T) %>% mutate(mod = "chrono") %>% filter(str_detect(term, "p30429|Intercept")) %>%
  mutate(term = str_remove(term, "p30429"),
         across(c(estimate, conf.low, conf.high), ~case_when(term == "(Intercept)" ~ 0, TRUE ~ .x)),
         p.value  = case_when(term == "(Intercept)" ~ NA_integer_, TRUE ~ p.value),
         term = case_when(term == "(Intercept)" ~ "Definitely morning",TRUE ~ term)),
  tidy(fit2, conf.int = T) %>% mutate(mod = "rmeq") %>% filter(str_detect(term, "rmeq_chronotype|Intercept")) %>%
    mutate(term = str_remove(term, "rmeq_chronotype"),
           across(c(estimate, conf.low, conf.high), ~case_when(term == "(Intercept)" ~ 0, TRUE ~ .x)),
           p.value  = case_when(term == "(Intercept)" ~ NA_integer_, TRUE ~ p.value),
           term = case_when(term == "(Intercept)" ~ "Definitely morning",
                            TRUE ~ term))) %>%
  mutate(FDR = p.adjust(p.value),
         in_mins = estimate * 60)

# term               estimate std.error statistic   p.value conf.low conf.high mod          FDR in_mins
# <chr>                 <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl> <chr>      <dbl>   <dbl>
#   1 Definitely morning    0        0.265       2.81 NA           0        0      chrono NA           0
# 2 Rather morning       -0.148    0.0226     -6.57  5.06e-11   -0.193   -0.104  chrono  1.01e-10   -8.90
# 3 Neither              -0.257    0.0342     -7.52  5.55e-14   -0.324   -0.190  chrono  1.67e-13  -15.4
# 4 Rather evening       -0.384    0.0264    -14.6   8.04e-48   -0.436   -0.332  chrono  6.43e-47  -23.0
# 5 Definitely evening   -0.499    0.0355    -14.1   1.01e-44   -0.568   -0.429  chrono  7.05e-44  -29.9
# 6 Definitely morning    0        0.305       2.86 NA           0        0      rmeq   NA           0
# 7 Rather morning       -0.177    0.0470     -3.76  1.69e- 4   -0.269   -0.0846 rmeq    1.69e- 4  -10.6
# 8 Neither              -0.410    0.0460     -8.93  4.90e-19   -0.501   -0.320  rmeq    2.45e-18  -24.6
# 9 Rather evening       -0.642    0.0574    -11.2   7.32e-29   -0.754   -0.529  rmeq    4.39e-28  -38.5
# 10 Definitely evening   -1.35     0.164      -8.22  2.15e-16   -1.67    -1.03   rmeq    8.59e-16  -81.1

p<- data %>%
  mutate(term = factor(term, levels = c("Definitely morning", "Rather morning", "Neither", "Rather evening", "Definitely evening")),
         term = relevel(term, ref = "Definitely morning")) %>%
  ggplot(aes(x = estimate, y = fct_rev(term), shape = FDR < 0.05, color = mod)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3, fill = "white", position = position_dodge(width = 0.5)) +
  scale_shape_manual(
    values = c(`TRUE` = 15, `FALSE` = 22),
    na.value = 20, guide = "none"
  ) +
  scale_color_manual(values = c("rmeq" = "orange", "chrono" = "lightblue"), labels = c("rmeq" = "rMEQ", "chrono" = "Single-item chronotype")) +
  labs(x = "CA (β, 95% CI)", y = "Chronotype", shape = NULL,  color = "Chronotype definition") +
  theme_classic(base_size = 10) +
  guides(alpha = "none") +
  theme(
    panel.grid.major = element_line(color = "gray"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    plot.title.position = "plot"
  )

ggsave("plots/FX_rmeq_chrono.png", p, width = 6, height = 5)


"One hears about 'morning-types' and 'evening-types.' Which one of these types do you consider yourself to be?"

#### CHECK IF PEOPLE THAT CHANGE GROUP HAD HIGHER DYSREGULATION

#Fluvial plot
install.packages("ggalluvial")
library(ggalluvial)

alluvial_data <- sleep %>% inner_join(sleep_new) %>%
  filter(!is.na(chronotype_assess) & !is.na(rmeq_chronotype)) %>%
  count(chronotype_assess, rmeq_chronotype)


ggplot(alluvial_data,
       aes(axis1 = chronotype_assess,
           axis2 = rmeq_chronotype,
           y = n)) +
  geom_alluvium(aes(fill = chronotype_assess), width = 1/12) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Assessment Centre",  "rMEQ-Based"),
                   expand = c(.05, .05)) +
  labs(
    y = "Number of Participants"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")



### rMEQ compared to blood based

df %>%
  filter(!is.na(rmeq_chronotype)) %>%
  ggplot(aes(x = factor(round(time_day, 0)), y = res, fill = rmeq_chronotype)) +
  geom_boxplot() +
  theme_minimal()

chronotype_avgs <- df %>%
  left_join(sleep) %>%
  group_by(chronotype_assess) %>%
  summarise(avg_res = mean(res, na.rm = TRUE))

df2 <- df %>%
  left_join(chronotype_avgs, by = "chronotype_assess") %>%
  mutate(diff_from_type = res - avg_res) %>%

df %>%
  ggplot(aes(x = res, y = rmeq_score)) +
  geom_smooth()

summary(lm(res ~ rmeq_chronotype, data = df))
summary(lm(res ~ chronotype_assess, data = df))
summary(lm(res ~ rmeq_score, data = df))
