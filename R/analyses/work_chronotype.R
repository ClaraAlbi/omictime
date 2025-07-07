
install.packages("broom")
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

covs <- data.table::fread("covariates.tsv") %>%
  mutate(across(c(2,3, 6, 7, 9), as.factor),
         bmi = `21002-0.0` / (`50-0.0`/100)^2,
         smoking = case_when(`20116-0.0` == "-3" ~ NA_character_,
                             TRUE ~ `20116-0.0`)) %>% select(-"50-0.0", -"21002-0.0", -"20116-0.0")

colnames(covs) <- c("eid", "sex", "birth_year",  "age_recruitment",  "assessment_centre", "month_attending", "bmi", "smoking")

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20")
colnames(gen_covs) <- c("eid", paste0("PC", 1:20))

time <- readRDS("time.rds")

l <- list.files("/mnt/project/biomarkers_3", full.names = T)

data_b <- tibble(f = l[str_detect(l, "predictions")]) %>%
  mutate(d = map(f, readRDS),
         type = stringr::str_extract(f, "(?<=predictions_)([^_]+)")) %>%
  unnest(d) %>%
  pivot_wider(id_cols = c(eid, time_day), names_from = type, values_from = contains("pred")) %>%
  mutate(gap = pred_lassox2_olink - time_day) %>%
  select(eid, time_day, gap, pred_lassox2_olink) %>%
  left_join(covs) %>%
  left_join(gen_covs)

sleep <- data.table::fread("/mnt/project/chronotype2.tsv")

data_b %>%
  left_join(sleep %>% select(eid, h_sleep = `1160-0.0`,
                             chrono = `1180-0.0`,
                             ever_insomnia = `1200-0.0`,
                             wakeup = `1170-0.0`)) %>%
  #mutate(wakeup = 5 - wakeup) %>%
  filter(ever_insomnia %in% 1:3) %>%
  filter(chrono %in% 1:4) %>%
  filter(wakeup %in% 1:4) %>%
  filter(h_sleep > 2 & h_sleep < 18) %>%
  #filter(time_day > 11 & time_day < 18) %>%
  mutate(dys = abs(gap) > 3) %>%
  ggplot(aes(x = time_day, y = gap, color = as.factor(chrono))) +
  geom_smooth(linewidth = 2) +
  labs(y = "Acceleration", color = "Chronotype \nMorning to Evening", x = "Recorded time") +
  scale_color_viridis_d(direction = -1) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.position = "bottom"
  )


job_vars <- data.table::fread("/mnt/project/job_vars.tsv") %>%
  mutate(night_shift = case_when(`3426-0.0` == 1 ~ "Never",
                                 `3426-0.0` == 2 ~ "Sometimes",
                                 `3426-0.0` == 3 ~ "Usually",
                                 `3426-0.0` == 4 ~ "Always"),
         sleep_h = `826-0.0`) %>%
  filter(`3426-0.0` %in% 1:4) %>%
  mutate(night_shift = as.factor(night_shift),
         night_shift = relevel(night_shift, ref = "Never"))


nu <- data_b %>%
  rename(pred_time = pred_lassox2_olink) %>%
  filter(!is.na(gap)) %>%
  left_join(job_vars) %>%
  group_by(night_shift) %>%
  summarise(n = n(),
            mean_gap = mean(abs(gap),na.rm = T),
            mean_h = mean(time_day, na.rm = T),
            mean_p = mean(pred_time, na.rm = T)) %>%
  mutate(night_shift = factor(night_shift, levels = c("Never", "Sometimes", "Usually", "Always"))) %>%
  filter(!is.na(night_shift))


NIGHT_DP<- bind_rows(broom::tidy(lm(abs(gap) ~ night_shift, data= data_b %>%
                                      left_join(job_vars)), conf.int = TRUE) %>% mutate(type = "Dysregulation"),
                     broom::tidy(lm(gap ~ night_shift, data = data_b %>%
                                      left_join(job_vars)), conf.int = TRUE) %>% mutate(type = "Acceleration")) %>%
  group_by(type) %>%
  mutate(m = case_when(term == "(Intercept)" ~ "yes",
                       term != "Intercept" ~ "no"),

         intercept = estimate[m == "yes"][1],  # Pull the intercept for each group
         predicted_mean = if_else(
           m == "yes",
           intercept,                         # Intercept stays as is
           intercept + estimate               # Add intercept to other terms
         ),
         conf_high_p = if_else(
           m == "yes",
           conf.high,                         # Intercept stays as is
           intercept + conf.high               # Add intercept to other terms
         ),
         conf_low_p = if_else(
           m == "yes",
           conf.low,                         # Intercept stays as is
           intercept + conf.low               # Add intercept to other terms
         ),
         night_shift = case_when(term == "(Intercept)" ~ "Never",
                                 term == "night_shiftAlways" ~ "Always",
                                 term == "night_shiftSometimes" ~ "Sometimes",
                                 term == "night_shiftUsually" ~ "Usually"),
         night_shift = factor(night_shift, levels = c("Never", "Sometimes", "Usually", "Always"))) %>%
  left_join(nu) %>%
  mutate(night_shift_label = paste0(night_shift, "\nn=", n),
         night_shift = factor(night_shift_label,
                              levels = paste0(c("Never", "Sometimes", "Usually", "Always"),
                                              "\nn=", nu$n[match(c("Never", "Sometimes", "Usually", "Always"), nu$night_shift)]), ordered = T)
  )

pwork <- NIGHT_DP %>% ggplot(aes(x = night_shift, y = predicted_mean, fill = night_shift)) +
  geom_col() +
  geom_errorbar(aes(ymin = conf_low_p, ymax = conf_high_p), width = 0.4) +
  scale_fill_viridis_d(option = "mako", direction = -1, begin = 0.2) +
  facet_wrap(~type, scales = "free") +
  labs(x = "Night Shift Work") +
  scale_x_discrete(labels = levels(NIGHT_DP$night_shift)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.title.y = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 0.7, size = 12)
  )

ggsave(filename = "work.png", plot = pwork, width = 6, height = 7)

