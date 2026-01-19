library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

time <- readRDS("/mnt/project/biomarkers/time.rds")

mh <- data.table::fread("/mnt/project/mh_time.csv")

chrono <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid, chrono = `1180-0.0`, qchrono = `1180-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ "Morning",
    chrono == 2 ~ "Morning",
    chrono == 3 ~ "Evening",
    chrono == 4 ~ "Evening",
    TRUE ~ NA_character_))

pred <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  filter(i == 0) %>%
  filter(!is.na(time_day)) %>%
  left_join(mh) %>%
  left_join(chrono) %>%
  left_join(covs) %>%
  mutate(
    p29017 = factor(
      p29017,
      levels = c(
        "Mood did not vary",
        "In the evening or at night",
        "In the morning",
        "Do not know"
      )
    )
  ) %>%
  filter(!is.na(p29017)) %>%
  filter(!is.na(chrono))

a <- pred %>% count(chrono, p29017)

pred$res <- residuals(lm(pred_mean ~ time_day, data = pred))

summary(lm(res~chrono, data = pred))

fit <- lm(res ~ p29017 + sex*age_recruitment, data = pred %>% filter(chrono == "Evening"))
fit2 <- lm(res ~ p29017 + sex*age_recruitment, data = pred %>% filter(chrono == "Morning"))
summary(fit)
summary(fit2)

pred %>% ggplot(aes(x = p29017, y = res, color = chrono)) +
  geom_boxplot()

ggplot(pred, aes(x = chrono, fill = p29017)) +
  geom_bar(position = "dodge")

ggplot(pred, aes(x = chrono, y = res, fill = p29017)) +
  geom_boxplot()


summary(lm(res~p29017, data = pred %>% filter(chrono == 'Definitely morning')))

summary(lm(res~p29017, data = pred %>% filter(chrono == 'Definitely evening')))

summary(lm(res~p29017, data = pred %>% filter(chrono == 'Rather morning')))

summary(lm(res~p29017, data = pred %>% filter(chrono == 'Rather evening')))

summary(lm(res~p29017, data = pred))
