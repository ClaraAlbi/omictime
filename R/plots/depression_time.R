library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

time <- readRDS("/mnt/project/biomarkers/time.rds")

mh <- data.table::fread("mh_time.csv")


pred <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  filter(i == 0) %>%
  filter(!is.na(time_day)) %>%
  left_join(mh) %>%
  left_join(sleep) %>%
  mutate(
    p29017 = factor(
      p29017,
      levels = c(
        "In the morning",
        "In the evening or at night",
        "Mood did not vary",
        "Do not know"
      )
    )
  ) %>%
  filter(!is.na(p29017)) %>%
  filter*(!is.na(chrono))

a <- pred %>% count(chrono, p29017)

pred$res <- residuals(lm(pred_mean ~ time_day, data = pred))

summary(lm(res~p29017*chrono, data = pred))

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
