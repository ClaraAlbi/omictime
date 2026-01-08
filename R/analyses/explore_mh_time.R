
library(ggplot2)
library(tidyr)
library(lubridate)
library(dplyr)

mh <- data.table::fread("/mnt/project/psychosocial_MH.csv")

time <- readRDS("/mnt/project/biomarkers/time.rds")

data <- mh %>%
  mutate(across(-eid, as.character)) %>%
  left_join(covs %>% select(eid, sex, age_recruitment)) %>%
  mutate(age_b = cut(
    age_recruitment,
    breaks = c(40, 50, 60, 70),
    right = FALSE,
    labels = c("40–49", "50–59", "60–69")
  )) %>% select(-age_recruitment) %>%
  pivot_longer(
    cols = -c(eid, sex, age_b),
    names_to = "name",
    values_to = "value"
  ) %>%
  filter(!is.na(value), value != "") %>%
  left_join(time, by = "eid") %>%
  filter(time_day > 9, time_day < 21) %>%
  mutate(t = factor(round(time_day))) %>%
  group_by(t, name, value, sex, age_b) %>%
  summarise(total = n(), .groups = "drop") %>%
  group_by(t, name, sex, age_b) %>%
  mutate(freq = total / sum(total)) %>%
  ungroup()


library(ggplot2)


data %>%
  filter(name == "p2080_i0") %>%
  filter(!is.na(age_b)) %>%
  ggplot(aes(x = t, y = freq, color = value, group = value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ sex + age_b, scales = "free_y") +
  labs(
    x = "Time (rounded hour)",
    y = "Frequency",
    color = "Value"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = "none"
  )


neu <- mh %>%
  left_join(covs) %>%
  mutate(age_b = cut(
    age_recruitment,
    breaks = c(40, 50, 60, 70),
    right = FALSE,
    labels = c("40–49", "50–59", "60–69")
  )) %>%
  filter(!is.na(age_b)) %>%
  #select(eid, p20127_i0) %>%
  left_join(time, by = "eid") %>%
  filter(time_day > 9, time_day < 21) %>%
  mutate(t = factor(round(time_day)))


neu %>% ggplot(aes(x = t, y = p20127_i0, fill = age_b)) +
  geom_boxplot() +
  facet_grid(~sex)

summary(lm(p20127_i0 ~ t, data = neu %>% filter(sex == "Female")))

