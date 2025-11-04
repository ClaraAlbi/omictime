


library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
install.packages("paletteer")

covs <- data.table::fread("/mnt/project/covariates.tsv") %>%
  select(eid, sex = 2)

time <- readRDS("/mnt/project/biomarkers/time.rds")

fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

cells <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_counts.rds")

p <- cells %>%
  #slice(1:1000) %>%
  left_join(covs) %>%
  left_join(time %>% select(eid, time_day)) %>%
  pivot_longer(-c(eid, sex, time_day)) %>%
  mutate(name = as.numeric(name), time_day = round(time_day, 0)) %>%
  left_join(fields %>% select(field_id, title), by = c("name" = "field_id")) %>%
  mutate(title = if_else(
    str_length(title) > 25,
    str_c(str_sub(title, 1, 25), "\n", str_sub(title, 26)),
    title
  )) %>%
  group_by(time_day, sex, title) %>%
  summarise(vm = mean(value, na.rm = T), vsd = sd(value, na.rm = T), n = n()) %>%
  ggplot(aes(x = time_day, y = vm, color = as.factor(sex))) + geom_point() +
  geom_errorbar(aes(ymin = vm - 2*vsd, ymax = vm + 2*vsd)) +
  labs(title = "UK Biobank blood cell counts by sex", color = "Sex", y = "Mean value") +
  facet_wrap(~title, scales = "free") + theme_classic()

ggsave("blood_cell_count_sex_specific_timeday_distribution.png", p, height = 14, width = 16)



p2 <- cells %>%
  #slice(1:1000) %>%
  left_join(covs) %>%
  left_join(time %>% select(eid, time_day)) %>%
  pivot_longer(-c(eid, sex, time_day)) %>%
  mutate(name = as.numeric(name), time_day = round(time_day, 0)) %>%
  left_join(fields %>% select(field_id, title), by = c("name" = "field_id")) %>%
  mutate(title = if_else(
    str_length(title) > 25,
    str_c(str_sub(title, 1, 25), "\n", str_sub(title, 26)),
    title
  )) %>%
  ggplot(aes(x = time_day, y = value, color = as.factor(sex))) + geom_smooth() +
  labs(title = "UK Biobank blood cell counts by sex", color = "Sex", y = "Mean value") +
  facet_wrap(~title, scales = "free") + theme_classic()

ggsave("blood_cell_count_sex_specific_timeday_smooth.png", p2, height = 14, width = 16)



