
install.packages("broom")
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

covs <- readRDS("/mnt/project/biomarkers/covs.rds")

oth <- fread("covariates.csv")

sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid, h_sleep = `1160-0.0`,
         chrono = `1180-0.0`,
         ever_insomnia = `1200-0.0`,
         wakeup = `1170-0.0`) %>%
  filter(chrono %in% 1:4)


df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  mutate(gap = pred_lasso - time_day) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  filter(i == 0)

data <- df %>%
  inner_join(covs) %>%
  inner_join(sleep) %>%
  mutate(Sex = factor(sex, labels = c("Female", "Male")))


data %>%
  ggplot(aes(x = age_recruitment, y = gap, color = factor(chrono))) + geom_smooth() + facet_grid(~sex)


data %>%
  ggplot(aes(x = as.factor(month_attending), y = gap)) + geom_violin()

summary(lm(abs(gap) ~ as.factor(month_attending), data = data))
