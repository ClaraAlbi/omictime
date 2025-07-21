
install.packages("broom")
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)



sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid, h_sleep = `1160-0.0`,
         chrono = `1180-0.0`,
         ever_insomnia = `1200-0.0`,
         wakeup = `1170-0.0`) %>%
  filter(chrono %in% 1:4)

df %>% left_join(sleep) %>%
  ggplot(aes(x = age_recruitment, y = res, color = factor(chrono))) + geom_smooth() +
  facet_grid(~sex)

df %>% ggplot(aes(x = age_recruitment, y = res, color = factor(sex))) + geom_smooth()
