library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  mutate(gap = pred_lasso - time_day) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  filter(i == 0)
df$res <- residuals(lm(pred_lasso ~ time_day, data = df))


prot <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_olink.rds") %>%
  select(eid, angptl1, spon2)

df %>%
  left_join(prot) %>%
  ggplot(aes(x = spon2, y = res)) +
  geom_point()
