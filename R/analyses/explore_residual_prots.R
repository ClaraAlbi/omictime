
library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(ggplot2)


### Check RES w/ protein values

prot <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_olink.rds") %>%
  select(eid, spon2, spink5, gdf15, angptl1, relt)

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2)))
df$res <- residuals(lm(pred_mean ~ time_day, data = df))
q_res <- ntile(df$res, n = 5)

df %>%
  left_join(prot) %>%
  ggplot(aes(x = time_day, y = spon2, color = factor(q_res), group = q_res)) +
  geom_smooth()


df %>%
  left_join(prot) %>%
  ggplot(aes(x = spon2, fill = factor(q_res))) +
  geom_density()
