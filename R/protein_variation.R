library(tidyr)
library(dplyr)
library(ggplot2)

olink <- data.table::fread("/mnt/project/olink_instance_0.csv")
#time <- readRDS("/mnt/project/biomarkers/time.rds")

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  ungroup()

df$res <- residuals(lm(pred_mean ~ time_day, data = df))
df$res_q <- factor(ntile(df$res, 5))

olink$tnr

d <- olink %>%
  left_join(df)

d %>%
  ggplot(aes(x = angptl1, y = res, color = factor(ntile(res, 5)))) + geom_point()


d %>%
  ggplot(aes(x = time_day, y = tnr, color = factor(ntile(res, 5)))) + geom_smooth()

d %>%
  pivot_longer(c(angptl1, tnr)) %>%
  ggplot(aes(x = time_day, y = value, color = name)) + geom_point()

d %>%
  pivot_longer(c(angptl1, tnr)) %>%
  ggplot(aes(x = time_day, y = value, color = res_q) + geom_point() +
           facet_grid(~name) +
           theme_minimal()


         d %>%
           mutate(r = as.factor(round(res, 0))) %>%
           #filter(r %in% c(-6, -5, -4, 4,5,6)) %>%
           ggplot(aes(x = angptl1, y = tnr, color = as.factor(round(res, 0)))) + geom_point()


         d %>%
           filter(eid %in% c(5807521, 4082402)) %>%
           select(eid, time_day, res, res_q, spon2, angptl1, hyal1, tnr, spink5, relt, gdf15)
