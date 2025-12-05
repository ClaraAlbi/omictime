
library(ggplot2)


olink_cohort <- data.table::fread("/mnt/project/olink_instance_0.csv")

sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid,
         chrono = `1180-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ 2,
    chrono == 2 ~ 1,
    chrono == -1~ 0,
    chrono == 3 ~ -1,
    chrono == 4 ~ -2))

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  select(eid, time_day, pred_mean)
df$res <- residuals(lm(pred_mean ~ time_day, data = df))

d <- olink_cohort %>%
  left_join(sleep) %>%
  left_join(df)

summary(lm(chrono ~ res, data = d))
summary(lm(chrono ~ spon2, data = d))
summary(lm(chrono ~ relt, data = d))
summary(lm(chrono ~ gdf15, data = d))


  ggplot(aes(x = chrono, y = spon2)) + geom_point()

