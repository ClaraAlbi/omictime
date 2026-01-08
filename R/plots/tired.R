

# 1	Not at all
# 2	Several days
# 3	More than half the days
# 4	Nearly every day

pred <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  filter(i == 0) %>%
  left_join(mh %>% select(eid, p2080_i0)) %>%
  left_join(covs)

pred$res <- residuals(lm(pred_mean ~ time_day, data = pred))

summary(lm(res ~ p2080_i0 + sex + age_recruitment, data = pred))

