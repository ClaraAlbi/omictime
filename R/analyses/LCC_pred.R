library(lme4)
install.packages("performance")

pred <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  mutate(y = lubridate::year(date_bsampling)) %>%
    group_by(eid) %>% filter(n() == 3)
pred$res <- residuals(lm(pred_mean ~ time_day, data = pred))

m_icc <- lmer(
  res ~ 1 + (1 | eid),
  data = pred,
  REML = TRUE
)

p <- performance::icc(m_icc)

m <- lmer(
  res ~ time_day + (1 | eid),
  data = pred
)
p2 <- performance::icc(m)

m_icc_adj <- lmer(res ~ y + (1|eid), data=pred, REML=TRUE)
performance::icc(m_icc_adj)

m2 <- lmer(
  res ~ time_day + y + (1 | eid),
  data = pred
)

p2 <- performance::icc(m2)
