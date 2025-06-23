library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lightgbm)
library(xgboost)
library(glmnet)

lgb1 <- lightgbm::lgb.load("data_share/cv.olink_lightgbm_cv1.rds")
xgb <- xgboost::xgb.load("data_share/cv.olink_xgb_cv1.rds")
lasso <- readRDS("data_share/cv.olink_lasso_cv1.rds")
lassox2 <- readRDS("data_share/cv.olink_lassox2_cv1.rds")

# File with time of blood sampling info
time_file <- data.table::fread("/mnt/project/blood_sampling.tsv")

# Raw proteomics file
olink_raw_file <- data.table::fread("/mnt/project/olink_instance_0.csv")
# Order of proteins in model
prot_order <- lasso$glmnet.fit$beta@Dimnames[[1]]

# Order dataset and scale all raw values
ckb_olink <- olink_raw_file %>%
  select(-glipr1) %>% # Remove this protein because it had too few values in UKB
  select(eid, any_of(tolower(prot_order))) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

# The UKB has 6 available times, each for a blood sample replicate. I just pick the latest one because the empty value is "1900-01-01"
time_day <- time_file %>%
  mutate(max_time = pmax(`3166-0.0`,`3166-0.1`,`3166-0.2`,`3166-0.3`,`3166-0.4`, `3166-0.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  filter(time_day >= 9 & time_day <= 20)

# Impute values needed to use glmnet
ckb_olink_imp <- glmnet::makeX(ckb_olink %>% select(-eid), na.impute = T)
ckb_olink_imp_x2 <- glmnet::makeX(ckb_olink %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

# Individual-lebel predictions
out_ckb <- tibble(eid = time_day$eid[match(ckb_olink$eid, time_day$eid)],
                      y_test = time_day$time_day[match(ckb_olink$eid, time_day$eid)],
                      pred_lgb = predict(lgb1, as.matrix(ckb_olink %>% select(-eid))),
                      pred_xgb = predict(xgb, as.matrix(ckb_olink %>% select(-eid))),
                      pred_lasso = predict(lasso, ckb_olink_imp)[,1],
                      pred_lassox2 = predict(lassox2, ckb_olink_imp_x2)[,1])

# Prediction accuracy estimation
pred <- out_ckb %>%
  pivot_longer(c(-y_test, -eid)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value, use = "complete.obs")^2)) %>%
  select(-data)

saveRDS(pred, "prediction_olink_ckb.rds")


### Validate chronotype associations if possible

sleep <- data.table::fread("/mnt/project/chronotype2.tsv")

plot <- out_ckb %>%
  mutate(gap = y_test - pred_lasso) %>%
  left_join(sleep %>% select(eid, chrono = `1180-0.0`)) %>%
  filter(chrono %in% 1:4) %>%
  ggplot(aes(x = y_test, y = gap, color = as.factor(chrono))) +
  geom_smooth(linewidth = 2) +
  labs(y = "Acceleration", color = "Chronotype \nMorning to Evening", x = "Recorded time") +
  scale_color_viridis_d(direction = -1) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        legend.position = "bottom"
  )

ggsave("chronotype_ckb_proteotime.png", plot, width = 10, height = 10)
