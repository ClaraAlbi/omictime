library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
install.packages("lightgbm")
library(lightgbm)
install.packages("xgboost")
install.packages("glmnet")


l <- list.files("/mnt/project/biomarkers_3", full.names = T)

preds_i0_olink <- tibble(f = l[str_detect(l, "predictions_i0")]) %>%
  mutate(d = map(f, readRDS)) %>%
  unnest(d)

preds_i0_nmr <- tibble(f = l[str_detect(l, "predictions_NMR")]) %>%
  mutate(d = map(f, readRDS)) %>%
  unnest(d)

### HISTOGRAMS

time_i0 <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  filter(eid %in% preds_i0_olink$eid)

time_i1 <- data.table::fread("blood_sampling_instance1.tsv") %>%
  filter(eid %in% preds_i0_nmr$eid) %>%
  mutate(max_time = pmax(`3166-1.0`,`3166-1.1`,`3166-1.2`,`3166-1.3`,`3166-1.4`, `3166-1.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  filter(time_day >= 9 & time_day <= 20)

light_band <- data.frame(
  xmin = 6,
  xmax = 20,
  ymin = -Inf,
  ymax = Inf
)

night_band <- data.frame(
  xmin = c(0, 20),
  xmax = c(6, 24),
  ymin = -Inf,
  ymax = Inf
)


i0_hist <- time_i0 %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

ggsave("plots/plot_histogram_i0_olink.png", i0_hist, width = 8, height = 8)


i1_hist <- time_i1 %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())


ggsave("plots/plot_histogram_i1_olink.png", i1_hist, width = 8, height = 8)


###

out_i0 <- preds_i0 %>%
  group_by(cv) %>%
  nest() %>%
  mutate(N = map_dbl(data, ~nrow(.x)),
         lgb = map_dbl(data, ~cor(.x$time_day, .x$pred_lgb)^2),
         xgboost = map_dbl(data, ~cor(.x$time_day, .x$pred_xgboost)^2),
         lasso = map_dbl(data, ~cor(.x$time_day, .x$pred_lasso)^2),
         lassox2 = map_dbl(data, ~cor(.x$time_day, .x$pred_lassox2)^2)) %>%
  select(-data)

# MODELS OLINK
lgb1 <- readRDS("/mnt/project/biomarkers_3/cv.i0_lightgbm_cv1.rds")
xgb <- readRDS("/mnt/project/biomarkers_3/cv.i0_xgb_cv1.rds")
lasso <- readRDS("/mnt/project/biomarkers_3/cv.i0_lasso_cv1.rds")
lassox2 <- readRDS("/mnt/project/biomarkers_3/cv.i0_lassox2_cv1.rds")

lightgbm::lgb.save(lgb1, "cv.i0_lightgbm_cv1.rds")
xgboost::xgb.save(xgb, "cv.i0_xgb_cv1.rds")


### validation

i2_data <- data.table::fread("/mnt/project/OLINK_i2.tsv") %>%
  mutate(across(-eid, ~scale(.x)[,1]))
i2_meta <- data.table::fread("/mnt/project/blood_sampling_instance2.tsv") %>%
  filter(!is.na(`3166-2.0`)) %>%
  filter(eid %in% i2_data$eid) %>%
  mutate(max_time = pmax(`3166-2.0`,`3166-2.1`,`3166-2.2`,`3166-2.3`,`3166-2.4`, `3166-2.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  filter(time_day >= 9 & time_day <= 20)

summary(i2_meta$y[!is.na(i2_meta$time_day)])

i2_imp <- glmnet::makeX(i2_data %>% select(-eid), na.impute = T)
i2_imp_x2 <- glmnet::makeX(i2_data %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i2 <- tibble(eid = i2_meta$eid[match(i2_data$eid, i2_meta$eid)],
                   y_test = i2_meta$time_day[match(i2_data$eid, i2_meta$eid)],
                 pred_lgb = predict(lgb.load("cv.i0_lightgbm_cv1.rds"), as.matrix(i2_data %>% select(-eid))),
                 pred_xgb = predict(xgboost::xgb.load("cv.i0_xgb_cv1.rds"), as.matrix(i2_data %>% select(-eid))),
                 pred_lasso = predict(lasso, i2_imp)[,1],
                 pred_lassox2 = predict(lassox2, i2_imp_x2)[,1])

out_i2 <- preds_i2 %>%
  pivot_longer(-c(y_test, eid)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value)^2)) %>%
  select(-data)

i2_hist <- i2_meta %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

ggsave("plots/plot_histogram_i2.png", i2_hist, width = 8, height = 8)



# i3
i3_data <- data.table::fread("/mnt/project/OLINK_i3.tsv") %>%
  mutate(across(-eid, ~scale(.x)[,1]))
i3_meta <- data.table::fread("/mnt/project/blood_sampling_instance3.tsv") %>%
  filter(!is.na(`3166-3.0`)) %>%
  filter(eid %in% i3_data$eid) %>%
  mutate(max_time = pmax(`3166-3.0`,`3166-3.1`,`3166-3.2`,`3166-3.3`,`3166-3.4`, `3166-3.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  filter(time_day >= 9 & time_day <= 20)

summary(i3_meta$y[!is.na(i3_meta$time_day)])

i3_imp <- glmnet::makeX(i3_data %>% select(-eid), na.impute = T)
i3_imp_x2 <- glmnet::makeX(i3_data %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i3 <- tibble(eid = i3_meta$eid[match(i3_data$eid, i3_meta$eid)],
                   y_test = i3_meta$time_day[match(i3_data$eid, i3_meta$eid)],
                pred_lgb = predict(lgb.load("cv.i0_lightgbm_cv1.rds"), as.matrix(i3_data %>% select(-eid))),
                pred_xgb = predict(xgboost::xgb.load("cv.i0_xgb_cv1.rds"), as.matrix(i3_data %>% select(-eid))),
                pred_lasso = predict(lasso, i3_imp)[,1],
                pred_lassox2 = predict(lassox2, i3_imp_x2)[,1])

out_i3 <- preds_i3 %>%
  pivot_longer(-c(y_test, eid)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value)^2)) %>%
  select(-data)

i3_hist <- i3_meta %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

ggsave("plots/plot_histogram_i3.png", i3_hist, width = 8, height = 8)

## AGExSEX plots
covs <- data.table::fread("/mnt/project/covariates.tsv") %>%
  select(eid, `31-0.0`, `34-0.0`, `21022-0.0`)
colnames(covs) <- c("eid", "Sex", "year_birth", "Age_baseline")

#data_i0 <-
preds_i0 %>%
  left_join(covs) %>%
  mutate(Age = Age_baseline, Sex = as.factor(Sex),
         gap = time_day - pred_lasso,
         absgap = abs(gap)) %>%
  pivot_longer(c(gap, absgap)) %>%
  #group_by(Age, Sex) %>% summarise(m_gap = mean(gap), n = n())
  ggplot(aes(x = Age, y = value, color = Sex)) +
  geom_smooth() +
  facet_wrap(~name, scales = "free")

covs_i2 <- preds_i2 %>%
  left_join(i2_meta) %>%
  left_join(covs) %>%
  mutate(Age = y - year_birth, Sex = as.factor(Sex),
         gap = y_test - pred_lasso,
         absgap = abs(gap))
covs_i2 %>%
  pivot_longer(c(gap, absgap)) %>%
  #group_by(Age, Sex) %>% summarise(m_gap = mean(gap), n = n())
  ggplot(aes(x = Age, y = value, color = Sex)) +
  geom_smooth() +
  facet_wrap(~name, scales = "free")

covs_i3 <- preds_i3 %>%
  left_join(i3_meta) %>%
  left_join(covs) %>%
  mutate(Age = y - year_birth, Sex = as.factor(Sex),
         gap = y_test - pred_lasso,
         absgap = abs(gap))

covs_i3 %>%
  pivot_longer(c(gap, absgap)) %>%
  #group_by(Age, Sex) %>% summarise(m_gap = mean(gap), n = n())
  ggplot(aes(x = Age, y = value, color = Sex)) +
  geom_smooth() +
  facet_wrap(~name, scales = "free")


i2i3 <- inner_join(covs_i3, covs_i2, by = c("eid", "Sex"))

i2i3 %>% ggplot(aes(x = pred_lasso.x, y = pred_lasso.y, color = Sex)) +
  geom_point()

i2i3 %>% ggplot(aes(x = y_test.x, y = y_test.y, color = Sex)) +
  geom_point()


# MODELS NMR
lgb1_nmr <- readRDS("/mnt/project/biomarkers_3/cv.NMR_lightgbm_cv1.rds")
xgb_nmr <- readRDS("/mnt/project/biomarkers_3/cv.NMR_xgb_cv1.rds")
lasso_nmr <- readRDS("/mnt/project/biomarkers_3/cv.NMR_lasso_cv1.rds")
lassox2_nmr <- readRDS("/mnt/project/biomarkers_3/cv.NMR_lassox2_cv1.rds")

lightgbm::lgb.save(lgb1_nmr, "cv.NMR_lightgbm_cv1.rds")
xgboost::xgb.save(xgb_nmr, "cv.NMR_xgb_cv1.rds")

summary(time_i1$y[!is.na(time_i1$time_day)])

nmr_i1 <- data.table::fread("nmr_nightingale_i1.tsv") %>%
  filter(!is.na(`20280-1.0`)) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

i1_imp <- glmnet::makeX(nmr_i1 %>% select(-eid), na.impute = T)
i1_imp_x2 <- glmnet::makeX(nmr_i1 %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i1 <- tibble(eid = time_i1$eid[match(nmr_i1$eid, time_i1$eid)],
                   y_test = time_i1$time_day[match(nmr_i1$eid, time_i1$eid)],
                   pred_lgb = predict(lgb.load("cv.NMR_lightgbm_cv1.rds"), as.matrix(nmr_i1 %>% select(-eid))),
                   pred_xgb = predict(xgboost::xgb.load("cv.NMR_xgb_cv1.rds"), as.matrix(nmr_i1 %>% select(-eid))),
                   pred_lasso = predict(lasso_nmr, i1_imp)[,1],
                   pred_lassox2 = predict(lassox2_nmr, i1_imp_x2)[,1])

out_i1 <- preds_i1 %>%
  pivot_longer(-c(y_test, eid)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value, use = "complete.obs")^2)) %>%
  select(-data)
