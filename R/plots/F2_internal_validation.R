
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
  filter(y > 2006)

time_i1 <- data.table::fread("/mnt/project/blood_sampling_instance1.tsv") %>%
  #filter(eid %in% preds_i0_nmr$eid) %>%
  mutate(max_time = pmax(`3166-1.0`,`3166-1.1`,`3166-1.2`,`3166-1.3`,`3166-1.4`, `3166-1.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date) %>%
  filter(y > 2006)

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
  #filter(eid %in% preds_i0_olink$eid) %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  ggtitle(label = "Initial assessment (2006-2010)", subtitle = paste0("n=", nrow(time_i0))) +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

ggsave("plots/plot_histogram_i0.png", i0_hist, width = 8, height = 8)


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
  labs(x = "Time of day") +
  ggtitle(label = "First repeat assessment (2012-13)", subtitle = paste0("n=", nrow(time_i1))) +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())


ggsave("plots/plot_histogram_i1.png", i1_hist, width = 8, height = 8)


###

out_i0 <- preds_i0_olink %>%
  group_by(cv) %>%
  nest() %>%
  mutate(N = map_dbl(data, ~nrow(.x)),
         lgb = map_dbl(data, ~cor(.x$time_day, .x$pred_lgb)^2),
         xgboost = map_dbl(data, ~cor(.x$time_day, .x$pred_xgboost)^2),
         lasso = map_dbl(data, ~cor(.x$time_day, .x$pred_lasso)^2),
         lassox2 = map_dbl(data, ~cor(.x$time_day, .x$pred_lassox2)^2)) %>%
  select(-data)

# MODELS OLINK
lgb1 <- lightgbm::lgb.load("data_share/cv.i0_lightgbm_cv1.rds")
xgb <- xgboost::xgb.load("data_share/cv.i0_xgb_cv1.rds")
lasso <- readRDS("data_share/cv.i0_lasso_cv1.rds")
lassox2 <- readRDS("data_share/cv.i0_lassox2_cv1.rds")

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
  rename(date_bsampling = date)

summary(i2_meta$y[!is.na(i2_meta$time_day)])

i2_imp <- glmnet::makeX(i2_data %>% select(-eid), na.impute = T)
i2_imp_x2 <- glmnet::makeX(i2_data %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i2 <- tibble(eid = i2_meta$eid[match(i2_data$eid, i2_meta$eid)],
                   y_test = i2_meta$time_day[match(i2_data$eid, i2_meta$eid)],
                   pred_lgb = predict(lgb1, as.matrix(i2_data %>% select(-eid))),
                   pred_xgboost = predict(xgb, as.matrix(i2_data %>% select(-eid))),
                   pred_lasso = predict(lasso, i2_imp)[,1],
                   pred_lassox2 = predict(lassox2, i2_imp_x2)[,1]) %>%
  left_join(i2_meta %>% select(eid, time_day, date_bsampling))

out_i2 <- preds_i2 %>%
  pivot_longer(-c(eid, y_test, time_day, date_bsampling)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value)^2),
         N = map_dbl(data, ~sum(!is.na(.x$value)))) %>%
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
  ggtitle(label ="Imaging (2014+)" , subtitle = paste0("n=", nrow(i2_meta))) +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
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
  rename(date_bsampling = date)

summary(i3_meta$y[!is.na(i3_meta$time_day)])

i3_imp <- glmnet::makeX(i3_data %>% select(-eid), na.impute = T)
i3_imp_x2 <- glmnet::makeX(i3_data %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i3 <- tibble(eid = i3_meta$eid[match(i3_data$eid, i3_meta$eid)],
                   y_test = i3_meta$time_day[match(i3_data$eid, i3_meta$eid)],
                   pred_lgb = predict(lgb1, as.matrix(i3_data %>% select(-eid))),
                   pred_xgboost = predict(xgb, as.matrix(i3_data %>% select(-eid))),
                   pred_lasso = predict(lasso, i3_imp)[,1],
                   pred_lassox2 = predict(lassox2, i3_imp_x2)[,1]) %>%
  left_join(i3_meta %>% select(eid, time_day, date_bsampling))

out_i3 <- preds_i3 %>%
  pivot_longer(-c(eid, y_test, time_day, date_bsampling)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value)^2),
         N = map_dbl(data, ~sum(!is.na(.x$value)))) %>%
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
  ggtitle(label = "First repeat imaging (2019+)", subtitle = paste0("n=", nrow(i3_meta))) +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

ggsave("plots/plot_histogram_i3.png", i3_hist, width = 8, height = 8)

df <- preds_i0_olink %>% mutate(i = 0) %>%
  left_join(time_i0 %>% select(eid, date_bsampling)) %>%
  select(-f) %>%
  bind_rows(preds_i2 %>% mutate(i = 2) %>% select(-y_test)) %>%
  bind_rows(preds_i3 %>% mutate(i = 3)  %>% select(-y_test)) %>%
  group_by(eid) %>% mutate(n = n())  %>% ungroup()

saveRDS(df, "olink_int_replication.rds")


# MODELS NMR


lgb1_nmr <- lightgbm::lgb.load("data_share/cv.NMR_lightgbm_cv1.rds")
xgb_nmr <- xgboost::xgb.load("data_share/cv.NMR_xgb_cv1.rds")
lasso_nmr <- readRDS("data_share/cv.NMR_lasso_cv1.rds")
lassox2_nmr <- readRDS("data_share/cv.NMR_lassox2_cv1.rds")

summary(time_i1$y[!is.na(time_i1$time_day)])

nmr_i1 <- data.table::fread("/mnt/project/nmr_nightingale_i1.tsv") %>%
  filter(!is.na(`20280-1.0`)) %>%
  mutate(across(-eid, ~scale(.x)[,1]))

i1_imp <- glmnet::makeX(nmr_i1 %>% select(-eid), na.impute = T)
i1_imp_x2 <- glmnet::makeX(nmr_i1 %>% select(-eid) %>% mutate(across(where(is.numeric), list(sq = ~ .^2), .names = "{.col}_sq")), na.impute = T)

preds_i1 <- tibble(eid = time_i1$eid[match(nmr_i1$eid, time_i1$eid)],
                   y_test = time_i1$time_day[match(nmr_i1$eid, time_i1$eid)],
                   pred_lgb = predict(lgb1_nmr, as.matrix(nmr_i1 %>% select(-eid))),
                   pred_xgboost = predict(xgb_nmr, as.matrix(nmr_i1 %>% select(-eid))),
                   pred_lasso = predict(lasso_nmr, i1_imp)[,1],
                   pred_lassox2 = predict(lassox2_nmr, i1_imp_x2)[,1]) %>%
  filter(!is.na(eid)) %>%
  left_join(time_i1 %>% select(eid, time_day, date_bsampling))

out_i1 <- preds_i1 %>%
  pivot_longer(-c(eid, y_test, time_day, date_bsampling)) %>%
  group_by(name) %>%
  nest() %>%
  mutate(r2 = map_dbl(data, ~cor(.x$y_test, .x$value, use = "complete.obs")^2),
         N = map_dbl(data, ~sum(!is.na(.x$value)))) %>%
  select(-data)


### NMR INTERNAL

df_nmr <- preds_i0_nmr %>% mutate(i = 0, time = time_day) %>%
  inner_join(time_i0 %>% select(eid, date_bsampling)) %>%
  select(-f) %>%
  bind_rows(preds_i1 %>% mutate(i = 1) %>% select(-y_test)) %>%
  group_by(eid) %>% mutate(n = n()) %>% ungroup()

saveRDS(df_nmr, "nmr_int_replication.rds")



### Histograms

library(cowplot)
plot_intval <- plot_grid(i0_hist, i1_hist, i2_hist, i3_hist, nrow = 2)

ggsave("plots/time_histograms.png", plot_intval, width = 10, height = 12)
