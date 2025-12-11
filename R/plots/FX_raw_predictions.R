library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)

l <- c(list.files("/mnt/project/circadian/results/models",
                  pattern = "predictions", full.names = TRUE))

# l <- c(list.files("/mnt/project/biomarkers_3",
#                   pattern = "predictions", full.names = TRUE)[-c(31:35, 1:5, 16:20)],
#        list.files("/mnt/project/biomarkers_3/covariate_res/MODELS",
#                   pattern = "predictions", full.names = TRUE))

preds_i0_olink <- tibble(f = l[str_detect(l, "tech_14")]) %>%
  mutate(d = map(f, readRDS)) %>%
  unnest(d) %>%
  rowwise() %>% mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  unnest()



data <- preds_i0_olink %>%
  pivot_longer(contains("pred")) %>%
  group_by(name, cv) %>%
  nest() %>%
  mutate(
    N    = map_dbl(data, ~ sum(!is.na(.x[[4]]))),
    r2 = map_dbl(data, ~ cor(.x$time_day, .x$value)^2)) %>%
  select(-data) %>%
  group_by(name) %>% summarise(m_r2 = mean(r2))

preds_i0_olink$res <- residuals(lm(pred_mean ~ time_day, data = preds_i0_olink))

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  filter(smoking != "-3") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
  )
gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20", `22006-0.0`, `22021-0.0`, `22000-0.0`)
colnames(gen_covs) <- c("eid", paste0("PC", 1:20), "is_white", "rel", "batch")


df <- preds_i0_olink %>%
  left_join(covs) %>%
  left_join(data.table::fread("/mnt/project/ancestry_new.csv")) %>%
  left_join(gen_covs %>% select(eid, contains("PC"))) %>%
  filter(age_recruitment > 39) %>%
  mutate(age_g = case_when(age_recruitment <= 50 ~ "40-50",
                           age_recruitment <= 60 & age_recruitment > 50 ~ "50-60",
                           age_recruitment <= 70 & age_recruitment > 60 ~ "60-70"))

ggplot(df, aes(x = age_recruitment, y = res, color = sex)) + geom_smooth() +
  theme_minimal() +
  facet_grid(~assessment_centre)

summary(lm(res~ age_recruitment*sex + as.factor(assessment_centre) + p30079, data =df))

summary(lm(res~ 0 + as.factor(p30079) + PC1 + PC2 + PC3, data =df))
