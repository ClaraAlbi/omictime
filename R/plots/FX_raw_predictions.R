

l <- c(list.files("/mnt/project/circadian/results/models/",
                  pattern = "predictions", full.names = TRUE))

preds_i0_olink <- tibble(f = l[str_detect(l, "male")]) %>%
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

preds_i0_olink %>%
  left_join(covs) %>%
  filter(age_recruitment > 39) %>%
  ggplot(aes(x = age_recruitment, y = res, color = sex )) + geom_smooth()

