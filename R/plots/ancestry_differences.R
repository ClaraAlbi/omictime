





pred <- readRDS("/mnt/project/olink_int_panels14.rds") %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds")) %>%
  #filter(i == 0) %>%
  filter(!is.na(time_day)) %>%
  mutate(date = as.POSIXct(date_bsampling, tz = "Europe/London"))

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
         assessment_centre = as.factor(assessment_centre)
  )

a <- data.table::fread("/mnt/project/ancestry_new.csv") %>%
  mutate(p30079 = case_when(p30079 == "" ~ NA_character_,
                            TRUE ~ p30079),
         p30079 = relevel(as.factor(p30079), ref = "European ancestry (EUR)"))

job_vars <- data.table::fread("/mnt/project/job_vars.tsv") %>%
  mutate(night_shift = case_when(`3426-0.0` == 1 ~ "Never",
                                 `3426-0.0` == 2 ~ "Sometimes",
                                 `3426-0.0` == 3 ~ "Usually",
                                 `3426-0.0` == 4 ~ "Always"),
         shift_work = case_when(`826-0.0` == 1 ~ "Never/rarely",
                                `826-0.0` == 2 ~ "Sometimes",
                                `826-0.0` == 3 ~ "Usually",
                                `826-0.0` == 4 ~ "Always"),
         employment = case_when(`6142-0.0` == 1	~ "Employed",
                                `6142-0.0` == 2	~ "Retired",
                                `6142-0.0` == 3 ~ "Home/family manager",
                                `6142-0.0` == 4 ~	"Disability",
                                `6142-0.0` == 5 ~	"Unemployed",
                                `6142-0.0` == 6 ~	"Voluntary work",
                                `6142-0.0` == 7 ~	"Student")) %>%
  mutate(night_shift = factor(night_shift, levels = c("Never", "Sometimes", "Usually", "Always")),
         shift_work = factor(shift_work, levels = c("Never/rarely", "Sometimes", "Usually", "Always")),
         employment = fct_relevel(employment, "Employed", "Retired"))


data <- pred %>%
  left_join(sleep) %>%
  left_join(job_vars) %>%
  left_join(covs)


data$res <- residuals(lm(pred_mean ~ time_day, data = data))

data %>% mutate(t = round(time_day, 0)) %>% filter(t == 15) %>% arrange(res, chrono) %>% slice(1:5) %>% pull(res, chrono)

data %>% arrange(res, chrono) %>% slice(1:5) %>% pull(res, chrono)

library(ggplot2)

data %>%
  ggplot(aes(x = age_recruitment, y = pred_mean, color = sex)) +
  geom_smooth()


tidy(lm(pred_mean ~ age_recruitment*sex, data = data))


d <- data %>%
  left_join(a) %>%
  group_by(cv, p30079) %>% nest() %>%
  filter(!is.na(p30079)) %>%
  mutate(mod = map_dbl(data, ~cor(.x$time_day, .x$pred_mean)^2)) %>%
  select(-data)

saveRDS(d, "data_share/prediction_by_ancestry.rds")
