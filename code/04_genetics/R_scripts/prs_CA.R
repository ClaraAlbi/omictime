source("config/paths.R")

library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(broom)
library(data.table)
library(lubridate)

ukbppp <- readRDS(project_path("biomarkers_res/olink_int_replication_v2.rds")) # %>%
  filter(i == 0)

a <- data.table::fread(project_path("ancestry_new.csv")) %>%
  filter(p30079 == "European ancestry (EUR)")

covs <- readRDS(project_path("biomarkers/covs.rds")) %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
         assessment_centre = as.factor(assessment_centre)
  )

job_vars <- data.table::fread(project_path("job_vars.tsv")) %>%
  mutate(night_shift = case_when(`3426-0.0` == 1 ~ "Never",
                                 `3426-0.0` == 2 ~ "Sometimes",
                                 `3426-0.0` == 3 ~ "Usually",
                                 `3426-0.0` == 4 ~ "Always"),
         shift_work = case_when(`826-0.0` == 1 ~ "Never/rarely",
                                `826-0.0` == 2 ~ "Sometimes",
                                `826-0.0` == 3 ~ "Usually",
                                `826-0.0` == 4 ~ "Always")) %>%
  mutate(night_shift = factor(night_shift, levels = c("Never", "Sometimes", "Usually", "Always")),
         shift_work = factor(shift_work, levels = c("Never/rarely", "Sometimes", "Usually", "Always")))

gen_covs <- data.table::fread(project_path("genetic_covs.tsv")) %>%
  select(eid, "22009-0.1":"22009-0.20", `22006-0.0`, `22021-0.0`, `22000-0.0`)
colnames(gen_covs) <- c("eid", paste0("PC", 1:20), "is_white", "rel", "batch")

a <- data.table::fread(project_path("ancestry_new.csv")) %>%
  mutate(p30079 = case_when(p30079 == "" ~ NA_character_,
                            TRUE ~ p30079),
         p30079 = relevel(as.factor(p30079), ref = "European ancestry (EUR)"))

sleep <- data.table::fread(project_path("chronotype2.tsv")) %>%
  select(eid,
         h_sleep = `1160-0.0`,
         chrono = `1180-0.0`,
         ever_insomnia = `1200-0.0`,
         wakeup = `1170-0.0`,
         snoring = `1210-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ "Definitely morning",
    chrono == 2 ~ "Rather morning",
    chrono == -1~ "Don't know",
    chrono == 3 ~ "Rather evening",
    chrono == 4 ~ "Definitely evening",
    TRUE ~ NA_character_),
    wakeup = case_when(
      wakeup == 1 ~ "Not at all easy",
      wakeup == 2 ~ "Not very easy",
      wakeup == -1 ~ "Do not know",
      wakeup == 3 ~ "Fairly easy",
      wakeup == 4 ~ "Very easy"),
    chrono = factor(chrono, levels = c("Definitely morning", "Rather morning", "Don't know", "Rather evening", "Definitely evening")),
    wakeup = factor(wakeup, levels = c("Very easy", "Fairly easy",  "Not very easy", "Not at all easy")))


dis2 <- bind_cols(fread(project_path("top_diseases_IEFG.csv")), fread(project_path("immune.csv")) %>% select(-eid)) %>%
  mutate(across(-eid, as_date))


dis_prev <- dis2 %>%
  select(eid, 1:10,
         contains("130894"),
         contains("130814"),
         contains("131848")) %>%
  left_join(
    readRDS(project_path("biomarkers/time.rds")) %>%
      select(eid, date_bsampling),
    by = "eid"
  ) %>%
  mutate(
    date_bsampling = as.Date(date_bsampling),
    across(contains("13"), as.Date)
  ) %>%
  mutate(
    across(
      contains("13"),
      list(
        prev = ~ if_else(!is.na(.x) & .x < date_bsampling, 1L, 0L),
        prev_lt5y = ~ if_else(
          !is.na(.x) &
            .x < date_bsampling &
            .x >= (date_bsampling - years(5)),
          1L, 0L
        )
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  select(eid, contains("prev"))



prs <- data.table::fread(project_path("PRS_CA_allchr.profile")) %>% mutate(eid = as.numeric(IID)) %>% select(eid, PRS_SUM) %>%
  left_join(data.table::fread(project_path("prs_prots_allchr.sscore")) %>% mutate(eid = as.numeric(IID)) %>% select(-FID, -IID), by = c("eid")) %>%
  filter(!eid %in% ukbppp$eid) %>%
  #¢filter(eid %in% a$eid) %>%
  mutate(CA_prs = scale(PRS_SUM)[,1]) %>%
  mutate(across(ends_with("AVG"), ~scale(.x)[,1])) %>%
  left_join(sleep) %>%
  left_join(covs) %>%
  left_join(gen_covs) %>%
  left_join(job_vars) %>%
  left_join(dis_prev)

hist(prs$CA_prs)

res <- prs %>%
  pivot_longer(starts_with("p13")) %>%
  group_by(name) %>%
  nest() %>%
  mutate(d = map(data, ~tidy(lm(paste0("value ~ ", paste0(c("CA_prs", "sex","age_recruitment" , paste0("PC", 1:20)), collapse = " + ")), data = .x)))) %>%
  select(name, d) %>% unnest(d)

""
m <- lm(paste0("CA_prs ~ ", paste0(c("chrono", "sex","age_recruitment" , paste0("PC", 1:20)), collapse = " + ")), data = prs)
tidy(m)
summary(m)
summary(lm(paste0("CA_prs ~ ", paste0(c("wakeup", "sex","age_recruitment" , paste0("PC", 1:20)), collapse = " + ")), data = prs))

summary(lm(CA_prs ~ night_shift, data = prs))
summary(lm(CA_prs ~ sex, data = prs))


tidy(lm(paste0(" ~ ", paste0(c("CA_prs", "sex","age_recruitment" , paste0("PC", 1:20)), collapse = " + ")), data = prs))





# PREDICTION INTERNAL
df <- ukbppp %>%
  left_join(a) %>%
  left_join(data.table::fread("PRS_CA_norep_allchr.profile") %>% mutate(eid = as.numeric(IID)) %>% select(eid, PRS_SUM) ) %>%
  filter(eid %in% a$eid)

df$res <- residuals(lm(pred_mean ~ time_day, data= df))

df %>%
  group_by(p30079) %>% nest() %>%
  filter(!is.na(p30079)) %>%
  mutate(r2_prs = map_dbl(data, ~cor(.x$res, .x$PRS_SUM, use = "complete.obs")^2),
         n = map_dbl(data, nrow)) %>%
  select(-data)

df$CA_prs<- scale(df$PRS_SUM)[,1]



df_i0 <- df %>% filter(n > 1) %>% filter(i == 0)

p0 <- cor(df_i0$res, df_i0$CA_prs, use = "complete.obs")^2
df_i2 <- df %>% filter(i == 2)

p2<-cor(df_i2$res, df_i2$CA_prs, use = "complete.obs")^2
df_i3 <- df %>% filter(i == 3)

p3<- cor(df_i3$res, df_i3$CA_prs, use = "complete.obs")^2


### Prediction


prs <- data.table::fread("PRS_CA_norep_allchr.profile") %>%  mutate(FID = as.integer(FID)) %>% select(eid = FID, PRS_SUM) %>%
  filter(eid %in% a$eid) %>%
  mutate(CA_prs = scale(PRS_SUM)[,1]) %>%
  left_join(sleep) %>%
  filter(!eid %in% ukbppp$eid)


tidy(lm(CA_prs ~ chrono, data = prs))

prs %>% group_by(chrono) %>% count()
