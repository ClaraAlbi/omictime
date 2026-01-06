
ukbppp <- readRDS("/mnt/project/olink_internal_time_predictions.rds") %>%
  filter(i == 0)

a <- data.table::fread("/mnt/project/ancestry_new.csv") %>%
  filter(p30079 == "European ancestry (EUR)")

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
         assessment_centre = as.factor(assessment_centre)
  )

job_vars <- data.table::fread("/mnt/project/job_vars.tsv") %>%
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

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20", `22006-0.0`, `22021-0.0`, `22000-0.0`)
colnames(gen_covs) <- c("eid", paste0("PC", 1:20), "is_white", "rel", "batch")


sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
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
    chrono = factor(chrono, levels = c("Definitely morning", "Rather morning", "Don't know", "Rather evening", "Definitely evening")))

prs <- data.table::fread("PRS_CA_allchr.profile") %>% mutate(eid = as.numeric(IID)) %>%
  filter(!FID %in% ukbppp$eid) %>%
  filter(FID %in% a$eid) %>%
  mutate(CA_prs = scale(PRS_SUM)[,1]) %>%
  left_join(sleep) %>%
  left_join(covs) %>%
  left_join(gen_covs) %>%
  left_join(job_vars)

hist(prs$CA_prs)

summary(lm(paste0("CA_prs ~ ", paste0(c("chrono", "sex","age_recruitment" , paste0("PC", 1:20)), collapse = " + ")), data = prs))
summary(lm(paste0("CA_prs ~ ", paste0(c("wakeup", "sex","age_recruitment" , paste0("PC", 1:20)), collapse = " + ")), data = prs))

summary(lm(CA_prs ~ night_shift, data = prs))
summary(lm(CA_prs ~ sex, data = prs))
