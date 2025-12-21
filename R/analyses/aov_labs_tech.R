library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(lubridate)

install.packages("janitor")
install.packages("broom")
install.packages('bit64')

type <- "labs_tech"

bioquem_cov <- data.table::fread("/mnt/project/biomarkers_assay_date.tsv") %>%
  janitor::clean_names()

bioquem <- data.table::fread("/mnt/project/biomarkers_2.tsv") %>%
  janitor::clean_names()

time <- readRDS("/mnt/project/biomarkers/time.rds")

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
         assessment_centre = as.factor(assessment_centre))

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  janitor::clean_names() %>%
  select(eid, x22000_0_0, x22006_0_0, x22021_0_0, x22009_0_1:x22009_0_20)
colnames(gen_covs) <- c("eid", "genotype_batch", "is_white", "rel", paste0("PC", 1:20))

dilution_factor <- data.table::fread("/mnt/project/dilution_factor.tsv") %>%
  rename(dilution_factor = 2)


#### MAKE TRANS FOR EACH BIOMARK INDIVIDUALLY
l2 <- list()
l3 <- list()

for (i in 2:ncol(bioquem)) {
  print(i)
  var <- colnames(bioquem)[i]
  numeric_name <- as.numeric(str_remove(str_remove(var, "_0_0"), "x"))
  assay_date <- numeric_name + 1

  sel <- paste0("x", c(assay_date), "_0_0")

  covs_c <- bioquem_cov %>% select(eid, any_of(sel))
  colnames(covs_c)  <- c("eid", "assay_date")

  d <- tibble(eid = bioquem$eid, raw = bioquem[[i]]) %>%
    drop_na() %>%
    left_join(gen_covs %>% select(eid, any_of(paste0("PC", 1:20)))) %>%
    left_join(covs %>% select(eid, sex, age_recruitment, assessment_centre, month_attending, bmi, smoking)) %>%
    left_join(covs_c) %>%
    left_join(time %>% select(eid, fasting, time_day)) %>%
    left_join(dilution_factor) %>%
    mutate(assay_date = as.factor(ntile(assay_date, 20)),
           dilution_factor = as.factor(ntile(dilution_factor, 20))) %>%
    filter(time_day > 0) %>%
    filter(fasting < 24) %>%
    mutate(across(c(sex, assessment_centre, age_recruitment, month_attending, smoking), as.factor),
           log_b = log(raw)) %>%
    select(eid,
           raw,
           log_b,
           time_day,
           fasting,
           assessment_centre,
           month_attending,
           any_of(paste0("PC", 1:20)),
           dilution_factor,
           bmi,
           smoking) %>%
    filter(if_all(where(is.factor), ~ . %in% names(which(table(.) >= 50)))) %>%
    mutate(across(where(is.factor), droplevels)) %>%
    filter(time_day >= 9 & time_day <= 20)

  m <- lm(log_b ~ ., data = d %>% select(-eid, -raw, -time_day, -bmi, -smoking), na.action = "na.exclude")
  res_rint <- residuals(m)

  m_time_only_cos <- lm(res_rint ~ cos(2 * pi * time_day / 24) + sin(2 * pi * time_day / 24), data = d)

  out <- tibble(eid = d$eid, raw = d$raw, res = res_rint)
  l2 <- c(l2, list(out))

  effects <- broom::tidy(m_time_only_cos , conf.int = TRUE) %>% mutate(phen = numeric_name)
  l3 <- c(l3, list(effects))

}


output_effects <- do.call(bind_rows, l3)
saveRDS(output_effects, glue("effects_{type}.rds"))

output_res <- tibble(res = l2) %>%
  mutate(phen = as.numeric(str_remove(str_remove(colnames(bioquem)[-1], "_0_0"), "x"))) %>%
  unnest(res) %>%
  pivot_wider(id_cols = eid, names_from = phen, values_from = res)

saveRDS(output_res, glue("res_{type}.rds"))
