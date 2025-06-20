library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(lubridate)

install.packages("janitor")
install.packages("broom")
install.packages('bit64')

type <- "labs"

bioquem_cov <- data.table::fread("/mnt/project/biomarkers_assay_date.tsv") %>%
  janitor::clean_names()

bioquem <- data.table::fread("/mnt/project/biomarkers_2.tsv") %>%
  janitor::clean_names()

time <- readRDS("/mnt/project/biomarkers/time.rds")

covs <- data.table::fread("/mnt/project/circadian/data/covariates.tsv") %>%
  janitor::clean_names() %>%
  mutate(across(c(2, 7), as.factor),
         bmi = `x21002_0_0` / (`x50_0_0`/100)^2)
colnames(covs) <- c("eid", "sex", "birth_year", "age_recruitment", "height", "assessment_centre", "month_attending", "weight", "smoking", "bmi")

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  janitor::clean_names() %>%
  select(eid, x22000_0_0, x22006_0_0, x22021_0_0, x22009_0_1:x22009_0_20)
colnames(gen_covs) <- c("eid", "genotype_batch", "is_white", "rel", paste0("PC", 1:20))

dilution_factor <- data.table::fread("/mnt/project/dilution_factor.tsv") %>%
  rename(dilution_factor = 2)


#### MAKE TRANS FOR EACH BIOMARK INDIVIDUALLY
l1 <- list()
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
    filter(raw != 0) %>%
    mutate(smoking = case_when(smoking == -3 ~ NA_integer_,
                               TRUE ~ smoking),
           across(c(sex, assessment_centre, age_recruitment, month_attending, smoking), as.factor),
           log_b = log(raw)) %>%
    select(eid,
           raw,
           log_b,
           time_day,
           sex, 
           age_recruitment, 
           fasting,
           assessment_centre,
           month_attending,
           any_of(paste0("PC", 1:20)),
           dilution_factor,
           bmi,
           smoking) %>%
    filter(if_all(where(is.factor), ~ . %in% names(which(table(.) >= 50)))) %>%
    mutate(across(where(is.factor), droplevels)) %>%
    filter(time_day >= 8 & time_day <= 20)
  
  mod <- aov(log_b ~ . + sex * age_recruitment, data = d %>% select(-eid, -raw))
  out_aov <- broom::tidy(mod) %>% mutate(pr2 = sumsq / sum(sumsq),
                                         phen = numeric_name)
  l1 <- c(l1, list(out_aov))
  
  m <- lm(log_b ~ . + sex * age_recruitment, data = d %>% select(-eid, -raw, -time_day, -bmi, -smoking), na.action = "na.exclude")
  res_rint <- residuals(m)
  
  m_time_only_cos <- lm(res_rint ~ cos(2 * pi * time_day / 24) + sin(2 * pi * time_day / 24), data = d)
  
  out <- tibble(eid = d$eid, raw = d$raw, res = res_rint)
  l2 <- c(l2, list(out))
  
  effects <- broom::tidy(m_time_only_cos , conf.int = TRUE) %>% mutate(phen = numeric_name)
  l3 <- c(l3, list(effects))
  
}

output_aov <- do.call(bind_rows, l1)
saveRDS(output_aov, glue("aov_{type}.rds"))

output_effects <- do.call(bind_rows, l3)
saveRDS(output_effects, glue("effects_{type}.rds"))

output_res <- tibble(res = l2) %>%
  mutate(phen = as.numeric(str_remove(str_remove(colnames(bioquem)[-1], "_0_0"), "x"))) %>%
  unnest(res) %>%
  pivot_wider(id_cols = eid, names_from = phen, values_from = res)

saveRDS(output_res, glue("res_{type}.rds"))

output_raw <- tibble(res = l2) %>%
  mutate(phen = as.numeric(str_remove(str_remove(colnames(bioquem)[-1], "_0_0"), "x"))) %>%
  unnest(res) %>%
  pivot_wider(id_cols = eid, names_from = phen, values_from = raw)

saveRDS(output_raw, glue("raw_{type}.rds"))


