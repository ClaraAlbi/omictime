library(tidyr)
library(lubridate)
library(dplyr)
library(glue)
library(purrr)
library(ggplot2)
install.packages("forcats")
library(forcats)
install.packages("janitor")
install.packages("broom")
library(stringr)

dx download file-J4y3QpQJk8K76J9q0Y30zVQP

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
         assessment_centre = as.factor(assessment_centre)
  )

pcs <- data.table::fread("/mnt/project/covariates.txt") %>%
  select(eid = 1, contains("PC"))


meds <- data.table::fread("/mnt/project/medications.tsv")

# https://www.sciencedirect.com/science/article/abs/pii/S1389945717302757?via%3Dihub
# Current medications were self-reported to the research nurse, and the participants were dichotomised according to whether they were taking:
#   sleep medication (sedatives and hypnotics),
#   any other psychotropic medication (mood stabilisers, antidepressants, and antipsychotics) or
#   antihypertensive medication (ACE inhibitors, angiotensin II antagonists, beta blockers, calcium channel blockers, and diuretics).

cohort <- meds %>%
  #slice(1:1000) %>%
  pivot_longer(-eid) %>%
  filter(!is.na(value)) %>%
  left_join(X41467_2019_9572_MOESM3_ESM_1_ %>% select(Category, `Medication ATC code`, `Coding a`, `Drug name`), by = c("value" = "Coding a")) %>%
  filter(str_detect(`Medication ATC code`, "N05C") | #sleep medication (sedatives and hypnotics)
           str_detect(`Medication ATC code`, "N03AG01") | # mood stabiliser
           str_detect(`Medication ATC code`, "N03AX09") |  # mood stabiliser
           str_detect(`Medication ATC code`, "N03AF01") |  # mood stabiliser
           str_detect(`Medication ATC code`, "N03AF02") |  # mood stabiliser
           str_detect(`Medication ATC code`, "N05AN01") |  # lithium

           str_detect(`Medication ATC code`, "N06A") | # N06A ANTIDEPRESSANTS


           str_detect(`Medication ATC code`, "C09A") | # ACE inhibitors, plain
           str_detect(`Medication ATC code`, "C09C") | # ANGIOTENSIN II RECEPTOR BLOCKERS (ARBs), PLAIN
           str_detect(`Medication ATC code`, "C07") | # BETA BLOCKING AGENTS
           str_detect(`Medication ATC code`, "C08") | #CALCIUM CHANNEL BLOCKERS
           str_detect(`Medication ATC code`, "C03")) #C03 DIURETICS

cohort_f <- cohort %>%
  mutate(
    N05C = str_detect(`Medication ATC code`, "^N05C"),
    N05A = str_detect(`Medication ATC code`, "^N05A"),
    N03 = str_detect(`Medication ATC code`, "^N03"),
    N06 = str_detect(`Medication ATC code`, "^N06"),
    C0 = str_detect(`Medication ATC code`, "^C0")
  ) %>%
  group_by(eid) %>%
  summarise(
    across(c(N05C, N05A, N03, N06, C0), ~ any(.x, na.rm = TRUE)),
    .groups = "drop"
  )
