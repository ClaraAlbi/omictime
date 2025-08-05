library(lubridate)
library(stringr)
library(dplyr)
library(dplyr)


olink <- readRDS("/mnt/project/olink_int_replication.rds") %>% filter(i == 0)
time <- readRDS("/mnt/project/biomarkers/time.rds")

list_diseases <- tribble(~p, ~names, ~class,
        "p130896", "Recurrent depression", "1. Neuro-psychiatric",
        "p130894", "Depressive episode", "1. Neuro-psychiatric",
        "p130892", "Bipolar disorder", "1. Neuro-psychiatric",
        "p130906", "Anxiety", "1. Neuro-psychiatric",
        "p130920", "Sleep - F51", "1. Neuro-psychiatric",
        "p131060", "Sleep - G47", "1. Neuro-psychiatric",
        "p131052", "Migraine", "1. Neuro-psychiatric",
        "p131048", "Epilepsy", "1. Neuro-psychiatric",
        "p130842", "Dementia", "1. Neuro-psychiatric",
        "p130836", "Alzheimer's disease (dementia)", "1. Neuro-psychiatric",
        "p131022", "Parksinson's disease", "1. Neuro-psychiatric",
        "p131286", "Hypertension", "2. Cardiometabolic",
        "p131298", "Acute Myocardial Infarction","2. Cardiometabolic",
        "p131380", "Atherosclerosis", "2. Cardiometabolic",
        "p131306", "Chronic ischaemic heart disease", "2. Cardiometabolic",
        "p130708", "Type 2 Diabetes", "2. Cardiometabolic",
        "p130792", "Obesity", "2. Cardiometabolic",
        "p130814", "Disorders of lipoprotein metabolism", "2. Cardiometabolic",
        "p131668", "Inflammatory liver diseases", "2. Cardiometabolic",
        "p131848", "Rheumatoid arthritis", "3. Immune",
        "p131494", "Asthma", "3. Immune",
        "p131626", "Crohns's disease", "3. Immune",
        "p131628", "Ulcerative colitis", "3. Immune",
        "p131742", "Psoriasis", "3. Immune",
        "date_breast", "Breast cancer", "4. Cancer",
        "date_lung", "Lung cancer","4. Cancer",
        "date_melanoma", "Melanoma", "4. Cancer")


# CANCER DIAGNOSES
cancer <- fread("/mnt/project/cancer.csv")

long_cancer <- cancer %>%
  pivot_longer(
    cols = matches("^(p40006|p40005)_i\\d+$"),
    names_to = c(".value", "idx"),
    names_pattern = "(p40006|p40005)_i(\\d+)"
  ) %>% rename(code = p40006, date = p40005)
long_cancer$code[which(long_cancer$code == "")] <- NA

target_codes <- c("C34", "C43", "C50")
filtered <- long_cancer %>%
  filter(str_detect(code, paste(target_codes, collapse = "|")))
labeled <- filtered %>%
  mutate(cancer = case_when(
    str_detect(code, "C34") ~ "lung",
    str_detect(code, "C43") ~ "melanoma",
    str_detect(code, "C50") ~ "breast"
  ))
earliest <- labeled %>%
  filter(!is.na(date)) %>%
  group_by(eid, cancer) %>%
  summarise(earliest_date = min(date), .groups = "drop") %>%
  pivot_wider(
    names_from = cancer,
    values_from = earliest_date,
    names_prefix = "date_"
  )

# NON-CANCER DIAGNOSES
dis2 <- bind_cols(fread("/mnt/project/top_diseases_IEFG.csv"), fread("immune.csv") %>% select(-eid)) %>%
  mutate(across(-eid, as_date))  %>%
  select(eid,any_of(list_diseases$p)) %>%
  left_join(earliest) %>%
  left_join(time %>% select(eid, date_bsampling), by = "eid") %>%
  filter(eid %in% olink$eid)

saveRDS(dis2, "diseases_circadian.rds")

# Total counts
fields <- fread("/mnt/project/Showcase metadata/field.tsv")


total_df <- dis2 %>%
  summarise(across(-c(eid, date_bsampling), ~sum(!is.na(.x)), .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "disease", values_to = "total")

# 3. Compute prevalent counts (< sampling date)
prev_df <- dis2 %>%
  summarise(across(-c(eid, date_bsampling), ~sum(.x < date_bsampling, na.rm = TRUE), .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "disease", values_to = "prevalent")

# 4. Compute incident counts (> sampling date)
inc_df <- dis2 %>%
  summarise(across(-c(eid, date_bsampling), ~sum(.x > date_bsampling, na.rm = TRUE), .names = "{.col}")) %>%
  pivot_longer(everything(), names_to = "disease", values_to = "incident")

# 5. Join all three
cs_pi <- total_df %>%
  left_join(prev_df, by = "disease") %>%
  left_join(inc_df, by = "disease") %>%
  mutate(
    type = map_chr(disease, function(x) {
      if (x %in% c("date_breast", "date_melanoma", "date_lung")) {x} else {
        fields$title[paste0("p", fields$field_id) == x][1] %||% x }}))

saveRDS(cs_pi, "counts_diseases_circadian.rds")


cs_pi <- readRDS("/mnt/project/counts_diseases_circadian.rds") %>%
  left_join(list_diseases, by = c("disease" = "p")) %>%
  select(disease, names, type, class, total, prevalent, incident)

library(data.table)
fwrite(cs_pi, "data_share/variables_diseases_circadian.csv")

