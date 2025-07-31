
library(stringr)

time <- readRDS("/mnt/project/biomarkers/time.rds")

dis2 <- fread("/mnt/project/top_diseases_IEFG.csv") %>%
  mutate(across(-eid, ~ ifelse(.x > "1902-02-02", as.Date(.x), NA)))  # keep only real dates

fields <- fread("/mnt/project/Showcase metadata/field.tsv")

cs <- dis2 %>%
  summarise(across(-eid, ~sum(!is.na(.x)))) %>%
  pivot_longer(everything()) %>%
  mutate(type = map_chr(name, ~fields$title[paste0("p",fields$field_id) == .x])) %>%
  filter(!str_detect(type, "other|elsewhere|pulm"))


cs_pi <- dis2 %>%
  left_join(time %>% select(eid, date_bsampling)) %>%
  mutate(across(-c(eid, date_bsampling), ~ case_when(.x < date_bsampling ~ 1,
                                               is.na(.x) ~ 0), .names = "{.col}_prevalent"),
         across(-c(eid, date_bsampling), ~ case_when(.x > date_bsampling ~ 1,
                                               is.na(.x) ~ 0), .names = "{.col}_incident")) %>%
  select(-contains("_prevalent_incident")) %>%
  summarise(across(-c(eid, date_bsampling), ~sum(!is.na(.x)))) %>%
  pivot_longer(everything())

