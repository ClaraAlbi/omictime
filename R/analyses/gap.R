
install.packages("broom")

covs <- readRDS("/mnt/project/biomarkers/covs.rds")

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20")
colnames(gen_covs) <- c("eid", paste0("PC", 1:20))

d2 <- d %>% left_join(covs) %>%
  left_join(gen_covs)

time <- readRDS("/mnt/project/biomarkers/time.rds")

fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

outcomes <- tribble(~field_id, ~phen,
                    "130894", "Depressive_episode",
                    "130896", "Recurrent_depression",
                    "130892", "Bipolar",
                    "130874", "Schizophrenia",
                    "130708", "T2D",
                    "130792", "Obesity",
                    "131286", "Hypertension",
                    "131306", "IHD",
                    "131380", "Atherosclerosis",
                    "131060", "Sleep disorders")

dis2 <- inner_join(data.table::fread("/mnt/project/vars_diseases_2.tsv"), data.table::fread("/mnt/project/vars_diseases.tsv")) %>%
  select(eid, contains(outcomes$field_id)) %>%
  left_join(time) %>%
  mutate(across(starts_with("13"), ~ case_when(.x < date_bsampling ~ 1,
                                               is.na(.x) ~ 0), .names = "{.col}_prevalent"),
         across(starts_with("13"), ~ case_when(.x > date_bsampling ~ 1,
                                               is.na(.x) ~ 0), .names = "{.col}_incident")) %>%
  select(-contains("_prevalent_incident")) %>%
  select(eid, contains("prevalent"))


bp <- data.table::fread("/mnt/project/vars_diseases_2.tsv") %>%
  select(eid, contains("130892"))




Absgap <- dis2 %>%
  pivot_longer(-eid) %>%
  group_by(name) %>%
  nest() %>%
  ungroup() %>%
  mutate(field_id = str_extract(name, "^\\d+")) %>%
  left_join(outcomes) %>%
  #filter(field_id == "130892") %>%
  mutate(data = map(data, ~inner_join(.x, d2) %>%
                      filter(!is.na(value)) %>%
                      filter(if_all(where(is.factor), ~ . %in% names(which(table(.) >= 50)))) %>%
                      mutate(across(where(is.factor), droplevels))),
         mod = map(data, ~broom::tidy(glm(value ~ absgap + time_day + sex + age_recruitment + smoking, data = .x %>% mutate(absgap = abs(res)) %>% select(-eid, -res), family = "binomial"), conf.int = F)),
         n = map_dbl(data, ~sum(.x$value == 1, na.rm = T))) %>%
  select(field_id, phen, mod, n) %>%
  unnest(mod)

saveRDS(Absgap, "disease_associations_absgap_ci.rds")

Gap <- dis2 %>%
  pivot_longer(-eid) %>%
  group_by(name) %>%
  nest() %>%
  ungroup() %>%
  mutate(field_id = str_extract(name, "^\\d+")) %>%
  left_join(outcomes) %>%
  #filter(field_id != "130792") %>%
  mutate(data = map(data, ~inner_join(.x, d2) %>%
                      filter(!is.na(value)) %>%
                      filter(if_all(where(is.factor), ~ . %in% names(which(table(.) >= 50)))) %>%
                      mutate(across(where(is.factor), droplevels))),
         mod = map(data, ~broom::tidy(glm(value ~ res + time_day + sex + age_recruitment + smoking, data = .x , family = "binomial"), conf.int = F)),
         n = map_dbl(data, ~sum(.x$value == 1, na.rm = T))) %>%
  select(field_id, phen, mod, n) %>%
  unnest(mod)

saveRDS(Gap, "disease_associations_gap_ci.rds")
