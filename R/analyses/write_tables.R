

### OUTPUT TABLES

### 


results <- readRDS("data_share/association_results_disease_CA.rds") %>%
  bind_rows(readRDS("data_share/association_results_disease_CM.rds"))

fields <- data.table::fread("data/field.tsv")

r <- results %>%
  filter(term %in% c("res", "abs(res)")) %>%
  mutate(expo = case_when(term == "res" ~ "Circadian Acceleration",
                          term == "abs(res)" ~ "Circadian Misalignment"),
         disease = outcome,
         outcome = str_remove(outcome, "-0.0_prevalent"),
         outcome = str_remove(outcome, "p"),
         field_id = as.numeric(str_remove(outcome, "_prevalent"))) %>%
  left_join(fields %>% select(field_id, title)) %>%
  mutate(family = str_extract(title, "(?<=Date )\\S+(?= first reported)"),
         family = str_sub(family, 1, 2),
         disorder = sub(".*\\((.*)\\).*", "\\1", title),
         disorder = str_remove(disorder, "\\)"),
         disorder = str_to_sentence(disorder)) %>%
  filter(!field_id  %in% c(130898, 130902, 130932, 130944, 130852)) %>%
  filter(!family %in% c("F4", "F5", "F6", "F1" )) %>%
  distinct(disease, disorder, field_id, model, expo, .keep_all = TRUE)

labels <- r %>%
  distinct(field_id, title, family, disease, disorder, n) %>%
  arrange(desc(n))

t_prev <- readRDS("data_share/stats_prevalent.rds") %>%
  inner_join(labels) %>%
  select(family, disorder, field_id, N, contains("age"), contains("sex"), contains("bmi"), contains("smoke"))


write_csv(t_prev, "tables/ST_prevalent.csv")

results <- readRDS("data_share/association_results_cox_disease_CA.rds") %>%
  bind_rows(readRDS("data_share/association_results_cox_disease_CM.rds"))

fields <- data.table::fread("data/field.tsv")

r <- results %>%
  filter(term %in% c("res", "abs(res)")) %>%
  mutate(expo = case_when(term == "res" ~ "Circadian Acceleration",
                          term == "abs(res)" ~ "Circadian Misalignment"),
         disease = paste0(outcome, "_incident"),
         disorder = str_remove(disorder, "\\)")) %>%
  filter(!field_id  %in% c(130898, 130902, 130932, 130944, 130852)) %>%
  filter(!family %in% c("F4", "F5", "F6", "F1" )) %>%
  distinct(disease, disorder, field_id, model, expo, .keep_all = TRUE)

labels <- r %>%
  distinct(field_id, title, family, disease, disorder, n) %>%
  arrange(desc(n))

t_inc <- readRDS("data_share/stats_incident.rds") %>%
  inner_join(labels) %>%
  select(family, disorder, field_id, N, contains("age"), contains("sex"), contains("bmi"), contains("smoke"))

write_csv(t_inc, "tables/ST_incident.csv")
