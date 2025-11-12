
library(purrr)
install.packages("broom")
library(broom)

covs <- readRDS("/mnt/project/clara/covs.rds")
gencovs <- readRDS("/mnt/project/clara/gencovs.rds")

sleep1 <- data.table::fread("/mnt/project/Safia/Sleep_Assessment_Centre_All.csv")

data <- readRDS("/mnt/project/clara/proteomic_time_prediction_results.rds")
data$res <- residuals(lm(pred_lasso ~ time_day, data = data))
sleep <- data.table::fread("/mnt/project/Safia/Sleep_Online_All.csv")

data_m <- data %>%
  left_join(sleep)

summaries <- data_m %>%
  select(-p30489, -p30488, -p30491, -p30490, -p32122, -p32121, -p32120, -p30494, -p30493, -p30492, -p32124, -p32123, -p32126, -p32125, -p32128, -p32127, -p32129, -p32130, -p32132, -p32131, -p32134, -p32133) %>%
  pivot_longer(-c(eid, fasting, time_day, pred_lasso, res), names_to = "name") %>%
  filter(value != "") %>%
  group_by(name) %>%
  count(value) %>%
  mutate(field_id = as.numeric(str_remove(name, "p")))


fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

pall <- summaries %>%
  left_join(fields %>% select(field_id, title, main_category))

nested_d <- data_m %>%
  select(eid, res, contains(paste0("p", pall$field_id[pall$main_category == 209]))) %>%
  left_join(covs %>% select(eid, sex, age_recruitment)) %>%
  left_join(gencovs %>% select(eid, paste0("PC", 1:10))) %>%
  pivot_longer(starts_with("p",ignore.case = F)) %>%
  filter(value != "") %>%
  filter(!is.na(value)) %>%
  group_by(name, value) %>%
  mutate(n = n()) %>% filter(n > 20) %>%
  ungroup() %>% group_by(name) %>%
  nest()


results <- nested_d %>%
  mutate(m = map(data, ~lm(res ~ as.factor(value) + sex + age_recruitment + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +PC7 + PC8 + PC9 + PC10, data = .x)),
         t = map(m, broom::tidy)) %>%
  select(name, t) %>%
  unnest()


resutls <- nested_d %>%
  mutate(m = map(data, ~lm(res ~ recode(factor(value), reference = "Definitely a 'morning' person") , data = .x)),
         t = map(m, broom::tidy)) %>%
  select(name, t) %>%
  unnest()

data_m$chronotype <- factor(data_m$p1180_i0)

data_m$chronotype = factor(data_m$chronotype, levels = c("Definitely a 'morning' person", "More a 'morning' than 'evening' person", "Do not know", "More an 'evening' than a 'morning' person", "Definitely an 'evening' person"))

summary(lm(res ~ chronotype, data = data_m))
