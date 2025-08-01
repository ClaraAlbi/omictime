install.packages("broom")
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)

time <- readRDS("/mnt/project/biomarkers/time.rds")
covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>% mutate(bmi = weight/((height/100)^2)) %>%
  left_join(time %>% select(eid, fasting)) %>%
  filter(fasting < 24)
gen_covs <- fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20")
colnames(gen_covs) <- c("eid", paste0("PC", 1:20))

data_b <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(i == 0) %>%
  left_join(covs) %>% left_join(gen_covs)
data_b$res <- residuals(lm(pred_lasso ~ time_day, data = data_b))
data_b$res_abs <- abs(data_b$res)
data_b$gap <- data_b$pred_lasso - data_b$time_day
data_b$gap_abs <- abs(data_b$gap)
data_b$ares_q <- ntile(data_b$res_abs, 5)
data_b$ares_q <- factor(data_b$ares_q, levels = c(1:5))


dis2 <- readRDS("diseases_circadian.rds")

prev_df <- dis2 %>%
  mutate(across(-c(eid, date_bsampling), ~case_when(.x < date_bsampling ~ 1,
                                                       TRUE ~ 0)))

Absgap <- prev_df %>%
  pivot_longer(-c(eid, date_bsampling)) %>%
  group_by(name) %>%
  nest() %>%
  ungroup() %>%
  mutate(data = map(data, ~.x %>%
                      inner_join(data_b) %>%
                      filter(!is.na(res)) %>%
                      filter(!is.na(value)) %>%
                      filter(if_all(where(is.factor), ~ . %in% names(which(table(.) >= 50)))) %>%
                      mutate(across(where(is.factor), droplevels))),
         mod1 = map(data, ~broom::tidy(glm(paste0("value ~ gap_abs + sex + age_recruitment + ", paste0("PC",1:20, collapse = "+")), data = .x, family = "binomial"), conf.int = FALSE)),
         mod2 = map(data, ~broom::tidy(glm(paste0("value ~ gap_abs + sex + age_recruitment + smoking + bmi + fasting +", paste0("PC",1:20, collapse = "+")), data = .x, family = "binomial"), conf.int = TRUE)),
         n = map_dbl(data, ~sum(.x$value == 1, na.rm = T))) %>%
  select(name, mod1, mod2, n)


Absgap %>%
  filter(n >= 50) %>%
  mutate(
    est1 = map_dbl(mod1, ~ .x %>% filter(term == "gap_abs") %>% pull(estimate)),
    est2 = map_dbl(mod2, ~ .x %>% filter(term == "gap_abs") %>% pull(estimate)),
    se1  = map_dbl(mod1, ~ .x %>% filter(term == "gap_abs") %>% pull(std.error)),
    se2  = map_dbl(mod2, ~ .x %>% filter(term == "gap_abs") %>% pull(std.error))) %>%
  select(name, est1, est2, se1, se2) %>%
  pivot_longer(
    cols = c(est1, est2, se1, se2),
    names_to = c(".value", "model"),
    names_pattern = "(est|se)(\\d)"
  ) %>%
  mutate(model = recode(model, "1" = "mod1", "2" = "mod2")) %>%
  left_join(list_diseases, by = c("name" = "p")) %>%
  ggplot(aes(x = exp(est), y = names, color = model)) +
  geom_point() +
  facet_grid(rows = vars(class), scales = "free") +
  theme_minimal() +
  theme(legend.position = "none")
