library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)
library(stringr)
install.packages("broom")

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         smoking = case_when(smoking == -3 ~NA, TRUE ~ smoking),
         smoking = as.factor(smoking),
         sex = as.factor(sex))

pcs <- data.table::fread("/mnt/project/covariates.txt") %>%
  select(eid = 1, contains("PC"))

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  ungroup()

data <- df %>%
  left_join(covs) %>%
  left_join(pcs) %>%
  mutate(h = round(time_day, 0)) %>%
  filter(age_recruitment > 39)

data$res <- residuals(lm(pred_mean ~ time_day, data = data))


# diagnosis table
sleep2 <- readRDS("/mnt/project/diseases_circadian.rds") #%>% select(eid, p131060, p130708, p130792, p130842, p131306)

outcomes <- tribble(~field_id, ~phen, ~class,
                    "130894", "Depressive_episode","Neuro-psychiatric",
                    "130896", "Recurrent_depression","Neuro-psychiatric",
                    "130892", "Bipolar disorder","Neuro-psychiatric",
                    "130874", "Schizophrenia","Neuro-psychiatric",
                    "130846", "Delirium", "Neuro-psychiatric",
                    "p131060", "Sleep - G47", "Neuro-psychiatric",
                    "130842", "Unspecified dementia", "Neuro-psychiatric",
                    "130836", "Dementia Alzheimer's","Neuro-psychiatric",
                    "p130708", "Type 2 Diabetes", "Cardiometabolic",
                    "p130792", "Obesity", "Cardiometabolic",
                    "p131306", "Ischaemic heart disease", "Cardiometabolic")

dis2 <- data.table::fread("/mnt/project/vars_diseases_2.tsv") %>%
  #select(eid, contains(outcomes$field_id)) %>%
  inner_join(sleep2) %>%
  #left_join(readRDS("/mnt/project/biomarkers/time.rds") %>% select(eid, date_bsampling)) %>%
  mutate(across(contains("13"), ~ case_when(.x < date_bsampling ~ 1,
                                               is.na(.x) ~ 0), .names = "{.col}_prevalent")) %>%
  select(eid, contains("prevalent"))

d_counts <-
  dis2 %>%
  filter(eid %in% data$eid) %>%
  pivot_longer(-eid) %>%
  group_by(name) %>% count(value) %>%
  filter(value == 1 & n > 40)

# Combine
data1 <- data %>%
  select(eid, res, sex, age_recruitment, smoking, bmi, any_of(paste0("PC", 1:10))) %>%
  left_join(dis2 %>% select(eid, contains(d_counts$name)))



vars <- d_counts$name

base_covars   <- c("sex","age_recruitment", paste0("PC", 1:10))

extra_covars <- c("smoking", "bmi")

results <- map_dfr(vars, function(v) {

  # formula for abs(res) only

  f_prev1 <- as.formula(
    paste0("`", v, "` ~ abs(res)")
  )

  # formula for abs(res) + covariates
  f_prev2 <- as.formula(
    paste0("`", v, "` ~ res + ",
           paste(base_covars, collapse = " + ")))
  f_prev3 <- as.formula(
    paste0("`", v, "` ~ res + ",
           paste(c(base_covars, extra_covars), collapse = " + ")))


  # fit models
  m1 <- glm(f_prev1, data = data1, family = binomial)
  m2 <- glm(f_prev2, data = data1, family = binomial)
  m3 <- glm(f_prev3, data = data1, family = binomial)

  # collect results
  bind_rows(
    broom::tidy(m1, conf.int = F, exponentiate = TRUE) %>%
      mutate(model = "Model1", outcome = v),

    broom::tidy(m2, conf.int = F, exponentiate = TRUE) %>%
      mutate(model = "Model2", outcome = v),
    broom::tidy(m3, conf.int = F, exponentiate = TRUE) %>%
      mutate(model = "Model3", outcome = v)
  )
})

saveRDS(results %>%
          left_join(d_counts %>% rename(outcome = name)), "data_share/association_results_disease_CA.rds")



###

library(tidyverse)
library(forcats)

results <- readRDS("data_share/association_results_disease_CA.rds")

fields <- data.table::fread("data/field.tsv")

r <- results %>%
  filter(term == "res") %>%
  mutate(outcome = str_remove(outcome, "-0.0_prevalent"),
         outcome = str_remove(outcome, "p"),
         field_id = as.numeric(str_remove(outcome, "_prevalent"))) %>%
  left_join(fields %>% select(field_id, title))


p_res <- results %>%
  filter(model != "mt_TIME_DAY") %>%
  filter(i == 1) %>%
  mutate(
    # Extract the core field id
    o = str_extract(outcome, "(?<=p_)[0-9]+(?=0\\.0_prevalent)|(?<=p_)p?[0-9]+"),
    model = factor(model, levels =c("m0_MISALIGNMENT", "m1_MISALIGNMENT + COV", "m1_MISALIGNMENT + COV + chrono"),
                   labels = c("Circadian Misalignment", "+ sex + age + 10PCs + BMI + smoking", "+ Chronotype"))
  ) %>%
  left_join(outcomes, by = c("o" = "field_id")) %>%
  mutate(outcome_clean = coalesce(phen, outcome),
         outcome_clean = paste0((outcome_clean), "\nn=", value),
         outcome_clean = fct_reorder(outcome_clean, value)) %>%
  filter(term != "(Intercept)") %>%
  filter(term %in% c("time_day", "abs(res)")) %>%
  ggplot(aes(x = outcome_clean,
             y = estimate,
             ymin = conf.low, ymax = conf.high,
             color = model, shape = model,
             alpha = p.value < 0.05)) +
  geom_pointrange(position = position_dodge(width = 0.6),
                  size = 1, fatten = 3) +
  #geom_text(data = . %>% filter(p.value < 0.05), aes(label = sprintf("%.0e", p.value)),
  #          angle = 45, hjust = -0.1,
  #          position = position_dodge(width = 0.6), size = 4) +
  coord_flip() +
  facet_grid(rows = vars(class), space = "free", scales = "free") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c( "#2ca02c","#9467bd","#ff7f0e")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3)) +
  #scale_alpha_discrete(values = c(0.8, 0.7)) +
  labs(y = "Odds Ratio (95% CI)",
       x = "Outcome") +
  theme_classic(base_size = 24) +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.box = "horizontal") +
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE), alpha = FALSE
  )


ggsave("plots/FX_diseases.png", p_res, width = 16, height = 9)

### BIPOLAR

data1$bp <- data1$p_1308920.0_prevalent
data1$ares_q <- ntile(abs(data1$res), 5)


m0 <- broom::tidy(glm("bp ~ res", data = data1 %>% filter(res < 0), family = binomial))
m1 <- broom::tidy(glm("bp ~ abs(res) + chrono", data = data1, family = binomial))
m2 <- broom::tidy(glm("bp ~ abs(res) + ever_insomnia", data = data1, family = binomial))
m3 <- broom::tidy(glm("bp ~ factor(ares_q) + chrono", data = data1, family = binomial))

data1 %>%
  filter(p_1308920.0_prevalent == 1) %>%
  ggplot(aes(x = chrono, y = ares_q, color = ever_insomnia)) +
  geom_jitter()

data1 %>%
  filter(!is.na(p_1308920.0_prevalent)) %>%
  ggplot(aes(x = res, color = p_1308920.0_prevalent)) +
  geom_density()

data1 %>%
  filter(!is.na(bp)) %>%
  ggplot(aes(x = chrono, y = abs(res), fill = bp)) +
  geom_boxplot() +
  theme_classic(base_size = 24) +
  labs(y = "Misalignment") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


