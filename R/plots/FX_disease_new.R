
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)
library(stringr)
install.packages("broom")

library(table1)

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = as.factor(sex))

job_vars <- data.table::fread("/mnt/project/job_vars.tsv") %>%
  mutate(night_shift = case_when(`3426-0.0` == 1 ~ "Never",
                                 `3426-0.0` == 2 ~ "Sometimes",
                                 `3426-0.0` == 3 ~ "Usually",
                                 `3426-0.0` == 4 ~ "Always")) %>%
  filter(`3426-0.0` %in% 1:4) %>%
  mutate(night_shift = as.factor(night_shift),
         night_shift = relevel(night_shift, ref = "Never"))

sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid,
         h_sleep = `1160-0.0`,
         chrono = `1180-0.0`,
         ever_insomnia = `1200-0.0`,
         wakeup = `1170-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ "Definitely morning",
    chrono == 2 ~ "Rather morning",
    chrono == -1~ "Don't know",
    chrono == 3 ~ "Rather evening",
    chrono == 4 ~ "Definitely evening",
    TRUE ~ NA_character_),
    chrono = factor(chrono, levels = c("Definitely morning", "Rather morning", "Don't know", "Rather evening", "Definitely evening")),
    ever_insomnia = case_when(ever_insomnia == 1 ~ "Never/rarely",
                              ever_insomnia == 2 ~ "Sometimes",
                              ever_insomnia == 3 ~ "Usually"))

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  #filter(time_day > 12 & time_day < 18) %>%
  filter(i == 0) %>%
  #filter(cv %in% 1:5) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  ungroup() %>%
  mutate(season = case_when(m %in% c("12", "01", "02") ~ "Winter",
                            m %in% c("03", "04", "05") ~ "Spring",
                            m %in% c("06", "07", "08") ~ "Summer",
                            m %in% c("09", "10", "11") ~ "Fall"),
         season = relevel(as.factor(season), ref = "Winter"))

data <- df %>%
  left_join(covs) %>%
  left_join(job_vars) %>%
  left_join(sleep) %>%
  filter(!is.na(chrono)) %>%
  mutate(h = round(time_day, 0)) %>%
  filter(age_recruitment > 39) %>%
  filter(h_sleep > 0)

data$res <- residuals(lm(pred_mean ~ time_day, data = data))
data$res_q <- ntile(data$res, 5)


# diagnosis table
sleep2 <- readRDS("/mnt/project/diseases_circadian.rds") %>% select(eid, p131060, p130708, p130792, p130842, p131306)

outcomes <- tribble(~field_id, ~phen, ~class,
                    "130894", "Depressive_episode","Neuro-psychiatric",
                    "130896", "Recurrent_depression","Neuro-psychiatric",
                    "130892", "Bipolar disorder","Neuro-psychiatric",
                    "130874", "Schizophrenia","Neuro-psychiatric",
                    "p131060", "Sleep - G47", "Neuro-psychiatric",
                    "p130842", "Dementia", "Neuro-psychiatric",
                    "p130708", "Type 2 Diabetes", "Cardiometabolic",
                    "p130792", "Obesity", "Cardiometabolic",
                    "p131306", "Ischaemic heart disease", "Cardiometabolic")

dis2 <- data.table::fread("/mnt/project/vars_diseases_2.tsv") %>%
  select(eid, contains(outcomes$field_id)) %>%
  inner_join(sleep2) %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds") %>% select(eid, date_bsampling)) %>%
  mutate(across(contains("13"), ~ case_when(.x < date_bsampling ~ 1,
                                               is.na(.x) ~ 0), .names = "{.col}_prevalent")) %>%
  select(-contains("_prevalent_incident")) %>%
  select(eid, contains("prevalent")) %>%
  rename_with(~ paste0("p_", .x), -eid) %>%
  rename_with(~ str_remove( .x, "-"), -eid) %>%
  mutate(across(-eid, as.factor))

data1 <- data %>%
  left_join(dis2)


pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

tab_desc <- table1::table1(~ abs(res) + time_day + age_recruitment + factor(sex) + chrono + h_sleep + ever_insomnia | p_1308920.0_prevalent,
                           data = data1 %>% filter(!is.na(p_1308920.0_prevalent)), overall=F,
                           render.cont = my_render_cont, extra.col = list(`P-value`=pvalue))


data1 %>%
  filter(!is.na(p_1308920.0_prevalent) & !is.na(ever_insomnia)) %>%
  ggplot(aes(x = chrono, y = abs(res), fill = p_1308920.0_prevalent)) +
  geom_boxplot() +
  theme_classic() +
  labs(y = "Circadian misalignment", fill = "Bipolar") +
  theme(text = element_text(size = 20))

data1 %>%
  filter(!is.na(p_1308920.0_prevalent) & !is.na(ever_insomnia)) %>%
  ggplot(aes(x = ever_insomnia, y = abs(res), fill = p_1308920.0_prevalent)) +
  geom_boxplot() +
  facet_grid(~chrono) +
  theme_classic() +
  labs(y = "Circadian misalignment", fill = "Bipolar") +
  theme(text = element_text(size = 30),
        axis.text.x = element_text(angle = 45, vjust = 0.5))

vars <- colnames(dis2)[-1]

base_covars   <- c("sex","age_recruitment")
extra_covars  <- c("bmi", "smoking")

results <- map_dfr(vars, function(v) {

  # formula for abs(res) only
  f0 <- as.formula(paste0("`", v, "` ~ time_day"))
  f_prev1 <- as.formula(
    paste0("`", v, "` ~ abs(res)")
  )

  # formula for abs(res) + covariates
  f_prev2 <- as.formula(
    paste0("`", v, "` ~ abs(res) + ",
           paste(c(base_covars, extra_covars), collapse = " + "))
  )

  # fit models
  mt <- glm(f0, data = data1, family = binomial)
  m0 <- glm(f_prev1, data = data1, family = binomial)
  m1 <- glm(f_prev2, data = data1, family = binomial)

  # collect results
  bind_rows(
    broom::tidy(mt, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(model = "mt_TIME_DAY", outcome = v),

    broom::tidy(m0, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(model = "m0_MISALIGNMENT", outcome = v),

    broom::tidy(m1, conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(model = "m1_MISALIGNMENT + COV", outcome = v)
  )
})

results
library(forcats)
results %>%
  mutate(
    # Extract the core field id
    o = str_extract(outcome, "(?<=p_)[0-9]+(?=0\\.0_prevalent)|(?<=p_)p?[0-9]+"),
    model = factor(model, levels =c("mt_TIME_DAY", "m0_MISALIGNMENT", "m1_MISALIGNMENT + COV"))
  ) %>%
  left_join(outcomes, by = c("o" = "field_id")) %>%
  mutate(outcome_clean = coalesce(phen, outcome)) %>%
  filter(term != "(Intercept)") %>%
  filter(term %in% c("time_day", "abs(res)")) %>%
  ggplot(aes(x = outcome_clean,
             y = estimate,
             ymin = conf.low, ymax = conf.high,
             color = model, shape = model)) +
  geom_pointrange(position = position_dodge(width = 0.6),
                  size = 1, fatten = 3) +
  geom_text(data = . %>% filter(p.value < 0.05), aes(label = sprintf("%.0e", p.value)),
            angle = 45, hjust = -0.1,
            position = position_dodge(width = 0.6), size = 4) +
  coord_flip() +
  facet_grid(rows = vars(class), space = "free", scales = "free") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c( "#2ca02c","#9467bd","#ff7f0e")) +
  labs(y = "Odds Ratio (95% CI)",
       x = "Outcome") +
  theme_classic(base_size = 16) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.box = "horizontal")
  guides(
    color = guide_legend(nrow = 3, byrow = TRUE),
    shape = guide_legend(nrow = 3, byrow = TRUE)
  )



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


