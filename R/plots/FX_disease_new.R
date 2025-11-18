library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(ggplot2)
library(stringr)
install.packages("broom")
library(ggh4x)

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
    paste0("`", v, "` ~ abs(res) + ",
           paste(base_covars, collapse = " + ")))
  f_prev3 <- as.formula(
    paste0("`", v, "` ~ abs(res) + ",
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
          left_join(d_counts %>% rename(outcome = name)), "data_share/association_results_disease_CM.rds")



###

library(tidyverse)
library(forcats)

results <- readRDS("data_share/association_results_disease_CA.rds") %>%
  bind_rows(readRDS("data_share/association_results_disease_CM.rds"))

fields <- data.table::fread("data/field.tsv")

r <- results %>%
  filter(term %in% c("res", "abs(res)")) %>%
  mutate(expo = case_when(term == "res" ~ "Circadian Acceleration",
                          term == "abs(res)" ~ "Circadian Misalignment"),
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
  distinct(field_id, model, expo, .keep_all = TRUE)

r %>%
  distinct(field_id, title, family, disorder, n) %>%
  arrange(desc(n))

p_res <-
  r %>%
  mutate(model = forcats::fct_rev(factor(model)),
         disorder = paste0(disorder, "\n", n)) %>%
  #slice(1:80) %>%
  ggplot(aes(
    x = disorder,
    y = estimate,
    ymin = estimate - std.error, ymax = estimate + std.error,
    color = model,
    alpha = p.adjust(p.value) < 0.05
  )) +
  # points with vertical error bars (dodge by model)
  geom_pointrange(position = position_dodge(width = 0.8), size = 1, fatten = 1.5) +
  coord_flip() +
  facet_nested(
    cols = vars(expo),
    rows = vars(family),
    scales = "free_y",
    space = "free_y"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_color_manual(values = c("#2ca02c", "#9467bd", "#ff7f0e")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  labs(y = "Odds ratio (95% CI)",
       x = NULL) +
  theme_classic(base_size = 14) +
  theme(
    # place legend inside plot at top-right
    legend.position = c(0.6, 1),
    legend.justification = c("right", "top"),
    strip.background = element_rect(fill = "antiquewhite2", color = "black", linewidth = 0.8),
    legend.title = element_blank(),
    legend.box = "vertical",
    panel.spacing = unit(0.75, "lines"),
    strip.text.y.left = element_text(angle = 0, face = "bold", vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE, reverse = TRUE))


ggsave("plots/F7_diseases_log.png", p_res, width = 10, height = 11)




summarise_covs <- function(disease_col) {
  df <- data1 %>% filter(.data[[disease_col]] == 1)

  n_total <- nrow(df)

  tibble(
    disease = disease_col,
    N = n_total,
    age_mean = mean(df$age_recruitment, na.rm = TRUE),
    age_sd   = sd(df$age_recruitment, na.rm = TRUE),
    age_median = median(df$age_recruitment, na.rm = TRUE),
    age_q1   = quantile(df$age_recruitment, 0.25, na.rm = TRUE),
    age_q3   = quantile(df$age_recruitment, 0.75, na.rm = TRUE),

    sex_0 = sum(df$sex == 0, na.rm = TRUE),
    sex_1 = sum(df$sex == 1, na.rm = TRUE),

    bmi_mean = mean(df$bmi, na.rm = TRUE),
    bmi_sd   = sd(df$bmi, na.rm = TRUE),
    bmi_median = median(df$bmi, na.rm = TRUE),
    bmi_q1   = quantile(df$bmi, 0.25, na.rm = TRUE),
    bmi_q3   = quantile(df$bmi, 0.75, na.rm = TRUE),
    bmi_missing = sum(is.na(df$bmi)),

    smoke_0 = sum(df$smoking == 0, na.rm = TRUE),
    smoke_1 = sum(df$smoking == 1, na.rm = TRUE),
    smoke_2 = sum(df$smoking == 2, na.rm = TRUE)
  )
}

disease_cols <- colnames(data1)[17:60]

# Covariates to summarise
covs <- c("age_recruitment", "sex", "bmi", "smoking")
cov_stats_tbl <- map_dfr(disease_cols, summarise_covs)

saveRDS(cov_stats_tbl, "data_share/stats_prevalent.rds")





