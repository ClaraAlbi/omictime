library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(lubridate)
#install.packages("survminer")
library(survminer)
library(ggpubr)
library(ggplot2)
library(stringr)

# administrative censor date
censor_date <- as.Date("2024-09-01")

# base cohort: death + baseline date
base_cohort <- fread("/mnt/project/clara/death.csv") %>%
  transmute(
    eid,
    death_date    = as.Date(p40000_i0)
  ) %>%
  left_join(
    readRDS("/mnt/project/biomarkers/time.rds") %>%
      mutate(date_bsampling = as.Date(date_bsampling)) %>%
      filter(date_bsampling > "1900-01-01") %>%
      select(eid, date_bsampling),
    by = "eid"
  )

# biomarker residuals + covs (sex, age, PCs, BMI, smoking, fasting)

bio_covs <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  separate(date_bsampling, into = c("y", "m", "d"), sep = "-", remove = T) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  select(eid, time_day, pred_mean) %>%
  left_join(covs) %>%
  ungroup() %>%

  mutate(
    res        = residuals(lm(pred_mean ~ time_day, data = cur_data())),
    res_q     = factor(ntile(res, 5), levels = 1:5),
    res_abs    = abs(res),
    ares_q     = factor(ntile(res_abs, 5), levels = 1:5),
  ) %>%
  # ensure covariates are correct type
  mutate(
    sex           = factor(sex),
    smoking       = factor(smoking),
    age_recruitment = age_recruitment,
    BMI           = weight/((height/100)^2)
  )


# diagnosis table

dis2 <- data.table::fread("/mnt/project/vars_diseases_2.tsv") %>%
  inner_join(readRDS("/mnt/project/diseases_circadian.rds") %>% select(-p130896, -p130894, -p130892, -p130842, -p130836, -p130906))

colnames(dis2) <- c("eid", paste0("p",str_remove(str_remove(colnames(dis2)[-1], "-0.0"), "p")))

dis2_inc <-  dis2 %>%
  left_join(readRDS("/mnt/project/biomarkers/time.rds") %>% select(eid, date_bsampling)) %>%
  mutate(across(-c(eid, date_bsampling), ~ case_when(.x > (ymd(date_bsampling) + 365.25) ~ 1,
                                            is.na(.x) ~ 0), .names = "{.col}_incident")) %>%
  select(eid, contains("incident"))

fields <- data.table::fread("/mnt/project/Showcase metadata/field.tsv")

d_counts <-
  dis2_inc %>%
  pivot_longer(-eid, names_to = "outcome") %>%
  group_by(outcome) %>% count(value) %>%
  filter(value == 1 & n > 40) %>%
  mutate(outcome = str_remove(outcome, "p"),
         field_id = as.numeric(str_remove(outcome, "_incident"))) %>%
  left_join(fields %>% select(field_id, title)) %>%
  mutate(family = str_extract(title, "(?<=Date )\\S+(?= first reported)"),
         family = str_sub(family, 1, 2),
         disorder = sub(".*\\((.*)\\).*", "\\1", title),
         disorder = str_to_sentence(disorder)) %>%
  filter(!field_id  %in% c(130898, 130902, 130932, 130944, 130852, 130848)) %>%
  filter(!family %in% c("F4", "F5", "F6", "F1" )) %>%
  distinct(field_id, .keep_all = TRUE)



### COX

disease_field <- colnames(dis2)[which(colnames(dis2) %in% paste0("p" ,d_counts$field_id))]

base_covars   <- c("sex","age_recruitment")
extra_covars <- c("smoking", "bmi")


results <- map_dfr(disease_field, function(v) {

  df_plot <- base_cohort %>%
    left_join(dis2 %>% select(eid, !!sym(v))
              %>% rename(diag_date =  !!sym(v)),
              by = "eid") %>%
    left_join(bio_covs %>% select(eid, res, sex, age_recruitment, smoking, bmi), by = "eid") %>%
    filter(is.na(death_date) | death_date > date_bsampling,
           is.na(diag_date)  | diag_date  > date_bsampling) %>%
    mutate(
      end_date = as.Date(pmin(diag_date, death_date, censor_date, na.rm = TRUE)),
      event    = as.integer(!is.na(diag_date) &
                              diag_date <= pmin(death_date, censor_date, na.rm = TRUE)),
      time_yrs = as.numeric(end_date - date_bsampling) / 365.25
    ) %>%
    filter(time_yrs > 1)

  # formula for abs(res) only

  f_cox0 <- as.formula("Surv(time_yrs, event) ~ abs(res)")
  f_cox1 <- as.formula(paste("Surv(time_yrs, event) ~ abs(res)", " + ",  paste(base_covars, collapse = " + ")))
  f_cox2 <- as.formula(paste("Surv(time_yrs, event) ~ abs(res)", " + ",  paste(c(base_covars, extra_covars), collapse = " + ")))



  # fit models
  m1 <- survival::coxph(f_cox0, data = df_plot, x = FALSE, y = FALSE)
  m2 <- survival::coxph(f_cox1, data = df_plot, x = FALSE, y = FALSE)
  m3 <- survival::coxph(f_cox2, data = df_plot, x = FALSE, y = FALSE)


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
          left_join(d_counts %>% mutate(outcome = paste0("p", field_id))), "data_share/association_results_cox_disease_CM.rds")



# 2. Function to build the dataset and return a survminer plot
make_ci_plot <- function(disease_field) {
  d <- d_counts$disorder[which(disease_field == paste0("p",d_counts$field_id))]

  df_plot <- base_cohort %>%
    left_join(dis2 %>% select(eid, !!sym(disease_field))
              %>% rename(diag_date =  !!sym(disease_field)),
              by = "eid") %>%
    left_join(bio_covs %>% select(eid, res_q, ares_q), by = "eid") %>%
    filter(is.na(death_date) | death_date > date_bsampling,
           is.na(diag_date)  | diag_date  > date_bsampling) %>%
    mutate(
      end_date = as.Date(pmin(diag_date, death_date, censor_date, na.rm = TRUE)),
      event    = as.integer(!is.na(diag_date) &
                              diag_date <= pmin(death_date, censor_date, na.rm = TRUE)),
      time_yrs = as.numeric(end_date - date_bsampling) / 365.25
    ) %>%
    filter(time_yrs > 1, ares_q %in% c("1","3","5"))

  # Fit survival curves
  fit <- survfit(Surv(time_yrs, event) ~ ares_q, data = df_plot)
  test <- survdiff(Surv(time_yrs, event) ~ ares_q, data = df_plot)
  pval <- 1 - pchisq(test$chisq, length(test$n) - 1)

  # Plot cumulative incidence = 1 – survival

  df_km_full <- survminer::surv_summary(fit, data = df_plot)

  # thin points by taking every nth row per strata (reduces plotting clutter)
  df_km_thin <- df_km_full %>%
    group_by(strata) %>%
    mutate(idx = row_number()) %>%
    filter(idx %% 4 == 1 | idx == max(idx)) %>%  # keep every 4th point + last
    ungroup()

  # ribbon: use the full summary to draw continuous ribbons (geom_ribbon doesn't need every point)
  df_ribbon <- df_km_full %>%
    group_by(strata) %>%
    arrange(time) %>%
    ungroup()

  ggplot() +
    # ribbon for CI (cumulative incidence = 1 - surv)
    geom_ribbon(data = df_ribbon,
                aes(x = time, ymin = 1 - upper, ymax = 1 - lower, fill = strata, group = strata),
                alpha = 0.12, colour = NA) +
    # step curve for cumulative incidence
    geom_step(data = df_km_full, aes(x = time, y = 1 - surv, color = strata, group = strata), size = 0.9) +
    # thinned points for emphasis
    geom_point(data = df_km_thin, aes(x = time, y = 1 - surv, color = strata), size = 0.8) +
    #scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
    scale_color_manual(values = c("steelblue", "forestgreen", "firebrick"), ) +
    scale_fill_manual(values = c("steelblue", "forestgreen", "firebrick")) +
    labs(
      title =  str_sub(d, 1, 10),
      subtitle = paste0("pval: ",round(pval, 3)),
      x = "y",
      y = "Cumulative incidence",
      color = "res_abs quintile",
      fill = NULL
    ) +
    coord_cartesian(xlim = c(0, max(df_km_full$time, na.rm = TRUE))) +
    theme_minimal(base_size = 8) +
    guides(fill = "none") +
    theme(
      #legend.position = "none",
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = alpha("white", 0.8), colour = NA),
      legend.title = element_text(size = 9),
      panel.grid.minor = element_blank(),
      axis.title.y = element_blank()
    )

}


selected_p <- colnames(dis2)[which(colnames(dis2) %in% paste0("p" ,d_counts$field_id))]

# 3. Generate one plot per disease
ci_plots <- map(selected_p, make_ci_plot)

install.packages("patchwork")
library(patchwork)

combined <- wrap_plots(ci_plots, ncol = 5, nrow = 6)
combined
# save
ggsave("plots/CI_ares_q.png", combined, width = 7, height = 10, dpi = 300)


legend <- get_legend(p)
library(grid)
grid.newpage()
grid.draw(legend)

# 3. Save to file
ggsave("plots/legend.png", legend, width = 3, height = 1.5, dpi = 300)




### PLOTS FINAL


results <- readRDS("data_share/association_results_cox_disease_CA.rds") %>%
  bind_rows(readRDS("data_share/association_results_cox_disease_CM.rds"))

fields <- data.table::fread("data/field.tsv")

r <- results %>%
  filter(!field_id %in% c(130888)) %>%
  filter(term %in% c("res", "abs(res)")) %>%
  mutate(expo = case_when(term == "res" ~ "Circadian Acceleration",
                          term == "abs(res)" ~ "Circadian Misalignment"),
         disorder = case_when(disorder == "" ~ "Delirium",
                              TRUE ~ disorder))


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
  labs(y = "HR (95% CI)",
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


ggsave("plots/F7_diseases_cox.png", p_res, width = 10, height = 11)





data1 <- bio_covs %>%
  left_join(dis2_inc)

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

disease_cols <- colnames(dis2_inc)[-1]

# Covariates to summarise
covs <- c("age_recruitment", "sex", "bmi", "smoking")
cov_stats_tbl <- map_dfr(disease_cols, summarise_covs)

saveRDS(cov_stats_tbl, "data_share/stats_incident.rds")




