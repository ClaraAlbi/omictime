install.packages("broom")
library(stringr)

### Symptoms

covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  filter(smoking != "-3") %>%
  mutate(bmi = weight/(height/100)^2,
         sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")) ,
         smoking = factor(smoking, levels = c(0,1,2), labels = c("Never", "Previous", "Current")),
  )

pcs <- data.table::fread("/mnt/project/covariates.txt") %>%
  select(eid = 1, contains("PC"))

sleep_mh <- data.table::fread("/mnt/project/physical.csv") %>% select(eid, p20532,p20534,p20533,p20535)

sleep_mh_t <- sleep_mh %>%
  #filter(eid %in% data$eid) %>%
  count(p20532,p20534,p20533,p20535)
  pivot_longer(c(p20534,p20533,p20535)) %>%
  count(p20532,name, value)

phens <- sleep_mh %>%
  #filter(p20532 %in% c("Yes", "No")) %>%
  unite("sleep_depression", c(p20532, p20534,p20533,p20535), sep = "_") %>%
  mutate(sleep_depression = case_when(sleep_depression == "No___" ~ "No",
                                      #sleep_depression == "Do not know___" ~ "Do not know",
                                      sleep_depression == "Yes_No_No_Yes" ~ "Waking too early (WTE)",
                                      sleep_depression == "Yes_No_Yes_No" ~ "Trouble falling asleep (TFA)",
                                      sleep_depression == "Yes_No_Yes_Yes" ~ "WE + TFA",
                                      sleep_depression == "Yes_Yes_No_No" ~ "Sleep too much (STM)",
                                      sleep_depression == "Yes_Yes_No_Yes" ~ "STM + WE",
                                      sleep_depression == "Yes_Yes_Yes_No" ~ "STM + TFA",
                                      sleep_depression == "Yes_Yes_Yes_Yes" ~ "ALL"),
         sleep_depression = relevel(factor(sleep_depression), ref = "No"))


data <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  inner_join(readRDS("/mnt/project/olink_int_replication.rds") %>% select(-date_bsampling)) %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  left_join(covs) %>%
  left_join(pcs) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  left_join(phens)


data$res <- residuals(lm(pred_mean ~ time_day, data = data))

vars <- "sleep_depression"

covars <- c("sex", "age_recruitment", paste0("PC", 1:10))

results <- map_dfr(vars, function(v) {
  adj_vars <- if (v %in% c("sex", "age_recruitment")) paste0("PC", 1:10) else covars

  # Combine predictor + covariates safely
  rhs <- paste(c(v, adj_vars), collapse = " + ")
  f <- as.formula(paste("abs(res) ~ ", rhs))

  fit <- lm(f, data = data)

  broom::tidy(fit) %>%
    filter(str_detect(term, paste0("^", v))) %>%
    mutate(predictor = v, reference = FALSE) %>%
    {
      if (is.factor(data[[v]])) {
        ref <- tibble(
          term = paste0(v, levels(data[[v]])[1]),
          estimate = 0, std.error = NA, statistic = NA, p.value = NA,
          predictor = v, reference = TRUE
        )
        bind_rows(ref, .)
      } else .
    }
})

saveRDS(results, "data_share/results_associations_sleep_depression_CM.rds")


library(table1)
my_render_cont <- function(x){
  with(
    stats.apply.rounding(stats.default(x)),
    c(
      "",
      `Mean (SD)` = sprintf("%s (%s)", MEAN, SD),
      `Median [Q1, Q3]` = sprintf("%s [%s, %s]",
                                  MEDIAN, Q1, Q3)
    )
  )
}
tab_desc <- table1::table1(~ age_recruitment + sex + sleep_depression,
                           data = data,
                           render.cont = my_render_cont)

