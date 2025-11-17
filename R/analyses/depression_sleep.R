install.packages("broom")
library(stringr)

### Symptoms

mh <- data.table::fread("/mnt/project/psychosocial_MH.csv")

sleep_mh <- data.table::fread("/mnt/project/physical.csv") %>% select(eid, p20532,p20534,p20533,p20535)

sleep_mh_t <- sleep_mh %>%
  filter(eid %in% data$eid) %>%
  count(p20532,p20534,p20533,p20535)
  pivot_longer(c(p20534,p20533,p20535)) %>%
  count(p20532,name, value)

phens <- sleep_mh %>%
  filter(p20532 %in% c("Yes", "No")) %>%
  unite("sleep_depression", c(p20532, p20534,p20533,p20535), sep = "_") %>%
  mutate(sleep_depression = case_when(sleep_depression == "No___" ~ "No",
                                      sleep_depression == "Yes_No_No_Yes" ~ "Waking too early",
                                      sleep_depression == "Yes_No_Yes_No" ~ "Trouble falling asleep",
                                      sleep_depression == "Yes_No_Yes_No" ~ "WE + TFA",
                                      sleep_depression == "Yes_Yes_No_No" ~ "Sleep too much",
                                      sleep_depression == "Yes_Yes_No_Yes" ~ "STM + WE",
                                      sleep_depression == "Yes_Yes_Yes_No" ~ "STM + TFA",
                                      sleep_depression == "Yes_Yes_Yes_Yes" ~ "ALL"),
         sleep_depression = relevel(factor(sleep_depression), ref = "No"))

  unite("trouble_falling_asleep", c(p20532,p20533), sep = "_", remove = F) %>%
  unite("waking_too_early", c(p20532,p20535), sep = "_", remove = F) %>%
  filter(p20532 %in% c("Yes", "No"))

data <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  inner_join(readRDS("/mnt/project/olink_int_replication.rds") %>% select(-date_bsampling)) %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  left_join(phens)


  mutate(across(c(p20534,p20533,p20535), ~factor(.x, levels = c("No", "Yes"))))

  left_join(mh) %>%
  mutate(across(p1920_i0:p2040_i0, ~factor(.x, levels = c("No", "Yes", "Do not know", "Prefer not to answer"))))
         across(
           where(~ !is.numeric(.x)),
           ~ case_when(
             .x == "" ~ NA_character_,
             TRUE ~ .x
           )
         ))

data$res <- residuals(lm(pred_mean ~ time_day, data = data))

vars <- colnames(mh)[c(2:14)]

vars <- c("p20534", "p20533", "p20535")

vars <- "sleep_depression"

results <- map_dfr(vars, function(v) {
  adj_vars <- ""

  # Combine predictor + covariates safely
  rhs <- paste(c(v, adj_vars), collapse = " + ")
  f <- as.formula(paste("abs(res) ~", v))

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
