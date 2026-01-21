library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)




phy <- data.table::fread("/mnt/project/physical.csv") %>%
  rowwise() %>%
  mutate(sbp = mean(c(p4080_i0_a0, p4080_i0_a1)),
         dbp = mean(c(p4079_i0_a0, p4079_i0_a1))) %>%
  ungroup()  %>% select(eid, sbp, dbp)

mh <- data.table::fread("/mnt/project/psychosocial_MH.csv") %>%
  select(eid, rec_dep = p2050_i0,
         rec_rest = p2070_i0,
         rec_tired = p2080_i0,
         rec_enth = p2060_i0) %>%
  filter(rec_rest != "" | rec_tired != "" | rec_enth != "") %>%
  mutate(rec_rest = factor(rec_rest, levels = c("Not at all", "Several days", "More than half the days", "Nearly every day")),
         rec_dep = factor(rec_dep, levels = c("Not at all", "Several days", "More than half the days", "Nearly every day")),
         rec_tired =factor(rec_tired, levels = c("Not at all", "Several days", "More than half the days", "Nearly every day")),
        rec_enth = factor(rec_enth, levels = c("Not at all", "Several days", "More than half the days", "Nearly every day")))


data <- data %>%
  mutate(res_extreme = as.integer(res > 2 | res < -2))

results <- map_dfr(vars, function(v) {
  adj_vars <- if (v %in% c("sex", "age_recruitment")) paste0("PC", 1:20) else covars

  rhs <- paste(c(v, adj_vars), collapse = " + ")
  f <- as.formula(paste("res_extreme ~", rhs))

  fit <- glm(f, data = data, family = binomial)

  tidy(fit) %>%
    filter(str_detect(term, paste0("^", v))) %>%
    mutate(
      predictor = v,
      reference = FALSE,
      odds_ratio = exp(estimate)
    ) %>%
    {
      if (is.factor(data[[v]])) {
        ref <- tibble(
          term = paste0(v, levels(data[[v]])[1]),
          estimate = 0,
          std.error = NA,
          statistic = NA,
          p.value = NA,
          odds_ratio = 1,
          predictor = v,
          reference = TRUE
        )
        bind_rows(ref, .)
      } else .
    }
})


data <- data %>%
  mutate(res_gt2 = as.integer(res > 2))

results_gt2 <- map_dfr(vars, function(v) {
  adj_vars <- if (v %in% c("sex", "age_recruitment")) paste0("PC", 1:20) else covars

  rhs <- paste(c(v, adj_vars), collapse = " + ")
  f <- as.formula(paste("res_gt2 ~", rhs))

  fit <- glm(f, data = data, family = binomial)

  tidy(fit) %>%
    filter(str_detect(term, paste0("^", v))) %>%
    mutate(
      predictor = v,
      tail = "res > 2",
      odds_ratio = exp(estimate),
      reference = FALSE
    )
})


data <- data %>%
  mutate(res_lt2 = as.integer(res < -2))

results_lt2 <- map_dfr(vars, function(v) {
  adj_vars <- if (v %in% c("sex", "age_recruitment", "chrono")) paste0("PC", 1:20) else covars

  rhs <- paste(c(v, adj_vars), collapse = " + ")
  f <- as.formula(paste("res_lt2 ~", rhs))

  fit <- glm(f, data = data, family = binomial)

  tidy(fit) %>%
    filter(str_detect(term, paste0("^", v))) %>%
    mutate(
      predictor = v,
      tail = "res < -2",
      odds_ratio = exp(estimate),
      reference = FALSE
    )
})

results <- bind_rows(results_gt2, results_lt2)
