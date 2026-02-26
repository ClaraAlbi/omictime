
vars <- c(
  "age_recruitment", "sex", "p30079", "TDI", "bmi", "smoking", "season",
  "day_type", "fri_sun", "autumnDST", "springDST", "chrono",
  "h_sleep", "wakeup", "ever_insomnia", "shift_work", "night_shift",
  "employment", "chrono_Nightshift", "has_prescription",
  "antihypertensive", "sleep_medication", "antidepressants",
  "mood_stabiliser", "lithium"
)


cat_vars  <- vars[sapply(data[vars], is.factor)]
cont_vars <- vars[sapply(data[vars], is.numeric)]

comparisons <- list(
  acc_vs_mid = c("Accelerated", "Middle"),
  del_vs_mid = c("Delayed", "Middle")
)

cat_results <- map_dfr(names(comparisons), function(comp_name) {

  groups <- comparisons[[comp_name]]

  data_sub <- data %>%
    filter(res_sd_group %in% groups) %>%
    mutate(res_sd_group = droplevels(res_sd_group))

  group_counts <- table(data_sub$res_sd_group)

  # Skip entire comparison if small n
  if(any(group_counts < 10)) {
    return(NULL)
  }

  map_dfr(cat_vars, function(v) {

    tab <- table(data_sub[[v]], data_sub$res_sd_group)

    # Skip variable if any cell count < 10 per group total
    if(any(colSums(tab) < 10)) {
      return(NULL)
    }

    test <- chisq.test(tab)


      test_name <- "Chi-square"
      df_val <- unname(test$parameter)
      stat_val <- unname(test$statistic)


    data.frame(
      variable = v,
      comparison = comp_name,
      test_n = test_name,
      statistic = stat_val,
      df = df_val,
      p_value = test$p.value
    )
  })
})



cat_summary <- purrr::map_dfr(cat_vars, function(v) {

  data %>%
    count(predictor = v,
          level = .data[[v]],
          res_sd_group,
          name = "n") %>%
    group_by(predictor, res_sd_group) %>%
    mutate(percent = 100 * n / sum(n)) %>%
    ungroup() %>%
    mutate(
      value = sprintf("%d (%.1f%%)", n, percent),
      predictor_label = pretty_predictor[predictor],
      group = predictor_group[predictor]
    ) %>%
    select(group, predictor, predictor_label,
           level, res_sd_group, value)
})


cont_results <- map_dfr(names(comparisons), function(comp_name) {

  groups <- comparisons[[comp_name]]

  data_sub <- data %>%
    filter(res_sd_group %in% groups)

  group_counts <- table(data_sub$res_sd_group)

  map_dfr(cont_vars, function(v) {

    counts_per_group <- data_sub %>%
      group_by(res_sd_group) %>%
      summarise(n = sum(!is.na(.data[[v]])), .groups = "drop")

    if(any(counts_per_group$n < 10)) {
      return(NULL)
    }

    test <- wilcox.test(
      as.formula(paste(v, "~ res_sd_group")),
      data = data_sub
    )

    data.frame(
      variable = v,
      comparison = comp_name,
      test_n = "Wilcoxon",
      statistic = unname(test$statistic),
      p_value = test$p.value
    )
  })
})


cont_summary <- purrr::map_dfr(cont_vars, function(v) {

  data %>%
    group_by(res_sd_group) %>%
    summarise(
      predictor = v,
      mean = mean(.data[[v]], na.rm = TRUE),
      sd = sd(.data[[v]], na.rm = TRUE),
      n = sum(!is.na(.data[[v]])),
      .groups = "drop"
    ) %>%
    mutate(
      value = sprintf("%.2f (%.2f)", mean, sd),
      predictor_label = pretty_predictor[predictor],
      group = predictor_group[predictor],
      level = NA_character_
    ) %>%
    select(group, predictor, predictor_label,
           level, res_sd_group, value)
})

desc_table <- bind_rows(cat_summary, cont_summary)

desc_wide <- desc_table %>%
  pivot_wider(
    names_from = res_sd_group,
    values_from = value
  ) %>%
  filter(predictor != "chrono_Nightshift")




assoc_results <- bind_rows(cat_results, cont_results) %>%
  rename(predictor = variable,
         p = p_value) %>%
  mutate(
    group = unname(predictor_group[predictor]),
    predictor_label = unname(pretty_predictor[predictor]),
    comparison = recode(comparison,
                        acc_vs_mid = "Accelerated vs Middle",
                        del_vs_mid = "Delayed vs Middle")
  ) %>%
  filter(predictor != "chrono_Nightshift") %>%
  mutate(p_adj = p.adjust(p))

test_wide <- assoc_results %>%
  select(predictor, comparison, p_adj, test_n) %>%
  mutate(test_label = paste0(test_n, ", p=", signif(p_adj, 3))) %>%
  select(-test_n, -p_adj) %>%
  pivot_wider(
    names_from = comparison,
    values_from = test_label
  )

final_table <- desc_wide %>%
  left_join(test_wide, by = "predictor") %>%
  arrange(
    factor(group,
           levels = c("Demographics",
                      "Seasonality",
                      "Sleep questionnaire",
                      "Employment",
                      "Medication"))
  )



f_t <- final_table %>%
  filter(!is.na(level)) %>%
  select(Data = group, Predictor = predictor_label, Category = level, everything())


writexl::write_xlsx(f_t, "tables/extremes_associations.xlsx")
saveRDS(f_t, "tables/extremes_associations.rds")

