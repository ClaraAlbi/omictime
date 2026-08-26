source("config/paths.R")

library(data.table)

res_dir <- bmrc_path("results/cojo")

cojo_files <- list.files(
  res_dir,
  pattern = "\\.jma\\.cojo$",
  full.names = TRUE
)


traits = c(
  "ACC_M10_TIME",
  "ACC_L5_TIME",
  "ACC_SLEEP_MIDP",
  "ACC_DIURNAL_INACT",
  "ACC_N_SLEEP_EPISODES",
  "ACC_SLEEP_DUR",
  "ACC_SLEEP_DUR_SD",
  "ACC_SLEEP_EFF",
  "chrono_v2",
  #"morning_person-5.1",
  #"chrono_b",
  #"RA_logistic",
  #"RA_linear"
)

cojo_files_matched <- cojo_files[
  sapply(traits, function(t)
    grepl(t, basename(cojo_files))
  ) |> apply(1, any)
]

cojo_list <- lapply(cojo_files_matched, fread)
names(cojo_list) <- sub("\\.jma\\.cojo$", "", basename(cojo_files_matched))

cojo_tbl <- cojo_list %>%
  imap_dfr(~ as_tibble(.x) %>% mutate(source = .y)) %>%
  mutate(source = str_remove(source, "_?chr[0-9XY]+$")) %>%
  group_by(source) %>% count()


saveRDS(cojo_tbl, "cojo/cojo_results_all.rds")


####

gwas_CA <- fread("data/GWAS_reslasso/res_lasso_GWAS_merged.fastGWA") %>%
  rename(
    beta_CA = BETA,
    se_CA   = SE,
    AF1_CA  = AF1,
    P_CA    = P,
    N_CA    = N
  )

cojo_by_exposure <- split(cojo_tbl, cojo_tbl$source)

mr_results <- list()

for (exposure_name in names(cojo_by_exposure)[7]) {
  
  res <- run_mr_one_exposure_from_tbl(
    exposure_tbl = cojo_by_exposure[[exposure_name]],
    exposure_name = exposure_name,
    outcome_gwas = gwas_CA,
    outcome_name = "CA",
    outcome_beta = "beta_CA",
    outcome_se = "se_CA",
    outcome_eaf = "AF1_CA",
    outcome_p = "P_CA",
    outcome_n = "N_CA"
  )
  
  if (!is.null(res)) {
    mr_results[[exposure_name]] <- res
  }
}


mr_CHRONO <- list()
for (exposure_name in names(cojo_by_exposure)[9]) {
  
  res <- run_mr_one_exposure_from_tbl(
    exposure_tbl = cojo_by_exposure[[exposure_name]],
    exposure_name = exposure_name,
    outcome_gwas = gwas_CA,
    outcome_name = "CA",
    outcome_beta = "beta_CA",
    outcome_se = "se_CA",
    outcome_eaf = "AF1_CA",
    outcome_p = "P_CA",
    outcome_n = "N_CA"
  )
  
  if (!is.null(res)) {
    mr_CHRONO[[exposure_name]] <- res
  }
}



mr_table_ca <- bind_rows(
  lapply(mr_results, function(x) x$mr)
) %>%
  select(-id.exposure, -id.outcome) %>%
  mutate(
    pval_FDR = p.adjust(pval, method = "BH"),
    method = case_when(
      method == "Inverse variance weighted (multiplicative random effects)" ~ "IVW",
      method == "Weighted median" ~ "Weighted median",
      method == "MR Egger" ~ "Egger",
      TRUE ~ method
    )
  )

egger_table <- map_dfr(
  mr_results,
  ~ if (!is.null(.x$egger_intercept)) .x$egger_intercept
) %>%
  transmute(
    exposure,
    outcome,
    egger_intercept = egger_intercept,
    egger_se = se,
    egger_pval = pval,
    egger_pval_FDR = p.adjust(pval, method = "BH")
  )

presso_table <- map_dfr(mr_results, extract_presso)

final_table_ca <- mr_table_ca %>%
  left_join(egger_table, by = c("exposure", "outcome")) %>%
  left_join(presso_table, by = c("exposure", "outcome"))


harm_table <- map_dfr(
  mr_results,
  function(x) {
    tryCatch({
      
      harm <- x$harm_snps
      if (is.null(harm) || nrow(harm) == 0) return(NULL)
      
      harm %>%
        mutate(
          exposure = as.character(x$mr$exposure[1]),
          outcome = as.character(x$mr$outcome[1])
        )
      
    }, error = function(e) NULL)
  }
)


p_mr1 <- harm_table %>%
  filter(!exposure %in% c("chrono_b"))
  left_join(mr_table_ca, by = "exposure") %>%
  filter(method == "IVW") %>%
  mutate(text = paste0("n=", nsnp, ", bIVW=", round(b, 2), ", pFDR=", formatC(pval_FDR, format = "e", digits = 0))) %>%
  ggplot(aes(x = beta.exposure, y = beta.outcome, color = mr_keep,
             ymin = beta.outcome - se.outcome,
             ymax = beta.outcome + se.outcome,
             xmin = beta.exposure - se.exposure,
             xmax = beta.exposure + se.exposure)) +
  geom_smooth(method = "lm", color = "gray") +
  geom_errorbar() +
  geom_errorbarh() +
  geom_point(alpha = 0.7, size = 1) +
  facet_wrap(~ exposure + text, scales = "free", ncol = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  labs(
    title = "Traits -> CA",
    y = expression("GWAS"~beta~"(CA)"),
    x = expression("GWAS"~beta~"(Trait)"), color = "Included \nin MR"
  ) +
  theme_classic() +
  theme(text = element_text(size = 8))

p_mr1
