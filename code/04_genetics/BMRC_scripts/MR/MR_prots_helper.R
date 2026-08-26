
run_mr_one_protein <- function(
    files_for_protein,
    outcome_gwas,
    outcome_name,
    outcome_beta,
    outcome_se,
    outcome_eaf,
    outcome_p,
    outcome_n
) {
  
  protein_name <- get_protein(files_for_protein[1])
  message("Running MR for ", protein_name, " -> ", outcome_name)
  
  # ------------------
  # Read & combine exposure instruments
  # ------------------
  exposure <- bind_rows(lapply(files_for_protein, fread))
  
  if (nrow(exposure) < 3) {
    message("Skipping ", protein_name, ": <3 instruments")
    return(NULL)
  }
  
  # ------------------
  # Join to outcome GWAS
  # ------------------
  dat <- exposure %>%
    left_join(outcome_gwas, by = "SNP") %>%
    filter(
      !is.na(.data[[outcome_beta]]),
      !is.na(.data[[outcome_se]])
    )
  
  if (nrow(dat) < 3) {
    message("Skipping ", protein_name, ": no overlap with outcome")
    return(NULL)
  }
  
  # ------------------
  # Format exposure
  # ------------------
  exposure_dat <- TwoSampleMR::format_data(
    dat %>% 
      as_tibble() %>% 
      mutate(
        Phenotype = protein_name),
    type = "exposure",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "refA",
    eaf_col = "freq",
    pval_col = "p",
    samplesize_col = "n"
  )
  
  # ------------------
  # Format outcome
  # ------------------
  outcome_dat <- TwoSampleMR::format_data(
    dat %>% 
      as_tibble() %>% 
      mutate(Phenotype = outcome_name),
    type = "outcome",
    snp_col = "SNP",
    beta_col = outcome_beta,
    se_col = outcome_se,
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = outcome_eaf,
    pval_col = outcome_p,
    samplesize_col = outcome_n
  )
  
  # ------------------
  # Harmonise
  # ------------------
  harm <- TwoSampleMR::harmonise_data(
    exposure_dat,
    outcome_dat,
    action = 2
  )
  
  if (nrow(harm) < 3) {
    message("Skipping ", protein_name, ": <3 harmonised SNPs")
    return(NULL)
  }
  
  # ------------------
  # Main MR
  # ------------------
  mr_res <- TwoSampleMR::mr(
    harm,
    method_list = c(
      "mr_ivw_mre",
      "mr_egger_regression",
      "mr_weighted_median"
    )
  ) %>%
    mutate(
      protein = protein_name,
      outcome = outcome_name,
      nsnp = nrow(harm)
    )
  
  # ------------------
  # MR-Egger intercept
  # ------------------
  egger_int <- tryCatch(
    {
      TwoSampleMR::mr_pleiotropy_test(harm) %>%
        mutate(
          protein = protein_name,
          outcome = outcome_name,
          nsnp = nrow(harm)
        )
    },
    error = function(e) NULL
  )
  
  # ------------------
  # MR-PRESSO (requires >=4 SNPs)
  # ------------------
  presso_res <- NULL
  
  if (nrow(harm) >= 4) {
    presso_res <- tryCatch(
      {
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          data = harm,
          NbDistribution = 1000,
          SignifThreshold = 0.05
        )
      },
      error = function(e) NULL
    )
  }
  
  list(
    harm_snps = harm,
    mr = mr_res,
    egger_intercept = egger_int,
    presso = presso_res
  )
}

extract_presso <- function(x, outlier_alpha = 0.05) {
  tryCatch({
    
    if (is.null(x$presso)) return(NULL)
    if (is.null(x$mr) || nrow(x$mr) == 0) return(NULL)
    
    res <- x$presso$`MR-PRESSO results`
    
    # helper to parse "<0.001"
    parse_p <- function(p) {
      if (length(p) == 0 || is.na(p)) return(NA_real_)
      p <- as.character(p)
      if (grepl("^<", p)) return(as.numeric(sub("<", "", p)))
      suppressWarnings(as.numeric(p))
    }
    
    # Global / distortion
    global_p <- parse_p(res$`Global Test`$Pvalue)
    distortion_p <- parse_p(res$`Distortion Test`$Pvalue)
    
    # Outliers
    outlier_tbl <- res$`Outlier Test`
    outlier_pvals <- sapply(outlier_tbl$Pvalue, parse_p)
    n_outliers <- sum(outlier_pvals < outlier_alpha, na.rm = TRUE)
    
    # Robust extraction of labels
    exposure <- if ("exposure" %in% names(x$mr)) x$mr$exposure[1] else NA_character_
    outcome  <- if ("outcome"  %in% names(x$mr)) x$mr$outcome[1]  else NA_character_
    
    data.frame(
      exposure = exposure,
      outcome = outcome,
      presso_global_p = global_p,
      presso_distortion_p = distortion_p,
      presso_outlier_n = n_outliers,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) NULL)
}



run_mr_exposure_to_protein <- function(
    file,
    exposure_gwas,
    exposure_name,
    exposure_beta,
    exposure_se,
    exposure_eaf,
    exposure_p,
    exposure_n,
    exposure_A1,
    exposure_A2
) {
  
  protein <- sub("\\.tsv$", "", basename(file))
  message("Running MR: ", exposure_name, " → ", protein)
  
  # ------------------
  # Read protein GWAS (OUTCOME)
  # ------------------
  protein_gwas <- fread(
    file,
    col.names = c("SNP", "A1", "A2", "freq", "beta", "se", "p", "N")
  )
  
  if (nrow(protein_gwas) < 3) return(NULL)
  
  # ------------------
  # Join to exposure GWAS
  # ------------------
  dat <- exposure_gwas %>%
    left_join(protein_gwas, by = "SNP") %>%
    filter(
      !is.na(.data[[exposure_beta]]),
      !is.na(.data[[exposure_se]]),
      !is.na(beta),
      !is.na(se)
    )
  
  if (nrow(dat) < 3) return(NULL)
  
  dat <- as.data.frame(dat)
  
  # ------------------
  # Format EXPOSURE
  # ------------------
  exposure_dat <- TwoSampleMR::format_data(
    dat,
    type = "exposure",
    snp_col = "SNP",
    beta_col = exposure_beta,
    se_col = exposure_se,
    effect_allele_col = exposure_A1,
    other_allele_col  = exposure_A2,
    eaf_col = exposure_eaf,
    pval_col = exposure_p,
    samplesize_col = exposure_n,
    phenotype_col = exposure_name
  )
  
  # ------------------
  # Format OUTCOME
  # ------------------
  outcome_dat <- TwoSampleMR::format_data(
    dat,
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col  = "A2",
    eaf_col = "freq",
    pval_col = "p",
    samplesize_col = "N",
    phenotype_col = protein
  )
  
  # ------------------
  # Harmonise
  # ------------------
  harm <- TwoSampleMR::harmonise_data(
    exposure_dat,
    outcome_dat,
    action = 2
  )
  
  if (nrow(harm) < 3) return(NULL)
  
  # ------------------
  # Main MR
  # ------------------
  mr_res <- TwoSampleMR::mr(
    harm,
    method_list = c(
      "mr_ivw_mre",
      "mr_weighted_median",
      "mr_egger_regression"
    )
  ) %>%
    mutate(
      exposure = exposure_name,
      outcome = protein,
      nsnp = nrow(harm)
    )
  
  # ------------------
  # MR-Egger intercept
  # ------------------
  egger_int <- tryCatch(
    {
      TwoSampleMR::mr_pleiotropy_test(harm) %>%
        mutate(
          exposure = exposure_name,
          outcome = protein,
          nsnp = nrow(harm)
        )
    },
    error = function(e) NULL
  )
  
  # ------------------
  # MR-PRESSO (≥4 SNPs)
  # ------------------
  presso_res <- NULL
  
  if (nrow(harm) >= 4) {
    presso_res <- tryCatch(
      {
        MRPRESSO::mr_presso(
          BetaOutcome = "beta.outcome",
          BetaExposure = "beta.exposure",
          SdOutcome = "se.outcome",
          SdExposure = "se.exposure",
          OUTLIERtest = TRUE,
          DISTORTIONtest = TRUE,
          data = harm,
          NbDistribution = 1000,
          SignifThreshold = 0.05
        )
      },
      error = function(e) NULL
    )
  }
  
  # ------------------
  # Return
  # ------------------
  list(
    harm_snps = harm,
    mr = mr_res,
    egger_intercept = egger_int,
    presso = presso_res
  )
}

