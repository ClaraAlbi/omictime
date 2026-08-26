source("config/paths.R")

library(tidyverse)
library(data.table)

# Load 

res_dir <- bmrc_path("results")

cojo_files <- list.files(
  res_dir,
  pattern = "\\.jma\\.cojo$",
  full.names = TRUE
)

get_protein <- function(x) {
  x <- basename(x)
  x <- sub("_chr[0-9XY]+\\.jma\\.cojo$", "", x)
  x <- sub("_OID.*$", "", x)
  x <- sub("_.*$", "", x)
  toupper(x)
}

target_prots <- c(
  "MYOC","HYAL1","EFNA1","TNR","FAS","RELT","CD276",
  "SPON2","SPINK5","GDF15","PGF","ANGPTL1",
  "KLK12","LGALS1","HS3ST3B1", "KLK13"
  #"IGFBP4", "SERPINE1", "TNFRSF12A", "HSPB6", "SEMA3F", "PLAT", "HGF", "MB", "ADAMTS15", "TNFRSF1A", "FURIN", "IL6"
)

files_by_protein <- split(cojo_files, get_protein(cojo_files))

files_by_protein <- files_by_protein[
  names(files_by_protein) %in% target_prots
]

sapply(files_by_protein, length)

###  CA

gwas_CA <- fread("data/GWAS_res_v4/res_v4_GWAS.fastGWA") %>%
  rename(
    beta_CA = BETA,
    se_CA   = SE,
    AF1_CA  = AF1,
    P_CA    = P,
    N_CA    = N
  )


mr_results_ca <- files_by_protein %>%
  map(~ run_mr_one_protein(
    files_for_protein = .x,
    outcome_gwas = gwas_CA,
    outcome_name = "CA",
    outcome_beta = "beta_CA",
    outcome_se = "se_CA",
    outcome_eaf = "AF1_CA",
    outcome_p = "P_CA",
    outcome_n = "N_CA"
  )) %>%
  compact()

#saveRDS(mr_results_ca, "results/mr/MR_proteins_CA.rds")
mr_results_ca <- readRDS("results/mr/MR_proteins_CA.rds")

### RESULTS

mr_table_ca <- bind_rows(
  lapply(mr_results_ca, function(x) x$mr)
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
  mr_results_ca,
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

presso_table <- map_dfr(mr_results_ca, extract_presso)

final_table_ca <- mr_table_ca %>%
  left_join(egger_table, by = c("exposure", "outcome")) %>%
  left_join(presso_table, by = c("exposure", "outcome"))


harm_table <- map_dfr(
  mr_results_ca,
  function(x) {
    tryCatch({
      
      harm <- x$harm_snps
      if (is.null(harm) || nrow(harm) == 0) return(NULL)
      
      harm %>%
        mutate(
          protein = as.character(x$mr$protein[1]),
          outcome = as.character(x$mr$outcome[1])
        )
      
    }, error = function(e) NULL)
  }
)



p_mr1 <- harm_table %>%
  left_join(mr_table_ca, by = "protein") %>%
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
  facet_wrap(~ protein + text, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  labs(
    title = "Protein -> CA",
    y = expression("GWAS"~beta~"(CA)"),
    x = expression("GWAS"~beta~"(Protein)"), color = "Included \nin MR"
  ) +
  theme_classic() +
  theme(text = element_text(size = 8))

ggsave("OUTPUTS/MR_prots_CA.png", p_mr1, width = 10, height = 8)





### Chrono

gwas_chrono <- fread("data/chronotype/chrono_v2_GWAS.fastGWA") %>%
  rename(
    beta_chrono = BETA,
    se_chrono   = SE,
    AF1_chrono  = AF1,
    P_chrono    = P,
    N_chrono    = N
  )


mr_results_chrono <- files_by_protein %>%
  map(~ run_mr_one_protein(
    files_for_protein = .x,
    outcome_gwas = gwas_chrono,
    outcome_name = "Chronotype",
    outcome_beta = "beta_chrono",
    outcome_se = "se_chrono",
    outcome_eaf = "AF1_chrono",
    outcome_p = "P_chrono",
    outcome_n = "N_chrono"
  )) %>%
  compact()

#saveRDS(mr_results_chrono, "results/mr/MR_proteins_chrono.rds")


mr_results_chrono <- readRDS("results/mr/MR_proteins_CA.rds")

mr_table_chrono <- bind_rows(
  lapply(mr_results_chrono, function(x) x$mr)
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
  mr_results_chrono,
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

presso_table <- map_dfr(mr_results_chrono, extract_presso)

final_table_chrono <- mr_table_chrono %>%
  left_join(egger_table, by = c("exposure", "outcome")) %>%
  left_join(presso_table, by = c("exposure", "outcome"))


harm_table <- map_dfr(
  mr_results_ca,
  function(x) {
    tryCatch({
      
      harm <- x$harm_snps
      if (is.null(harm) || nrow(harm) == 0) return(NULL)
      
      harm %>%
        mutate(
          protein = as.character(x$mr$protein[1]),
          outcome = as.character(x$mr$outcome[1])
        )
      
    }, error = function(e) NULL)
  }
)



p_mr2 <- harm_table %>%
  left_join(mr_table_ca, by = "protein") %>%
  filter(method == "Inverse variance weighted (multiplicative random effects)") %>%
  mutate(text = paste0("n=", nsnp, ", bIVW=", round(b, 2), ", pFDR=", formatC(pval_FDR, format = "e", digits = 0))) %>%
  ggplot(aes(x = beta.exposure, y = beta.outcome, color = mr_keep)) +
  geom_point(alpha = 0.8) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = FALSE,
    colour = scales::alpha("black", 0.2)
  ) +
  facet_wrap(~ protein + text, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  #scale_x_continuous(limits = c(-0.3, 0.3)) +
  #scale_y_continuous(limits = c(-0.3, 0.3)) +
  labs(
    x = expression("GWAS"~beta~"(Protein)"),
    y = expression("GWAS"~beta~"(Chronotype)"), color = "Included \nin MR"
  ) +
  theme_classic() +
  theme(text = element_text(size = 8))

ggsave("OUTPUTS/MR_prots_chrono.png", p_mr2, width = 10, height = 8)




# SYN2 (chr3)
# Start: 12,045,876  − 1,000,000 = 11,045,876  
# End:   12,232,900  + 1,000,000 = 13,232,900  
# Region: chr3: 11,045,876–13,232,900
# 
# COLEC10 (chr8)
# Start: 120,007,691 − 1,000,000 = 119,007,691  
# End:   120,118,821 + 1,000,000 = 121,118,821  
# Region: chr8: 119,007,691–121,118,821
# 
# NLRP12 (chr19)
# Start: 54,296,857 − 1,000,000 = 53,296,857  
# End:   54,327,648 + 1,000,000 = 55,327,648  
# Region: chr19: 53,296,857–55,327,648


cis_syn2 <- comb %>%
  filter(Chr == 3 & bp >= 11045876 & bp <= 13232900)

cis_colec10 <- comb %>%
  filter(Chr == 8 & bp >= 119007691 & bp <= 121118821)

cis_nlrp12 <- comb %>%
  filter(Chr == 19 & bp >= 53296857 & bp <= 55327648)

cis_c1orf220 <- comb %>%
  filter(
    Chr == 1,
    bp >= 177511887,
    bp <= 179518024
  )
