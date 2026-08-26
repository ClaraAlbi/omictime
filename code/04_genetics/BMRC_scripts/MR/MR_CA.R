library(data.table)
library(dplyr)
library(purrr)

gwas_CA <- fread("data/GWAS_res_v4/res_v4_GWAS.fastGWA") %>%
  rename(
    SNP = SNP,
    beta_CA = BETA,
    se_CA   = SE,
    AF1_CA  = AF1,
    P_CA    = P,
    N_CA    = N,
    A1_CA   = A1,
    A2_CA   = A2
  ) %>%
  as_tibble()

mr_files <- list.files(
  "protein_CA_overlap",
  pattern = "\\.tsv$",
  full.names = TRUE
)

target_prots <- c(
  "MYOC","HYAL1","EFNA1","TNR","FAS","RELT","CD276",
  "SPON2","SPINK5","GDF15","PGF","ANGPTL1",
  "KLK12","LGALS1","HS3ST3B1", "KLK13"
)

get_protein <- function(x) {
  x <- basename(x)
  x <- sub("_chr[0-9XY]+\\.jma\\.cojo$", "", x)
  x <- sub("_OID.*$", "", x)
  x <- sub("_.*$", "", x)
  toupper(x)
}


mr_results_CA <- lapply(
  mr_files[which(get_protein(mr_files) %in% target_prots)],
  run_mr_exposure_to_protein,
  exposure_gwas = gwas_CA,
  exposure_name = "CA",
  exposure_beta = "beta_CA",
  exposure_se   = "se_CA",
  exposure_eaf  = "AF1_CA",
  exposure_p    = "P_CA",
  exposure_n    = "N_CA",
  exposure_A1   = "A1_CA",
  exposure_A2   = "A2_CA"
)

#saveRDS(mr_results_CA, "results/mr/MR_CA_proteins.rds")

mr_results_CA <- readRDS("results/mr/MR_CA_proteins.rds")

mr_table_ca <- bind_rows(
  lapply(mr_results_CA, function(x) x$mr)
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
  mr_results_CA,
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

presso_table <- map_dfr(mr_results_CA, extract_presso)

final_table_ca <- mr_table_ca %>%
  left_join(egger_table, by = c("exposure", "outcome")) %>%
  left_join(presso_table, by = c("exposure", "outcome"))



harm_table <- map_dfr(
  mr_results_CA,
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


p_mrX <- harm_table %>%
  left_join(mr_table_ca, by = c("exposure", "outcome")) %>%
  filter(method == "IVW") %>%
  mutate(text = paste0("n=", nsnp, ", bIVW=", round(b, 2), ", pFDR=", formatC(pval_FDR, format = "e", digits = 0)),
         protein = map_chr(outcome, ~str_split(.x, "_")[[1]][1])) %>%
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
    title = "CA -> Protein",
    y = expression("GWAS"~beta~"(Protein)"),
    x = expression("GWAS"~beta~"(CA)"), color = "Included \nin MR"
  ) +
  theme_classic() +
  theme(text = element_text(size = 8))

ggsave("OUTPUTS/MR_CA_prots.png", p_mrX, width = 10, height = 8)





