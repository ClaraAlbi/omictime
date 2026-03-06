#remotes::install_github("lcpilling/ukbrapR")
library(ukbrapR)
library(tidyr)
library(dplyr)
install.packages("broom")
library(purrr)


diag <- data.table::fread("/mnt/project/sleep_phase.csv")

ckd <- ukbrapR:::codes_df_ckd

ukbrapR:::codes_df_test

bp <- tribble(~condition, ~vocab_id, ~code,
              "bp", "ICD10","F31.2",
              # "bp", "ICD9", "296",
              # "bp","Read2", "E11..11",
              # "bp", "CTV3","E1176",
              "sleep", "ICD10", "G47",
              "insomnias", "ICD10", "G470",
              "hypersomnias", "ICD10", "G471",
              "sleep-wake schedule", "ICD10", "G472",
              "delayed", "ICD10", "G4721",
              "adv", "ICD10", "G4722",
              "apnea", "ICD10", "G473",
              "narcolepsy", "ICD10", "G474",
              "other", "ICD10" , "G478",
              "unspecified", "ICD10", "G479",
              "sleep", "ukb_noncancer", "1123")


# G47.0 Disorders of initiating and maintaining sleep [insomnias]113
# G47.1 Disorders of excessive somnolence [hypersomnias]137
# G47.2 Disorders of the sleep-wake schedule35
# G47.3 Sleep apnoea3847
# G47.4 Narcolepsy and cataplexy24
# G47.8 Other sleep disorders140
# G47.9 Sleep disorder, unspecified

diagnosis_list <- get_diagnoses(bp)
diagnosis_df <- get_df(diagnosis_list, group_by="condition")



phase %>%
  left_join(diagnosis_df) %>%
  left_join(job_vars) %>%
  group_by(`sleep-wake schedule_bin_prev`) %>% summarise(time_mean = mean(pred_mean, na.rm = T),
                                                         time_sd = sd(pred_mean, na.rm = T),
                                                    ca_mean = mean(res, na.rm = T),
                                                    ca_sd = sd(res, na.rm = T), n = n())

res <-
phase %>%
  left_join(diagnosis_df) %>%
  left_join(covs) %>%
  left_join(gen_covs) %>%
  pivot_longer(ends_with("_bin")) %>%
  group_by(name) %>% summarise(cases = sum(value == 1, na.rm = T)) %>% arrange(desc(cases))
  group_by(name) %>% nest() %>% ungroup() %>%
  #slice(1) %>%
  mutate(m = map(data, ~broom::tidy(lm(value ~ res + sex + age_recruitment + PC1 + PC2 + chrono, data = .x)))) %>% select(-data) %>% unnest(m)



summary(glm(sleep_bin_prev ~ res, data = phase %>% left_join(diagnosis_df)))

summary(lm(res ~ `sleep-wake schedule_bin_prev`, data = phase %>% left_join(diagnosis_df)))

summary(glm(apnea_bin_prev ~ res, data = phase %>% left_join(diagnosis_df)))


dis <- c(131850, # Date M06 first reported (other rheumatoid arthritis)
         131848, # Date M05 first reported (seropositive rheumatoid arthritis)

         131286, #	Date I10 first reported (essential (primary) hypertension)
         130814, #	Date E78 first reported (disorders of lipoprotein metabolism and other lipidaemias)

         130894, # 	Date F32 first reported (depressive episode)
         130896, # Date F33 first reported (recurrent depressive disorder)
         131306, # Date I25 first reported (chronic ischaemic heart disease)
         130792, #	Date E66 first reported (obesity)
         130708, #	Date E11 first reported (non-insulin-dependent diabetes mellitus)

         130842, #	Date F03 first reported (unspecified dementia)
         130892, #	Date F31 first reported (bipolar affective disorder)
         130874) #	Date F20 first reported (schizophrenia)

dis2 <- left_join(fread("/mnt/project/top_diseases_IEFG.csv"), fread("/mnt/project/immune.csv")) %>%
  select(eid, any_of(paste0("p",dis)))

outcome_labels <- c(
  ra_date = "Rheumatoid arthritis",
  hyp     = "Hypertension",
  dys     = "Dyslipidaemia",
  dep_ep  = "Depression (episode)",
  chd     = "Chronic heart disease",
  obe     = "Obesity",
  t2d     = "Type 2 diabetes",
  dem     = "Dementia",
  dep_rec = "Depression (recurrent)",
  bip     = "Bipolar disorder",
  scz     = "Schizophrenia",
  liv     = "Chronic liver disease"
)

saveRDS(dis2, "diseases_circadian.rds")



