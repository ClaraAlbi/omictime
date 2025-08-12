library(tidyverse)
source("R/analyses/split_sex_age_r2.R")

fields <- data.table::fread("data/field.tsv")

# VARIANCE RESULTS

df_r2 <- bind_rows(readRDS("data/aov_labs.rds") %>% mutate(type = "Biochemistry") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("data/aov_counts.rds") %>% mutate(type = "Cell_counts") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("data/aov_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("data/aov_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = toupper(phen))
  ) %>%
  mutate(color_var = case_when(type == "Proteomics-Olink" ~ "#76B041",
                               type == "Metabolomics-NMR" ~ "#2374AB",
                               type == "Cell_counts" ~ "#8F3985",
                               type == "Biochemistry" ~ "#E85F5C"),
         type_clean = case_when(type == "Proteomics-Olink" ~ "Proteins",
                                type == "Metabolomics-NMR" ~ "Metabolites",
                                type == "Cell_counts" ~ "Cell counts",
                                type == "Biochemistry" ~ "Biochemistry"),
         title = case_when(title == "White blood cell (leukocyte) count" ~ "Leukocyte count",
                           title == "Phospholipids to Total Lipids in Small HDL percentage" ~ "Phosphlipid ratio SHDL",
                           title == "Cholesterol to Total Lipids in Very Large HDL percentage" ~ "Cholesterol ratio VLHDL",
                           title == "Phospholipids to Total Lipids in Very Large HDL percentage" ~ "Phospholipids ratio VLHDL",
                           title == "Cholesterol to Total Lipids in Small HDL percentage" ~ "Cholesterol ratio SHLD",
                           title == "Spectrometer-corrected alanine" ~ "Alanine",
                           TRUE ~ title))

# Reformat labels and r2

df_adj <- split_sex_age_r2(df_r2)

df_comb <- df_adj %>%
  group_by(type, phen, color_var, type_clean ) %>%
  ungroup() %>%
  distinct(phen, .keep_all = TRUE) %>%
  arrange(desc(pr2)) %>%
  select(phen) %>%
  inner_join(df_adj, by = "phen") %>%
  filter(term != "Residuals") %>%
  group_by(phen, type, term, title, color_var, type_clean ) %>%
  summarise(t_r2 = sum(pr2), pvalue = min(p.value)) %>%
  mutate(term = ifelse(grepl("^PC[0-9]+$", term), "PCs", term)) %>%
  group_by(phen, type, title, term, color_var, type_clean ) %>%
  summarise(pr2 = sum(t_r2), p.value = min(pvalue), .groups = "drop")

saveRDS(df_comb, "data/combined_variance.rds")

######### Save to github

#
# variance_table <- df_comb %>%
#   select(FID = phen, Name = title, Type = type_clean, Covariate = term, R2 = pr2, pvalue = p.value) %>%
#   mutate(R2 = round(R2, 5))
#   #group_by(Type, FID, Covariate) %>% arrange(desc(R2)) %>%
#   mutate(across(c(R2, pvalue), ~ ))
#
#
# mutate(across(where(is.numeric), ~ sprintf("%.15g", .x)))
#

#### EFFECTS

df_effects <- bind_rows(readRDS("data/effects_labs.rds") %>% mutate(type = "Biochemistry") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                        readRDS("data/effects_counts.rds") %>% mutate(type = "Cell_counts") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                        readRDS("data/effects_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("data/effects_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = toupper(phen))
  ) %>%
  mutate(color_var = case_when(type == "Proteomics-Olink" ~ "#76B041",
                               type == "Metabolomics-NMR" ~ "#2374AB",
                               type == "Cell_counts" ~ "#8F3985",
                               type == "Biochemistry" ~ "#E85F5C"),
         type_clean = case_when(type == "Proteomics-Olink" ~ "Proteins",
                                type == "Metabolomics-NMR" ~ "Metabolites",
                                type == "Cell_counts" ~ "Cell counts",
                                type == "Biochemistry" ~ "Biochemistry"),
         term = case_when(term == "(Intercept)" ~ "mesor",
                          term == "cos(2 * pi * time_day/24)" ~ "beta_cos1",
                          term == "sin(2 * pi * time_day/24)" ~ "beta_sin1"),
         title = case_when(title == "White blood cell (leukocyte) count" ~ "Leukocyte count",
                           title == "Phospholipids to Total Lipids in Small HDL percentage" ~ "Phosphlipid ratio SHDL",
                           title == "Cholesterol to Total Lipids in Very Large HDL percentage" ~ "Cholesterol ratio VLHDL",
                           title == "Phospholipids to Total Lipids in Very Large HDL percentage" ~ "Phospholipids ratio VLHDL",
                           title == "Cholesterol to Total Lipids in Small HDL percentage" ~ "Cholesterol ratio SHLD",
                           title == "Spectrometer-corrected alanine" ~ "Alanine",
                           TRUE ~ title)) %>%
  pivot_wider(id_cols = c(phen, color_var, type_clean, title), values_from = c(estimate, p.value, std.error), names_from = term) %>%
  mutate(amplitude_24hfreq = sqrt(estimate_beta_cos1^2 + estimate_beta_sin1^2),
         acrophase_24hfreq = (atan2(estimate_beta_sin1, estimate_beta_cos1) / (2 * pi) * 24 + 24) %% 24,
         q = as.integer(round(acrophase_24hfreq, 0)),
         pvalue_h = pmin(p.value_beta_cos1, p.value_beta_sin1, na.rm = TRUE))

saveRDS(df_effects, "data/combined_effects.rds")

