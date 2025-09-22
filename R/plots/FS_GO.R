library(tidyverse)

assay <- data.table::fread("data/olink_assay.dat") %>%
  mutate(Assay = toupper(Assay),
         Assay = str_replace(Assay, "-", "_"))

df_effects <- readRDS("data/combined_effects.rds") %>% mutate(phen = toupper(phen)) %>%
  left_join(assay, by = c("phen" = "Assay")) %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05) %>%
  filter(type_clean == "Proteins")

go <- data.table::fread("data/explore_ukb (1).csv") %>%
  rename(UniProt = `UniProt ID`) %>%
  mutate(
    go_bp = str_trim(`Biological Process`))  %>%
  separate_rows(go_bp, sep = ",\\s*") %>%
  filter(go_bp != "")

background <- go %>%
  #filter(UniProt %in% assay$UniProt) %>%
  group_by(go_bp) %>%
  #filter(str_detect(go_bp, "circadi")) %>%
  summarise(l = list(Gene),  n = n(), type = "bg")

rhythmic <- go %>%
  #filter(str_detect(go_bp, "circadi")) %>%
  filter(UniProt %in% df_effects$UniProt) %>% group_by(go_bp) %>%
  summarise(l = list(Gene),  n = n(), type = "r")

# Enrichment

proc <- bind_rows(background, rhythmic) %>%
  select(-l) %>%
  pivot_wider(names_from = type, values_from = n)

