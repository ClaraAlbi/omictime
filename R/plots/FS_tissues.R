library(tidyverse)

fields <- data.table::fread("data/field.tsv")

df_effects <- bind_rows(readRDS("data/effects_labs.rds") %>% mutate(type = "Biochemistry") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                        readRDS("data/effects_counts.rds") %>% mutate(type = "Cell_counts") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                        readRDS("data/effects_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                          left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("data/effects_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen)
  ) %>%
  mutate(color_var = case_when(type == "Proteomics-Olink" ~ "#76B041",
                               type == "Metabolomics-NMR" ~ "#2374AB",
                               type == "Cell_counts" ~ "#8F3985",
                               type == "Biochemistry" ~ "#E85F5C"),
         term = case_when(term == "(Intercept)" ~ "mesor",
                          term == "cos(2 * pi * time_day/24)" ~ "beta_cos1",
                          term == "sin(2 * pi * time_day/24)" ~ "beta_sin1")) %>%
  pivot_wider(id_cols = c(phen, color_var, type, title), values_from = c(estimate, p.value, std.error), names_from = term) %>%
  mutate(amplitude_24hfreq = sqrt(estimate_beta_cos1^2 + estimate_beta_sin1^2),
         acrophase_24hfreq = (atan2(estimate_beta_sin1, estimate_beta_cos1) / (2 * pi) * 24 + 24) %% 24,
         q = as.integer(round(acrophase_24hfreq, 0)))

df_r2 <- bind_rows(readRDS("data/aov_labs.rds") %>% mutate(type = "Biochemistry") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("data/aov_counts.rds") %>% mutate(type = "Cell_counts") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("data/aov_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("data/aov_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen)) %>%
  filter(term == "time_day")

tissues <- data.table::fread("~/Downloads/explore_ukb.csv") %>%
  rename(tissue_info = `Tissue Specificity`) %>%
  mutate(
    tissue_info = str_trim(tissue_info),
    category = ifelse(str_detect(tissue_info, ":"), str_extract(tissue_info, "^[^:]+"), tissue_info),
    tissues = ifelse(str_detect(tissue_info, ":"), str_extract(tissue_info, "(?<=: ).*"), NA)
  ) %>%
  separate_rows(tissues, sep = ",\\s*") %>%
  mutate(tissues = ifelse(is.na(tissues), category, tissues)) %>%
  select(Gene, `Protein Name`, category, tissue = tissues) %>%
  mutate(tissue = ifelse(str_detect(tissue, ":"), str_extract(tissue, "(?<=: ).*"), tissue), Gene = tolower(Gene))

df <- df_effects %>%
  left_join(df_r2) %>%
  filter(amplitude_24hfreq > 0.1 & pr2 > 0.01 & p.value < 0.05*3000) %>%
  inner_join(tissues, by = c("phen" = "Gene")) %>%
  filter(!tissue %in% c("Not detected", "Low tissue specificity", ""))

ptissue <- ggplot(df, aes(x = acrophase_24hfreq, y = amplitude_24hfreq, color = category, label = phen)) +
  geom_point() +
  coord_polar() +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  facet_wrap(~tissue, ncol = 6) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  labs(x = "Acrophase", y = "Amplitude", color = "GTEX category") +
  theme_minimal() +
  theme(legend.position = "bottom", text = element_text(size = 14))

ggsave("plots/tissue_enrichments.png", ptissue, width = 8, height = 10)



