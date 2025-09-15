library(tidyverse)

fields <- data.table::fread("data/field.tsv")

df_effects <- bind_rows(readRDS("data/effects_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen)
  ) %>%
  mutate(term = case_when(term == "(Intercept)" ~ "mesor",
                          term == "cos(2 * pi * time_day/24)" ~ "beta_cos1",
                          term == "sin(2 * pi * time_day/24)" ~ "beta_sin1")) %>%
  pivot_wider(id_cols = c(phen, type, title), values_from = c(estimate, p.value, std.error), names_from = term) %>%
  mutate(amplitude_24hfreq = sqrt(estimate_beta_cos1^2 + estimate_beta_sin1^2),
         acrophase_24hfreq = (atan2(estimate_beta_sin1, estimate_beta_cos1) / (2 * pi) * 24 + 24) %% 24,
         q = as.integer(round(acrophase_24hfreq, 0)))

df_r2 <- bind_rows(readRDS("data/aov_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen)) %>%
  filter(term == "time_day") %>%
  left_join(df_effects) %>%
  mutate(p_val = p.adjust(p.value)) %>%
  filter(amplitude_24hfreq > 0.1 & pr2 > 0.01 & p_val < 0.05)

tissues <- data.table::fread("data/explore_ukb.csv") %>%
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

df <- df_r2 %>%
  inner_join(tissues, by = c("phen" = "Gene")) %>%
  filter(!tissue %in% c("Not detected", "Low tissue specificity", "")) %>%
  group_by(tissue) %>% mutate(n = n()) %>% ungroup() %>%
  mutate(category = factor(category, levels = c("Tissue enriched", "Group enriched", "Tissue enhanced")),
         tissue = fct_reorder(factor(tissue), n, .desc = TRUE))

df %>% arrange(desc(n)) %>% count(tissue)
# tissue                n
# <fct>             <int>
#   1 liver                13
# 2 intestine            11
# 3 lymphoid tissue       9
# 4 adipose tissue        8
# 5 esophagus             7
# 6 brain                 6
# 7 pancreas              6

df %>% filter(category == "Tissue enriched") %>%
  filter(amplitude_24hfreq > 0.2) %>%
  select(phen, tissue) %>%
  arrange(desc(tissue)) %>% mutate(phen = toupper(phen))

ptissue <- ggplot(df, aes(x = acrophase_24hfreq, y = amplitude_24hfreq, color = category, label = phen)) +
  geom_point() +
  coord_polar() +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
  facet_wrap(~tissue, ncol = 6) +
  scale_color_manual(values = c( "#FC4E07","#00AFBB", "#E7B800")) +
  labs(x = "Acrophase", y = "Amplitude", color = "GTEX category") +
  theme_minimal() +
  theme(legend.position = "bottom", text = element_text(size = 14))

ggsave("plots/tissue_enrichments.png", ptissue, width = 9, height = 11)



