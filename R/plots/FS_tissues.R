library(tidyverse)

fields <- data.table::fread("data/field.tsv")

df_effects <- readRDS("data/combined_effects.rds") %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05) %>%
  filter(type_clean == "Proteins")

df_r2 <- readRDS("data/combined_variance.rds") %>%
  filter(term == "time_day") %>%
  mutate(pval = p.adjust(p.value)) %>%
  filter(pval < 0.05) %>%
  filter(type_clean == "Proteins")

prot_set <- df_effects %>%
  inner_join(df_r2) %>%
  filter(t_r2 > 0.01) %>%
  filter(amplitude_24hfreq > 0.1)

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

n_tissues <- tissues %>%
  filter(category %in% c("Tissue enriched", "Group enriched"))  %>%
  mutate(is_time = Gene %in% df_effects$phen) %>%
  group_by(tissue, is_time) %>%
  count() %>%
  filter(!tissue %in% c("Low tissue specificity", "", "Low tissue specificity, Low tissue specificity", "Not detected")) %>%
  pivot_wider(names_from = is_time, values_from = n, names_prefix = "is_time") %>%
  mutate(n = is_timeTRUE + is_timeFALSE,
         frac_time = is_timeTRUE/n)

sum(n_tissues$n)

df <- prot_set %>%
  inner_join(tissues, by = c("phen" = "Gene")) %>%
  filter(category %in% c("Tissue enriched", "Group enriched")) %>%
  group_by(tissue) %>% mutate(n = n()) %>% ungroup() %>%
  mutate(category = factor(category, levels = c("Tissue enriched", "Group enriched")),
         tissue = fct_reorder(factor(tissue), n, .desc = TRUE))


df %>% arrange(desc(n)) %>% count(tissue)

# tissue                  n
# 1 intestine             8
# 2 brain                 4
# 3 pancreas              4
# 4 adrenal gland         3
# 5 liver                 3
# 6 lymphoid tissue       3
# 7 parathyroid gland     3
# 8 pituitary gland       3
# 9 skeletal muscle       3

df %>% filter(category == "Tissue enriched") %>%
  #filter(amplitude_24hfreq > 0.2) %>%
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

ggplot(df, aes(x = acrophase_24hfreq, y = amplitude_24hfreq, shape = category, color = tissue)) +
  geom_point() +
  geom_text(data = df[df$amplitude_24hfreq > 0.3,],
            aes(label = paste0(phen, "_", tissue))) +
  coord_polar() +
  theme(legend.position = "none")

