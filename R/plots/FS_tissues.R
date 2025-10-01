library(tidyverse)

fields <- data.table::fread("data/field.tsv")

assay <- data.table::fread("data/olink_assay.dat") %>%
  mutate(Assay = toupper(Assay),
         Assay = str_replace(Assay, "-", "_"))

df_effects <- readRDS("data/combined_effects.rds") %>%
  mutate(phen = toupper(phen)) %>%
  left_join(assay, by = c("phen" = "Assay")) %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05) %>%
  filter(type_clean == "Proteins")

df_r2 <- readRDS("data/combined_variance.rds") %>%
  mutate(phen = toupper(phen)) %>%
  filter(term == "time_day") %>%
  left_join(assay, by = c("phen" = "Assay")) %>%
  mutate(pval = p.adjust(p.value)) %>%
  filter(pval < 0.05) %>%
  filter(type_clean == "Proteins")

tissues <- data.table::fread("data/explore_ukb.csv") %>%
  dplyr::rename(UniProt = `UniProt ID`,
         tissue_info = `Tissue Specificity`) %>%
  mutate(
    tissue_info = str_trim(tissue_info),
    category = ifelse(str_detect(tissue_info, ":"), str_extract(tissue_info, "^[^:]+"), tissue_info),
    tissues = ifelse(str_detect(tissue_info, ":"), str_extract(tissue_info, "(?<=: ).*"), NA)
  ) %>%
  separate_rows(tissues, sep = ",\\s*") %>%
  mutate(tissues = ifelse(is.na(tissues), category, tissues)) %>%
  dplyr::select(Gene, UniProt, `Protein Name`, category, tissue = tissues) %>%
  mutate(tissue = ifelse(str_detect(tissue, ":"), str_extract(tissue, "(?<=: ).*"), tissue)) %>%
  filter(UniProt %in% assay$UniProt)

# unique prots
unique_prots <- tissues %>%
  filter(category %in% c("Tissue enriched", "Group enriched")) %>% pull(Gene) %>% unique()

unique_prots_all <- assay %>% pull(UniProt) %>% unique()


n_tissues <- tissues %>%
  filter(category %in% c("Tissue enriched", "Group enriched")) %>%
  mutate(is_time = UniProt %in% df_effects$UniProt) %>%
  group_by(tissue, is_time) %>%
  count() %>%
  filter(!tissue %in% c("Low tissue specificity", "", "Low tissue specificity, Low tissue specificity", "Not detected")) %>%
  tidyr::pivot_wider(
    names_from  = is_time,
    values_from = n,
    names_prefix = "is_time",
    values_fill = 0
  ) %>%
  mutate(
    n         = is_timeTRUE + is_timeFALSE,
    frac_time = if_else(n > 0, is_timeTRUE / n, NA_real_)
  ) %>%
  filter(!is.na(frac_time)) %>%
  filter(is_timeFALSE > 0) %>%
  filter(is_timeTRUE > 0)

lvl <- n_tissues %>%
  arrange(desc(n), desc(frac_time)) %>%
  pull(tissue)

p_enrich <- n_tissues %>%
  pivot_longer(c(is_timeFALSE, is_timeTRUE), names_to = "name", values_to = "value") %>%
  mutate(tissue = factor(tissue, levels = rev(lvl)),
         name = case_when(name == "is_timeFALSE" ~ "No",
                          name == "is_timeTRUE" ~ "Yes")) %>%
  ggplot(aes(x = tissue, y = value, fill = name)) +
  geom_col() +
  labs(fill = "Rhythmic", y = "Tissue", x = "") +
  coord_flip() +
  scale_fill_manual(values = c("#E7B800","#00AFBB")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("plots/FS_tissue_enrichments.png", p_enrich, width = 6, height = 7)


prot_set <- df_effects %>%
  inner_join(df_r2) %>%
  filter(t_r2 > 0.01) %>%
  filter(amplitude_24hfreq > 0.1)

df <- prot_set %>%
  inner_join(tissues, by = c("phen" = "Gene")) %>%
  filter(category %in% c("Tissue enriched")) %>%
  group_by(tissue) %>% mutate(n = n()) %>% ungroup() %>%
  mutate(category = factor(category, levels = c("Tissue enriched", "Group enriched")))
         #tissue = case_when(n < 3 ~ "Other",
                            #TRUE ~ tissue))

df %>% group_by(tissue) %>% count() %>% arrange(desc(n))

ptissue <- ggplot(df, aes(x = acrophase_24hfreq, y = amplitude_24hfreq, color = tissue, label = paste0(toupper(phen), "_", tissue))) +
  geom_point(size = 3) +
  coord_polar() +
  ggrepel::geom_text_repel(size = 3.5, max.overlaps = 20, box.padding = 0.5) +
  #scale_color_manual(values = c( "#FC4E07","#00AFBB", "#E7B800")) +
  labs(x = "Acrophase", y = "Amplitude", color = "GTEX tissue") +
  theme_minimal() +
  guides(color = guide_legend(ncol = 3)) +
  theme(legend.position = "none", text = element_text(size = 12))

ggsave("plots/FS_tissue_enrichments_circle.png", ptissue, width = 6, height = 6)

ggplot(df, aes(x = acrophase_24hfreq, y = amplitude_24hfreq, shape = category, color = tissue)) +
  geom_point() +
  geom_text(data = df[df$amplitude_24hfreq > 0.3,],
            aes(label = paste0(phen, "_", tissue))) +
  coord_polar() +
  theme(legend.position = "none")

