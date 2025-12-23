library(tidyverse)

fields <- data.table::fread("data/field.tsv")

assay <- data.table::fread("data/olink_assay.dat") %>%
  mutate(Assay = toupper(Assay),
         Assay = str_replace(Assay, "-", "_"))

df_effects <- readRDS("data/combined_effects.rds") %>%
  mutate(phen = toupper(phen)) %>%
  left_join(assay, by = c("phen" = "Assay")) %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05) #%>%
  filter(type_clean == "Proteins")

df_r2 <- readRDS("data/combined_variance.rds") %>%
  mutate(phen = toupper(phen)) %>%
  filter(term == "time_day") %>%
  left_join(assay, by = c("phen" = "Assay")) %>%
  mutate(pval = p.adjust(p.value)) %>%
  filter(pval < 0.05) #%>%
  filter(type_clean == "Proteins")

length(which(!df_effects$phen %in% df_r2$phen))

length(which(!df_r2$phen %in% df_effects$phen))

secretome <- data.table::fread("../Downloads/proteinatlas.tsv") %>%
  inner_join(df_effects, by = c("Uniprot" = "UniProt")) %>%
  inner_join(df_r2, by = c("Uniprot" = "UniProt"))

list_p <- secretome %>%
  group_by(`Secretome location`) %>%
  summarise(n=n(),t = list(Gene))

n <- secretome %>%
  #filter(amplitude_24hfreq > 0.1) %>%
  filter(`Secretome location` != "") #%>%
  filter(`Secretome location` == "Secreted in brain")

secretome %>%
  count(`Secretome location`) %>% arrange(desc(n))

secretome %>%
  #filter(amplitude_24hfreq > 0.1) %>%
  #filter(t_r2 > 0.01) %>%
  #filter(`Secretome location` != "") %>%
  group_by(`Secretome location`) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(
    sec = paste0(`Secretome location`, "\nn=", n),
    sec = fct_reorder(sec, -n)
  ) %>%
  ggplot(aes(x = acrophase_24hfreq, color = sec)) +
  geom_density() +
  facet_wrap(~sec) +
  theme_classic() +
  theme(legend.position = "none")

secretome %>%
  ggplot(aes(x = acrophase_24hfreq, y = amplitude_24hfreq, color = `Secretome location`)) +
  geom_point() +
  coord_polar()
### 

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

lookout_effects <- df_effects %>%
  inner_join(df_r2) %>%
  filter(amplitude_24hfreq > 0.1) %>%
  filter(t_r2 > 0.01) %>%
  left_join(tissues) %>%
  count(category)


n_tissues <- tissues %>%
  #left_join(data.table::fread("../Downloads/rna_tissue_gtex_tissues.tsv"), by = c("tissue" ="Tissue"))
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
  left_join(tissues, by = c("phen" = "Gene")) %>%
  #filter(category %in% c("Tissue enriched", "Group enriched", NA)) %>%
  group_by(phen, acrophase_24hfreq, amplitude_24hfreq) %>% summarise(p = paste(tissue, collapse = ", ")) %>%
  mutate(p = paste0("(", p, ")")) %>%
  unite(col = "t", phen, p, sep = " ", remove = F)

unique(df$phen)

df %>% group_by(tissue) %>% count() %>% arrange(desc(n))

ptissue <- ggplot(df, aes(x = acrophase_24hfreq, y = amplitude_24hfreq, label = t)) +
  geom_point(size = 1) +
  coord_polar() +
  ggrepel::geom_label_repel(size = 2, max.overlaps = 50, box.padding = 0.1, alpha = 0.7) +
  #scale_color_manual(values = c( "#FC4E07","#00AFBB", "#E7B800")) +
  labs(x = "Acrophase", y = "Amplitude", color = "GTEX tissue") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.6)) +
  guides(color = guide_legend(ncol = 3)) +
  theme(legend.position = "none", text = element_text(size = 12))

ggsave("plots/FS_tissue_enrichments_circle.png", ptissue, width = 7, height = 7)



p_per <- prot_set %>%
  left_join(tissues, by = c("phen" = "Gene")) %>%
  group_by(tissue) %>%
  mutate(n = n()) %>% ungroup() %>%
  mutate(tissue_2 = fct_reorder(str_to_sentence(paste0(tissue, " n=",n)), -n)) %>%
  #select(acrophase_24hfreq, amplitude_24hfreq, phen, tissue_2, category) %>%
  ggplot(aes(x = acrophase_24hfreq, y = amplitude_24hfreq, label = phen, color = category)) +
  geom_point(size = 1) +
  coord_polar() +
  scale_y_continuous(limits = c(0, 0.8), n.breaks = 4) +
  scale_x_continuous(breaks = c(24, 6, 12, 18)) +
  labs(x = "Acrophase", y = "Amplitude", color = "GTEx category") +
  ggrepel::geom_label_repel(size = 2,
                            max.overlaps = 70,
                            box.padding = 0.1, label.padding = 0.1,
                            alpha = 0.7,
                            ) +
  facet_wrap(~tissue_2, ncol = 5) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3)
  )

ggsave("plots/FS_tissue_enrichments_circle2.png", p_per, width = 10, height = 15)





p_sec <- prot_set %>%
  left_join(secretome, by = c("UniProt" = "Uniprot")) %>%
  mutate(`Secretome location`  = case_when(`Secretome location` == "" ~ "Not assigned",
                                           TRUE ~ `Secretome location`)) %>%
  group_by(`Secretome location`) %>%
  mutate(n = n()) %>% ungroup() %>%
  mutate(tissue_2 = fct_reorder(str_to_sentence(paste0(`Secretome location`, " n=",n)), -n)) %>%
  ggplot(aes(x = acrophase_24hfreq, y = amplitude_24hfreq, label = phen, color = `Secretome location`)) +
  geom_point(size = 1) +
  coord_polar() +
  scale_y_continuous(limits = c(0, 0.6), n.breaks = 4) +
  scale_x_continuous(breaks = c(24, 6, 12, 18)) +
  labs(x = "Acrophase", y = "Amplitude", color = "Secretome category") +
  ggrepel::geom_label_repel(size = 2,
                            max.overlaps = 50,
                            box.padding = 0.1, label.padding = 0.1,
                            alpha = 0.7) +
  facet_wrap(~tissue_2, ncol = 3) +
  theme_classic() +
  theme(legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3)
  )
ggsave("plots/FS_tissue_secretions_circle.png", p_sec, width = 10, height = 15)



