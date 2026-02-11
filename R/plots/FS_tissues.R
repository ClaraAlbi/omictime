library(tidyverse)

fields <- data.table::fread("data/field.tsv")

assay <- data.table::fread("data/olink_assay.dat") %>%
  mutate(Assay = toupper(Assay),
         Assay = str_replace(Assay, "-", "_"))

df_r2 <- readRDS("data/combined_variance.rds") %>%
  mutate(phen = toupper(phen)) %>%
  mutate(pfdr = p.adjust(p.value)) %>%
  filter(pfdr < 0.05) %>%
  filter(term == "time_day") %>%
  filter(t_r2 > 0.01)

df_effects <- readRDS("data/combined_effects.rds") %>%
  mutate(phen = toupper(phen)) %>%
  mutate(pfdr = p.adjust(pvalue_h)) %>%
  filter(pfdr < 0.05) %>%
  filter(type_clean == "Proteins")

d <- inner_join(df_effects, df_r2, by = "phen") %>%
  left_join(assay, by = c("phen" = "Assay"))


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
  filter(UniProt %in% assay$UniProt) #%>% #count(category)
  filter(category %in% c("Tissue enriched", "Group enriched")) #%>%
  filter(!tissue %in% c("Low tissue specificity", "", "Low tissue specificity, Low tissue specificity", "Not detected"))

length(unique(tissues$UniProt))

length(unique(tissues$UniProt[tissues$Gene %in% df_effects$phen]))

lookout_effects <- tissues %>%
  filter(Gene %in% df_effects$phen) %>%
  count(category)

d %>%
  left_join(tissues) %>%
  group_by(tissue)


n_tissues <- tissues %>%
  mutate(is_time = Gene %in% df_effects$phen) %>%
  group_by(tissue, is_time) %>%
  count() %>%
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

plot_d <- n_tissues %>%
  pivot_longer(c(is_timeFALSE, is_timeTRUE), names_to = "name", values_to = "value") %>%
  mutate(tissue = factor(tissue, levels = rev(lvl)),
         name = case_when(name == "is_timeFALSE" ~ "No",
                          name == "is_timeTRUE" ~ "Yes"))

p_enrich <- plot_d %>%
  ggplot(aes(x = tissue, y = value, fill = name)) +
  geom_col() +
  geom_text(data = plot_d %>% filter(name == "Yes"), aes(y = value + 10, label = round(frac_time, 2))) +
  labs(fill = "Rhythmic", x = "GTEx tissue", y = "Number of proteins (Tissue enriched & Group enriched)") +
  coord_flip() +
  scale_fill_manual(values = c("#E7B800","#00AFBB")) +
  theme_minimal() +
  theme(legend.position = "bottom", axis.title.x = element_text(size = 10))

ggsave("plots/FS_tissue_enrichments.png", p_enrich, width = 6, height = 7)


prot_set <- d

df <- prot_set %>%
  filter(type_clean.x == "Proteins") %>%
  left_join(tissues, by = c("phen" = "Gene")) %>%
  filter(!is.na(tissue)) %>%
  filter(category %in% c("Tissue enriched", "Group enriched")) %>% #%>% count(tissue) %>% arrange(desc(n))
  group_by(phen,`Protein Name`, acrophase_24hfreq, amplitude_24hfreq, category) %>% summarise(p = paste(tissue, collapse = ", ")) %>%
  mutate(p2 = paste0("(", p, ")")) %>%
  unite(col = "t", phen, p2, sep = " ", remove = F) %>%
  select(acrophase_24hfreq, amplitude_24hfreq, t)

unique(df$phen)

df %>% group_by(p) %>% count() %>% arrange(desc(n))

ptissue <- ggplot(df, aes(x = acrophase_24hfreq, y = amplitude_24hfreq, label = t)) +
  geom_point(size = 1) +
  coord_polar() +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 50, box.padding = 0.1, alpha = 0.7) +
  #scale_color_manual(values = c( "#FC4E07","#00AFBB", "#E7B800")) +
  labs(x = "Acrophase", y = "Amplitude", color = "GTEX tissue") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 24), breaks = 0:23,
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  guides(color = guide_legend(ncol = 3)) +
  theme(legend.position = "none", text = element_text(size = 12))
ptissue
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
p_per
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



