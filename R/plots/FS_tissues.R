library(tidyverse)
library(ggtext)
library(tidytext)

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






d <- inner_join(df_effects, df_r2, by = "phen") %>%
  left_join(assay, by = c("phen" = "Assay"))


prot_set <- d

df <- prot_set %>%
  filter(type_clean.x == "Proteins") %>%
  left_join(tissues, by = c("phen" = "Gene")) %>%
  filter(!is.na(tissue)) %>%
  filter(category %in% c("Tissue enriched", "Group enriched")) %>% #%>% count(tissue) %>% arrange(desc(n))
  #group_by(phen,`Protein Name`, acrophase_24hfreq, amplitude_24hfreq, category) %>% summarise(p = paste(tissue, collapse = ",\n")) %>%
  #mutate(p2 = paste0("(", p, ")")) %>%
  #unite(col = "t", phen, p2, sep = " ", remove = F) %>%
  select(phen, acrophase_24hfreq, amplitude_24hfreq, tissue, t_r2)

unique(df$phen)

df %>% group_by(phen) %>% count() %>% arrange(desc(n))

ptissue <- ggplot(df, aes(x = acrophase_24hfreq, y = amplitude_24hfreq, label = phen, color = tissue, size = t_r2)) +
  geom_point() +
  #coord_polar() +
  ggrepel::geom_text_repel(size = 2.5, max.overlaps = 50, box.padding = 0.1, alpha = 0.7) +
  #scale_color_manual(values = c( "#FC4E07","#00AFBB", "#E7B800")) +
  labs(x = "Acrophase", y = "Amplitude", color = "GTEX tissue") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 24), breaks = 0:23,
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0.15, 0.6)) +
  guides(color = guide_legend(ncol = 7)) +
  theme_classic(base_size = 10) +
  theme(panel.grid.major = element_line(color = "gray"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.direction = "horizontal",
        #legend.justification = c("center", "top"),
        legend.background = element_rect(
          color = "black", fill = "white", linewidth = 0.2
        ),
        #axis.text.x = element_text(size = 14),
        #axis.title = element_text(size = 16),
        axis.line = element_blank())


ptissue
ggsave("plots/FS_tissue_enrichments_circle.png", ptissue, width = 7, height = 7)



system_map <- c(
  # CNS
  "Brain" = "CNS",
  "Choroid plexus" = "CNS",
  "Retina" = "CNS",

  "Adrenal gland" = "Endocrine",
  "Pituitary gland" = "Endocrine",
  "Thyroid gland" = "Endocrine",
  "Parathyroid gland" = "Endocrine",

  # GI & Metabolic
  "Intestine" = "Gastro",
  "Stomach 1" = "Gastro",
  "Esophagus" = "Gastro",
  "Liver" = "Gastro",
  "Gallbladder" = "Gastro",
  "Pancreas" = "Gastro",
  "Adipose tissue" = "Gastro",
  "Tongue" = "Gastro",

  # Endocrine

  # Other
  "Lymphoid tissue" = "Other",
  "Bone marrow" = "Other",
  "Kidney" = "Other",
  "Lung" = "Other",
  "Epididymis" = "Other",
  "Heart muscle" = "Other",
  "Skeletal muscle" = "Other"
)

# authoritative orders
tissue_order <- names(system_map)
system_order <- unique(system_map)

df2 <- df %>%
  mutate(
    tissue = str_to_sentence(tissue),
    system = recode(tissue, !!!system_map)
  )

# -----------------------------
# 3. BUILD ALTERNATING COLORS
#    (RESTART WITHIN EACH SYSTEM)
# -----------------------------

# build a label table from the authoritative tissue_order so colors and labels are consistent
label_tbl <- tibble(
  tissue = tissue_order,
  system = system_map[tissue_order]
) %>%
  group_by(system) %>%
  mutate(
    # alternate colors restarting within each system
    col = rep(c("#000000", "#6e6e6e"), length.out = n())
  ) %>%
  ungroup()

# join colors back to plotting data
df2 <- df2 %>%
  left_join(label_tbl, by = c("tissue", "system")) %>%
  mutate(
    tissue = factor(tissue, levels = tissue_order),
    system = factor(system, levels = system_order)
  ) %>% group_by(phen) %>% mutate(n = n(), a = (n == 1),
                                  a = factor(a, levels = c(FALSE, TRUE), labels = c("Group", "Single")
                                  ))


# Create colored axis labels (HTML) using ggtext::element_markdown for rendering
axis_labels <- label_tbl %>%
  mutate(
    label_html = paste0("<span style='color:", col, "'>", tissue, "</span>")
  )

# create a named vector mapping tissue -> label_html (used by scale_y_discrete)
axis_label_vector <- setNames(axis_labels$label_html, axis_labels$tissue)



p_tissue <- df2 %>%
  ggplot(aes(
    x = acrophase_24hfreq,
    y = tissue,
    label = phen,
    color = col,
    shape = a
  )) +

  # If your light_band/night_band use the same x-range for all facets and ymins/ymax
  # are provided as factor names or numeric positions, ensure they match your coordinate system.
  # The examples below assume light_band/night_band are already prepared with ymin/ymax on the
  # same scale as the y factor (if not, see the numeric-mapping helper above).
  geom_rect(data = light_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow",
            alpha = 0.3,
            inherit.aes = FALSE) +

  geom_rect(data = night_band,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightblue",
            alpha = 0.2,
            inherit.aes = FALSE) +

  geom_point(size = 3, alpha = 0.9) +

  ggrepel::geom_label_repel(
    size = 2,
    max.overlaps = 70,
    box.padding = 0.1,
    label.padding = 0.1,
    alpha = 0.75
  ) +

  scale_x_continuous(
    limits = c(0, 24),
    breaks = 0:23,
    expand = c(0, 0)
  ) +

  # Force ggplot to display tissues in EXACT order defined in tissue_order,
  # top-to-bottom. (We reverse because discrete y is drawn bottom->top.)
  scale_y_discrete(
    position = "right",
    #limits = rev(tissue_order),
    labels = axis_label_vector
  ) +

  scale_color_identity() +


  facet_grid(
    rows = vars(system),
    scales = "free",
    space = "free",
    switch = "y"
  ) +

  labs(
    x = "Acrophase",
    y = "GTEx tissue",
    title = "H", shape = "Enrichment"
  ) +

  theme_classic() +
  theme(
    plot.title = element_text(face = "bold"),

    panel.border = element_rect(colour = "black", fill = NA),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),

    legend.position = "bottom",
    legend.box.spacing = unit(2, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),

    # HTML colored axis labels:
    axis.text.y.right = ggtext::element_markdown(),
    axis.text.y = ggtext::element_markdown(),

    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.right = element_text(face = "bold")
  )

p_tissue


