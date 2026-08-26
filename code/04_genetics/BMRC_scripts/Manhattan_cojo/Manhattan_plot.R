#!/usr/bin/env Rscript
# Minimal Manhattan plot (no highlights, no labels)


# --- User-editable settings ---
output_prefix   <- "gwas_manhattan_v5"
png_w <- 10; png_h <- 4; dpi <- 300
keep_p_threshold <- 1e-6    # always keep variants with p < this
sample_fraction   <- 0.02   # fraction of remaining variants to keep
per_chr_min_keep  <- 500    # guarantee at least this many points per chr (if possible)

gene_list <- data.table::fread("cojo/CAannotate.variant_function", header = F) %>%
  separate(V3,
           into = c("chr", "POS", "pos2", "A1", "A2"),
           sep = " ") %>%
  select(-pos2) %>%
  mutate(CHR = str_remove(chr, "chr"), POS = as.numeric(POS)) %>%
  select(V1, V2, CHR, POS)


# -------------------------------

library(data.table); library(dplyr); library(ggplot2); library(scales); library(ggrepel)

# read
gwas <- data.table::fread("data/GWAS_res_v5/CA_v5.fastGWA") %>%
  mutate(CHR = as.character(CHR)) %>%
  bind_rows(data.table::fread("data/GWAS_res/res_chrX.fastGWA") %>% select(-BETA_F, -BETA_M, -SE_M, -SE_F)) %>%
  mutate(
    CHR_LABEL = as.character(CHR),
    CHR_NUM = if_else(toupper(CHR_LABEL) == "X", 23L, as.integer(CHR_LABEL))
  ) %>%
  filter(!is.na(CHR_NUM), !is.na(POS), !is.na(P)) %>%
  mutate(IN_GENE_LIST = paste(CHR, POS) %in%
           paste(gene_list$CHR, gene_list$POS))

keep_mask <- gwas$P < keep_p_threshold
gwas$TO_KEEP <- keep_mask

gwas[gwas$CHR == "X",]



# downsample per chromosome (simplest form)
if (sample_fraction < 1) {
  gwas <- gwas %>%
    group_by(CHR_NUM) %>%
    mutate(
      KEEP_ALWAYS = TO_KEEP | IN_GENE_LIST,   # <- key change
      KEEP_RANDOM = !KEEP_ALWAYS & (runif(n()) < sample_fraction),
      TO_KEEP = KEEP_ALWAYS | KEEP_RANDOM
    ) %>%
    ungroup()
} else {
  gwas$TO_KEEP <- TRUE
}

# create plotting dataset
plot_df <- gwas %>%
  filter(TO_KEEP) %>%
  mutate(MINUSP = -log10(pmax(P, 1e-300)))

# compute chromosome offsets for cumulative BP
chr_info <- plot_df %>%
  group_by(CHR_NUM) %>%
  summarise(maxpos = max(POS)) %>%
  arrange(CHR_NUM) %>%
  mutate(
    offset = cumsum(lag(maxpos + 1, default = 0)),
    center = offset + maxpos / 2
  )

plot_df <- plot_df %>%
  left_join(chr_info, by = "CHR_NUM") %>%
  mutate(
    BPcum = POS + offset,
    CHR_COL = factor(CHR_NUM %% 2)
  )

gene_list2 <- gene_list %>%
  left_join(plot_df) %>%
  mutate(
    MINUSP = -log10(pmax(P, 1e-300)),
    GENE_CLEAN = case_when(
      V1 %in% c("intronic", "exonic") ~ V2,
      
      V1 %in% c("intergenic", "upstream", "upstream;downstream", "UTR3") ~
        V2 %>%
        str_remove_all("\\s*\\(.*?\\)") %>%  # remove (dist=xxxx)
        str_replace_all(",\\s*", ",") %>%     # clean spacing
        str_remove("[,;].*") %>%
        paste0("*"),
      
      TRUE ~ V2
    )
  )

gene_list2[20,]$MINUSP <- 12
gene_list2[11,]$MINUSP <- 19.3
gene_list2[12,]$MINUSP <- NA
gene_list2[5,]$MINUSP <- NA

gene_list2[6,]$GENE_CLEAN <- "ANGPTL1*"


# Manhattan plot
p <- ggplot(plot_df, aes(x = BPcum, y = MINUSP)) +
  geom_point(aes(color = CHR_COL), size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("#009FB7", "#7A6263"), guide = "none") +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed") +
  scale_x_continuous(
    breaks = chr_info$center[c(1:12, 14, 16, 18, 20, 22, 23)],
    labels = c(1:12, 14, 16, 18, 20, 22, "X"),
    expand = c(0.01, 0.01)
  ) +
  geom_label_repel(min.segment.length = 0,
    data = gene_list2,
    aes(x = BPcum, label = GENE_CLEAN, color = CHR_COL, y = MINUSP),
    size = 2.5
  ) +
  labs(x = "Chromosome", y = expression(-log[10](italic(P)))) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size  = 8),
    axis.text.y = element_text(size  = 8),
    axis.title = element_text(size  = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), plot.margin = margin(1,1,1,1)
  )

ggsave(paste0("OUTPUTS/",output_prefix, ".png"), p, width = png_w, height = png_h, dpi = dpi)
message("Saved plot: ", paste0(output_prefix, ".png"))
