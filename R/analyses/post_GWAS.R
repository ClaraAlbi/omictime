library(tidyverse)

gwas <- data.table::fread("~/Downloads/gwas_res.res.glm.linear.gz") %>%
  filter(TEST == "ADD")

circadian_genes <- data.table::fread("~/Downloads/merged_final_collected_CR_geneset_space_removed_annotated_with_gene_name_add_Science_46tissues_baboon_paper_mouse_SCN_related.txt") %>%
  filter(label == "GOCRpath_human_NA")

library(biomaRt)
ensembl37 <- useEnsembl(biomart = "genes",
                        dataset = "hsapiens_gene_ensembl",
                        version = 75)
genes <- circadian_genes$Gene.name
res <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
             filters = "hgnc_symbol", values = genes, mart = ensembl37)

gene_seq <- data.table::fread("~/genes.gtf")

gene_seq$ensembl_id <- sub('.*gene_id "([^"]+)".*', '\\1', gene_seq$V9)
gene_seq$gene_name <- sub('.*gene_name "([^"]+)".*', '\\1', gene_seq$V9)

circadian_coords <- circadian_genes %>%
  left_join(gene_seq, by = c("Human.gene.stable.ID" = "ensembl_id")) %>%
  mutate(CHR = as.numeric(V1))


result <- data.frame()
for (i in 1:nrow(sig)) {
  snp_id  <- sig$ID[i]

  if (!snp_id %in% result$ID) {
    snp_chr <- sig$`CHROM`[i]
    snp_pos <- sig$POS[i]

    genes_on_chr <- gene_seq[as.numeric(gene_seq$V1) == snp_chr, ]
    genes_on_chr$distance <- pmin(abs(snp_pos - genes_on_chr$V4), abs(snp_pos - genes_on_chr$V5))
    closest_gene <- genes_on_chr[which.min(genes_on_chr$distance), ]

    result <- rbind(result, data.frame(
      snp = snp_id,
      chr = snp_chr,
      pos = snp_pos,
      gene_id = closest_gene$gene_name
    ))
  }
}

a <- circadian_genes %>%
  inner_join(result, by = c( "Gene.name"= "gene_id"))

i <- 3902473
snp_id  <- time$ID[i]
snp_chr <- time$`#CHROM`[i]
snp_pos <- time$POS[i]

genes_on_chr <- gene_seq[as.numeric(gene_seq$V1) == snp_chr, ]
genes_on_chr$distance <- pmin(abs(snp_pos - genes_on_chr$V4), abs(snp_pos - genes_on_chr$V5))
closest_gene <- genes_on_chr[which.min(genes_on_chr$distance), ]



########################################################################

sig <- gwas %>%
  filter(P < 1e-4)


# Compute cumulative BP position
chr_info <- sig %>%
  group_by(`CHROM`) %>%
  summarize(chr_len = max(POS)) %>%
  arrange(`CHROM`) %>%
  mutate(chr_start = lag(cumsum(as.numeric(chr_len)), default = 0))

sig <- sig %>%
  left_join(chr_info, by = "CHROM") %>%
  mutate(BP_cum = as.numeric(POS) + chr_start)

# Chromosome label positions
axis_df <- chr_info %>%
  mutate(center = chr_start + chr_len / 2)

# Plot
res <- sig %>%
  left_join(result, by = c("CHROM" = "chr", "POS" = "pos"))

res %>%
  ggplot(aes(x = BP_cum, y = -log10(P))) +
  geom_point(size = 1.2, alpha = 0.75) +
  geom_text(
    data            = subset(res, P < 1e-7),
    aes(x = BP_cum, y = -log10(P), label = gene_id),
    vjust           = -0.5,
    hjust = -0.2,
    size            = 2,
    check_overlap   = TRUE
  ) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$`CHROM`) +
  labs(x = "Chromosome", y = "-log10(p)") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )
