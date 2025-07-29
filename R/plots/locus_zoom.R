# 1. Install/load
library(tidyverse)

if (!requireNamespace("Gviz", quietly=TRUE)) BiocManager::install("Gviz")
if (!requireNamespace("biomaRt", quietly=TRUE)) BiocManager::install("biomaRt")
library(Gviz); library(biomaRt)

rsids <- data.table::fread("~/Downloads/snp.info") %>% as_tibble()
gwas <- data.table::fread("~/Downloads/gwas_res.res.glm.linear.gz") %>% as_tibble() %>%
  inner_join(rsids %>% dplyr::select(Chrom, PhysPos, ID), by = c("CHROM" = "Chrom", "POS" = "PhysPos")) %>%
  dplyr::rename(
    chrom = CHROM,
    pos   = POS,
    rsid  = ID.y,
    pval  = P
  ) %>%
  dplyr::mutate(
    chrom = as.integer(chrom),
    logp  = -log10(pval)
  )

if (!requireNamespace("ggbio", quietly=TRUE)) BiocManager::install("ggbio")
if (!requireNamespace("GenomicFeatures", quietly=TRUE)) BiocManager::install("GenomicFeatures")
library(ggbio); library(GenomicFeatures)
library(ggplot2); library(patchwork)


txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.107.gtf.gz", format="gtf")
