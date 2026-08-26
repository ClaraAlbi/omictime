library(tidyverse)

assay <- data.table::fread("data/olink_assay.dat") %>%
  mutate(Assay = toupper(Assay),
         Assay = str_replace(Assay, "-", "_"))

df_effects <- readRDS("data/combined_effects.rds") %>% mutate(phen = toupper(phen)) %>%
  left_join(assay, by = c("phen" = "Assay")) %>%
  mutate(pval_h = p.adjust(pvalue_h)) %>%
  filter(pval_h < 0.05) %>%
  filter(type_clean == "Proteins")

write(df_effects$UniProt, "data/list_rhythmic.txt")

df_effects %>%
  mutate(Panel = str_remove(Panel, " II")) %>%
  group_by(Panel) %>% count()

go <- data.table::fread("data/explore_ukb (1).csv") %>%
  dplyr::rename(UniProt = `UniProt ID`) %>%
  mutate(
    go_bp = str_trim(`Biological Process`))  %>%
  separate_rows(go_bp, sep = ",\\s*") %>%
  filter(go_bp != "")

background <- go %>%
  #filter(UniProt %in% assay$UniProt) %>%
  group_by(go_bp) %>%
  #filter(str_detect(go_bp, "circadi")) %>%
  summarise(l = list(Gene),  n = n(), type = "bg")

rhythmic <- go %>%
  #filter(str_detect(go_bp, "circadi")) %>%
  filter(UniProt %in% df_effects$UniProt) %>% group_by(go_bp) %>%
  summarise(l = list(Gene),  n = n(), type = "r")


#####
library(clusterProfiler)
library(org.Hs.eg.db)

#BiocManager::install("org.Hs.eg.db")
bg_entrez <- bitr(assay$UniProt,
                  fromType = "UNIPROT",
                  toType   = "ENTREZID",
                  OrgDb    = org.Hs.eg.db)

sig_entrez <- bitr(df_effects$UniProt,
                   fromType = "UNIPROT",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

ego <- enrichGO(
  gene          = sig_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",        # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_mf <- enrichGO(
  gene          = sig_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",        # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ekegg <- enrichKEGG(
  gene          = sig_entrez$ENTREZID,
  universe      = bg_entrez$ENTREZID,
  organism      = "hsa",
  keyType       = "kegg",
  pvalueCutoff  = 0.05
)

all_top <- bind_rows(ego@result %>%
  as.data.frame(), ekegg@result %>%
    as.data.frame()) %>%
  arrange(pvalue)

write_csv(all_top, "data/pathway_enrichment.csv")

