

batch <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_batch_number.dat")
panel <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_assay.dat") %>%
  mutate(prot = tolower(Assay), Panel = str_remove(Panel, " ")) %>% select(prot, Panel)

panels_14 <- panel$prot[panel$Panel %in% c("Cardiometabolic", "Inflammation", "Neurology", "Oncology")]

prots <- data.table::fread("/mnt/project/olink_instance_0.csv") %>% as_tibble() %>%
  select(contains(panels_14))

row_na_counts <- rowSums(is.na(prots))
prots_f <- prots[row_na_counts < ncol(prots)/3,]

col_na_counts <- colSums(is.na(prots_f))/nrow(prots_f)

prots_qc <- prots_f[,col_na_counts < 0.1]

#toupper(names(which(col_na_counts > 0.1)))
# [1] "ADAM12"    "ADAM9"     "ADAMTS1"   "ADAMTS4"   "ADAMTSL2"  "ADAMTSL4"  "ADAMTSL5"  "TADA3"     "ACADM"
# [10] "AIF1L"     "AKR1B10"   "APPL2"     "NAPRT"     "ARNTL"     "FARSA"     "C2ORF69"   "CCDC28A"   "CDC25A"
# [19] "CDC26"     "CLEC2L"    "DCDC2C"    "DNAJC21"   "DOC2B"     "ERC2"      "GIPC2"     "HDDC2"     "IGLC2"
# [28] "MAMDC2"    "MUC2"      "MYBPC2"    "NPC2"      "PLXDC2"    "SLC28A1"   "TSC22D1"   "VWC2L"     "ABCA2"
# [37] "CD164L2"   "CD226"     "CD5L"      "PCDH12"    "CDH22"     "CDH23"     "CEACAM16"  "CEACAM18"  "CEACAM19"
# [46] "SCPEP1"    "CSF1R"     "CSF3R"     "CTSS"      "DNAJB14"   "EGFLAM"    "MEGF11"    "VEGFB"     "CSF2RB"
# [55] "FGF3"      "MTIF3"     "RNF31"     "TRAF3"     "TRAF3IP2"  "FGF7"      "ZNF75D"    "C1QTNF9"   "FGF9"
# [64] "IGSF9"     "AFAP1"     "MFAP3L"    "MFAP4"     "TFAP2A"    "FGF20"     "FSTL1"     "C1GALT1C1" "GALNT5"
# [73] "ITGAL"     "LGALS3BP"  "GLRX5"     "GMPR2"     "NGRN"      "HDGFL2"    "HGFAC"     "IL13RA2"   "IL20RB"
# [82] "IL21R"     "IL22"      "IL25"      "IL2RG"     "SERPING1"  "SPRING1"   "KLK15"     "LATS1"     "LDLRAP1"
# [91] "RILPL2"    "IGHMBP2"   "LAMB1"     "MBL2"      "NUMB"      "MMP15"     "MSLNL"     "NGFR"      "NPM1"
# [100] "SPAG1"     "PCOLCE"    "PROCR"     "PRSS22"    "TPSG1"     "ERP29"     "LRP2"      "LRP2BP"    "PGLYRP2"
# [109] "SELENOP"   "SSBP1"     "TACSTD2"   "C1QTNF5"   "C1QTNF6"   "TNFAIP2"   "TNFAIP8L2" "TNFRSF17"  "TNFSF8"
# [118] "TP53BP1"   "TP53I3"    "TRIM58"    "TSPAN15"   "WASHC3"    "WASL"

saveRDS(prots_qc, "olink_instance_0_QC_14_panels.rds")
