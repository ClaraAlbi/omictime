library(dplyr)
library(tidyr)
library(stringr)

### PROTEIN QC CODE

prots <- data.table::fread("/mnt/project/olink_instance_0.csv") %>% as_tibble()

row_na_counts <- rowSums(is.na(prots))
prots_f <- prots[row_na_counts < ncol(prots)/3,]

col_na_counts <- colSums(is.na(prots_f))/nrow(prots_f)

prots_qc <- prots_f[,col_na_counts < 0.1]
# toupper(names(which(col_na_counts > 0.1)))
# [1] "AMY2B"   "CST1"    "CTSS"    "GLIPR1"  "NPM1"    "PCOLCE"  "TACSTD2"

# Remove batch 0 (too small)

olink_times <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_processing_start_date.dat")
batch <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_batch_number.dat")
panel <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_assay.dat") %>%
  mutate(prot = tolower(Assay), Panel = str_remove(Panel, " ")) %>% select(prot, Panel)

metadata <- data.table::fread("/mnt/project/olink_instance_0_meta.csv") %>%
  rename(num_prots = p30900_i0,
         PlateID = p30901_i0,
         well_id = p30902_i0,
         ppp_sel = p30903_i0) %>%
  filter(!is.na(PlateID)) %>%
  left_join(batch) %>%
  left_join(olink_times %>% pivot_wider(names_from = Panel, values_from = Processing_StartDate)) %>%
  mutate(ppp_sel = case_when(ppp_sel == "Yes" ~ 1,
                             TRUE ~ 0))

table(metadata$Batch)
# 0     1     2     3     4     5     6     7
# 69  1691  5679 13048  8484  8299  8230  7516

prots_out <- prots_qc %>%
  left_join(metadata) %>%
  filter(!Batch %in% c(0)) %>%
  select(everything(), -any_of(colnames(metadata)[-1]))
table(prots_out$Batch)
# 1     2     3     4     5     6     7
# 1691  5679 13048  8484  8299  8230  7516

saveRDS(prots_out, "olink_instance_0_QC.rds")



### use only Panels 1-4

panels_14 <- panel$prot[panel$Panel %in% c("Cardiometabolic", "Inflammation", "Neurology", "Oncology")]

saveRDS(panels_14, "olink_panels_1to4.rds")

prots_14 <- data.table::fread("/mnt/project/olink_instance_0.csv") %>% as_tibble() %>%
  select(eid, any_of(panels_14))

col_na_counts <- colSums(is.na(prots_14))/nrow(prots_14)

prots_14_qc <- prots_14[,col_na_counts < 0.1]
# toupper(names(which(col_na_counts > 0.1)))
# [1] "CTSS"    "NPM1"    "PCOLCE"  "TACSTD2"

row_na_counts <- rowSums(is.na(prots_14_qc))
prots_14_qc <- prots_14_qc[row_na_counts < ncol(prots_14_qc)/3,]
dim(prots_14_qc)
prots_14_qc_out <- prots_14_qc %>%
  left_join(metadata) %>%
  filter(!Batch %in% c(0)) %>%
  select(everything(), -any_of(colnames(metadata)[-1]))
table(prots_14_qc_out$Batch)
# 1     2     3     4     5     6     7
# 1671  5633 12834  8116  8190  8147  7444

saveRDS(prots_14_qc_out, "olink_instance_0_QC_panels14.rds")
