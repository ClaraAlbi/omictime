

batch <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_batch_number.dat")
panel <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_assay.dat") %>%
  mutate(prot = tolower(Assay), Panel = str_remove(Panel, " ")) %>% select(prot, Panel)

meta <- data.table::fread("/mnt/project/olink_instance_0_meta.csv") %>% filter(!is.na(p30900_i0))

panels_14 <- panel$prot[panel$Panel %in% c("Cardiometabolic", "Inflammation", "Neurology", "Oncology")]

panels_all <- panel

prots <- data.table::fread("/mnt/project/olink_instance_0.csv") %>% as_tibble()

prots_2 <- prots %>%
  select(eid, any_of(panels_14))

row_na_counts <- rowSums(is.na(prots_2))

ids_s <- which(row_na_counts > 1536)

prots_f <- prots_2[row_na_counts < ncol(prots_2)/3,]

col_na_counts <- colSums(is.na(prots_f))/nrow(prots_f)

prots_qc <- prots_f[,col_na_counts < 0.1]

#toupper(names(which(col_na_counts > 0.1)))
#[1] "CTSS"    "NPM1"    "PCOLCE"  "TACSTD2"

saveRDS(prots_qc, "olink_instance_0_QC_14_panels.rds")
