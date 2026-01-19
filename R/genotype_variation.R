
target_prots <- tribble(~prot, ~snp,
  "MYOC", "rs2234926",
  "HYAL1", "rs116482870",
  "EFNA1", "rs4390169",
  "TNR", "rs7515728",
  "FAS", "rs3781202",
  "RELT", "rs12099129",
  "CD276", "rs74026250",
  "SPON2", "rs7678661",
  "SPINK5", "rs2400477",
  "GDF15", "rs1054221",
  "PGF", 'rs6574205',
  "ANGPTL1", "rs1018669",
  "KLK12", "rs2569459",
  "KLK13", "rs2569459",
  "LGALS1", "rs62236671",
  "HS3ST3B1","rs62053895"
)

rsid <- data.table::fread("/mnt/project/CA_cojo_combined_rs.txt")

olink_t <- readRDS("/mnt/project/biomarkers_res/res_olink_tech_14panels.rds")


olink_long <- olink_t %>%
  rename_with(tolower) %>%
  select(eid, any_of(tolower(target_prots$prot))) %>%
  pivot_longer(
    -eid,
    names_to = "prot",
    values_to = "q"
  )

rs_long <- rsid %>%
  rename(eid = FID) %>%
  select(eid, contains(target_prots$snp)) %>%
  pivot_longer(
    -eid,
    names_to = "snp",
    values_to = "genotype"
  ) %>%
  filter(!is.na(genotype)) %>%
  mutate(snp = sub("_.*$", "", snp))


paired_data <- target_prots %>% mutate(prot = tolower(prot)) %>%
  left_join(olink_long, by = "prot") %>%
  left_join(rs_long, by = c("eid", "snp"))



paired_data %>%
  slice(1:1000000) %>%
  left_join(time %>% select(eid, time_day)) %>%
  filter(time_day > 9 & time_day < 20) %>%
  filter(!is.na(genotype)) %>%
  ggplot(aes(x = time_day, y = q, color = as.factor(genotype))) +
  geom_smooth() +
  facet_wrap(~prot+snp,scales = "free") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")


geno_mat <- paired_data %>%
  select(eid, snp, genotype) %>%
  distinct() %>%
  pivot_wider(
    names_from = snp,
    values_from = genotype
  ) %>%
  column_to_rownames("eid") %>%
  as.matrix()

geno_mat[is.na(geno_mat)] <- apply(geno_mat, 2, mean, na.rm = TRUE)


geno_scaled <- scale(geno_mat)

pca <- prcomp(geno_scaled, center = TRUE, scale. = FALSE)

pca_df <- as.data.frame(pca$x[, 1:2])
pca_df$eid <- rownames(pca_df)

ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(alpha = 0.6) +
  theme_classic(base_size = 14)



chr1 <- data.table::fread("/mnt/project/snps/subset_cojo_pQTL_chr1.raw") %>%
  filter(IID %in% olink$eid) %>%
  mutate(rs57742792_T = round(rs57742792_T, 0))

chr4 <- data.table::fread("subset_cojo_pQTL_chr4.raw") %>%
  filter(IID %in% olink$eid) %>%
  mutate(rs28856334_G = round(rs28856334_G, 0))

chr5 <- data.table::fread("subset_cojo_pQTL_chr5.raw") %>%
  filter(IID %in% olink$eid) %>%
  mutate(rs28856334_G = round(rs28856334_G, 0))

chr11 <- data.table::fread("subset_cojo_pQTLchr11.raw") %>%
  filter(IID %in% olink$eid) %>%
  mutate(rs12099129_T = round(rs12099129_T, 0))



df %>%
  left_join(chr1, by = c("eid" = "FID")) %>%
  filter(!is.na(rs57742792_T)) %>%
  ggplot(aes(x = factor(rs57742792_T), y = res)) + geom_boxplot() + theme_minimal()


df %>%
  left_join(olink) %>%
  left_join(chr1, by = c("eid" = "FID")) %>%
  mutate(rs7515728_T = as.factor(round(rs7515728_T, 0))) %>%
  filter(!is.na(rs7515728_T)) %>%
  filter(time_day > 10 & time_day < 13) %>%
  ggplot(aes(x = rs7515728_T, y = res, color = rs7515728_T)) + geom_boxplot() + theme_minimal()


df %>%
  left_join(olink) %>%
  filter(!is.na(rs7515728_T)) %>%
  ggplot(aes(x = rs7515728_T, y = tnr, color = rs7515728_T)) +
  geom_boxplot() + theme_minimal() +
  geom_

#summary(lm(res ~ rs7515728_T, data = dd1 %>% filter(t)))

# CHR 4

df %>%
  left_join(olink) %>%
  left_join(chr4, by = c("eid" = "FID")) %>%
  mutate(rs28856334_G = as.factor(rs28856334_G)) %>%
  filter(!is.na(rs28856334_G)) %>%
  ggplot(aes(x = time_day, y = res, color = rs28856334_G)) + geom_smooth() + theme_minimal()

df %>%
  left_join(olink) %>%
  left_join(chr4, by = c("eid" = "FID")) %>%
  mutate(rs28856334_G = as.factor(rs28856334_G)) %>%
  filter(!is.na(rs28856334_G)) %>%
  ggplot(aes(x = rs28856334_G, y = res, color = rs28856334_G)) + geom_boxplot() + theme_minimal()


df %>%
  left_join(olink) %>%
  left_join(chr4, by = c("eid" = "FID")) %>%
  mutate(rs28856334_G = as.factor(rs28856334_G)) %>%
  filter(!is.na(rs28856334_G)) %>%
  ggplot(aes(x = time_day, y = spon2, color = rs28856334_G)) + geom_smooth() + theme_minimal()


df %>%
  left_join(olink) %>%
  left_join(chr4, by = c("eid" = "FID")) %>%
  mutate(rs28856334_G = as.factor(rs28856334_G)) %>%
  filter(!is.na(rs28856334_G)) %>%
  ggplot(aes(x = rs28856334_G, y = spon2, color = rs28856334_G)) + geom_boxplot() + theme_minimal()



# CHR 11


# CHR 4

df %>%
  left_join(olink) %>%
  left_join(chr11, by = c("eid" = "FID")) %>%
  mutate(rs12099129_T = as.factor(rs12099129_T)) %>%
  filter(!is.na(rs12099129_T)) %>%
  ggplot(aes(x = time_day, y = res, color = rs12099129_T)) + geom_smooth() + theme_minimal()

df %>%
  left_join(olink) %>%
  left_join(chr11, by = c("eid" = "FID")) %>%
  mutate(rs12099129_T = as.factor(rs12099129_T)) %>%
  filter(!is.na(rs12099129_T)) %>%
  ggplot(aes(x = rs12099129_T, y = res, color = rs12099129_T)) + geom_boxplot() + theme_minimal()


df %>%
  left_join(olink) %>%
  left_join(chr11, by = c("eid" = "FID")) %>%
  mutate(rs12099129_T = as.factor(rs12099129_T)) %>%
  filter(!is.na(rs12099129_T)) %>%
  ggplot(aes(x = time_day, y = relt, color = rs12099129_T)) + geom_smooth() + theme_minimal()


df %>%
  left_join(olink) %>%
  left_join(chr11, by = c("eid" = "FID")) %>%
  mutate(rs12099129_T = as.factor(rs12099129_T)) %>%
  filter(!is.na(rs12099129_T)) %>%
  ggplot(aes(x = rs12099129_T, y = relt, color = rs12099129_T)) + geom_boxplot() + theme_minimal()



