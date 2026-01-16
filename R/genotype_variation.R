
target_prots <- c(
  "MYOC","HYAL1","EFNA1","TNR","FAS","RELT","CD276",
  "SPON2","SPINK5","GDF15","PGF","ANGPTL1",
  "KLK12","LGALS1","HS3ST3B1", "KLK13"
)


rsid <- data.table::fread("CA_cojo_combined_rs.txt")

olink_t <- readRDS("/mnt/project/biomarkers_res/res_olink_tech_14panels.rds")


olink_t %>%
  rename_with(tolower) %>%
  select(eid, any_of(tolower(target_prots))) %>%
  left_join(time) %>%
  filter(time_day > 9 & time_day < 20) %>%
  left_join(rsid, by = c("eid" = "FID")) %>%

  pivot_longer(2:17) %>%
  ggplot(aes(x = time_day, y = value, color = name)) +
  geom_smooth() +
  theme_classic(base_size = 16)





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



