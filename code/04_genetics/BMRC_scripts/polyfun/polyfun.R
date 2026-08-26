source("config/paths.R")


files <- list.files(bmrc_path("results"), "jma", full.names = T)

cojo_CA <- bind_rows(lapply(files[str_detect(files, "cojo_res_v4")], data.table::fread)) %>%
  arrange(Chr) %>%
  mutate(id = row_number())

finemap <- data.table::fread("results/polyfun/polyfun_CA.txt.gz")

f <- finemap %>%
  group_by(CREDIBLE_SET) %>%
  mutate(sum_PIP = sum(PIP), n = n()) %>%
  filter(sum_PIP >= 0.95)  %>%
  arrange(n) %>%
  filter(n < 10)

data.table::fwrite(f, "OUTPUTS/TX_credible_set_FINEMAP.txt", sep = " ", row.names = F)



ldsc <- data.table::fread(bmrc_path("results/polyfun/Sldsc.results"))


ann <- data.table::fread("credible_set_annotations.txt") %>%
  group_by(`#Uploaded_variation`) %>% 
  summarise(Consequence = list(Consequence), Gene = list(SYMBOL), type = list(BIOTYPE))

finemap2 <- finemap %>%
  mutate(.finemap_row = row_number())   # unique row id to re-join later

finemap2 <- finemap %>%
  rename(Chr = CHR, finemap_BP = BP, finemap_SNP = SNP) %>%
  mutate(.finemap_row = row_number())   # unique row id to re-join later

cojo2 <- cojo_CA %>%
  rename(Chr = Chr, cojo_BP = bp, cojo_id = id)

# 1) Find all candidate matches by chromosome, compute distance, filter to 10kb
matches <- finemap2 %>%
  inner_join(cojo2, by = c("Chr")) %>%
  mutate(dist = abs(finemap_BP - cojo_BP)) %>%
  filter(dist <= 5000000)

nearest_matches <- matches %>%
  group_by(.finemap_row) %>%
  slice_min(dist, with_ties = FALSE) %>%
  ungroup() %>%
  select(.finemap_row, cojo_id, cojo_BP, dist)

# 3) Left-join back to the original finemap to keep all finemap rows (NA for no match)
finemap_mapped <- finemap2 %>%
  left_join(nearest_matches, by = ".finemap_row") %>%
  select(-.finemap_row)

df_cs <- finemap_mapped %>%
  group_by(cojo_id, CREDIBLE_SET) %>%                
  arrange(desc(PIP), .by_group = TRUE) %>%
  mutate(
    high_pip_gt_0_95 = PIP > 0.95,
    pip_cumsum = cumsum(PIP),
    in_95_credible_set = pip_cumsum >= 0.95 ,
    n = n()
  ) %>%
  ungroup() %>%
  filter(n < 6)

df_cs %>% filter(in_95_credible_set) %>% pull(SNP) %>% write("credible_sets.txt")

df_cs %>%
  group_by(cojo_id, CREDIBLE_SET) %>% count()

w_ann <- df_cs %>%
  left_join(ann, by = c("finemap_SNP" = "#Uploaded_variation")) 


p <- cojo_CA %>%
  left_join(ann, by = c("SNP" = "#Uploaded_variation")) 
  
