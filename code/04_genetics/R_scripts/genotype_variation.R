source("config/paths.R")


library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(purrr)
library(ggplot2)
library(broom)


target_prots2 <- tribble(
  ~prot, ~snp,
  "EFNA1", "rs4390169",
  "MYOC", "rs2234926",
  "TNR", "rs7515728",
  "C1orf220", "rs2761462",
  "ANGPTL1", "rs1018669",
  "ANGPTL1", "rs34036521",
  "FAS", "rs3781202",
  "RELT", "rs12099129",
  "PGF", "rs6574205",
  "HS3ST3A1", "rs62053895",
  "HS3ST3B1", "rs62053895",
  "GDF15", "rs3787023",
  "GDF15", "rs1054221",
  "KLK12", "rs3760744",
  "NLRP12", "rs34436714",
  "LGALS1", "rs62236671",
  "SYN2", "rs184262",
  "DOCK3", "rs115981680",
  "SPON2", "rs7678661",
  "SPINK5", "rs2400477",
  "COLEC10", "rs2450067"
) %>%
  mutate(
    prot_lower = str_to_lower(prot)
  )

olink_t <- readRDS(project_path("biomarkers_res/res_olink_tech_14panels.rds")) %>%
  select(eid, any_of(target_prots2$prot_lower))

files <- list.files(project_path("raw_snps/"),
                    pattern = "\\.raw$",
                    full.names = TRUE)

dt_list <- lapply(files, fread)

combined <- dt_list[[1]]

for (i in 2:length(dt_list)) {
  combined[, (names(dt_list[[i]])[-(1:6)]) :=
             dt_list[[i]][, -(1:6), with = FALSE]]
}

#saveRDS(combined %>% select(eid = FID, starts_with("rs")), "top_snps_CA.rds")

geno_cols <- colnames(combined)[-(1:6)]

geno_map <- tibble(
  geno_col = geno_cols,
  snp = sub("_.*", "", geno_cols)
)


pair_map <- geno_map %>%
  inner_join(target_prots2, by = "snp") %>% filter(prot_lower %in% colnames(olink_t)) #%>%
  filter(snp != "rs3787023")


res <- readRDS(project_path("olink_int_panels14.rds")) %>%
  select(eid, time_day, pred_mean)
res$res <- residuals(lm(pred_mean ~ time_day, data = res))

analysis_df <- combined %>%
  rename(eid = IID) %>%
  inner_join(olink_t, by = "eid") %>%
  inner_join(res) %>%
  mutate(q_res = ntile(res, 5))


plot_data <- pair_map %>%
  mutate(
    data = purrr::map2(geno_col, prot_lower, ~
                         analysis_df %>%
                         select(
                           eid,
                           time_day,
                           q_res,
                           res,
                           genotype = all_of(.x),
                           protein  = all_of(.y)
                         )
    )
  ) %>%
  tidyr::unnest(data) %>%
  mutate(genotype = factor(genotype), q_res = factor(q_res))



p1 <- plot_data %>%
  filter(!is.na(genotype)) %>%
  #slice(1:5000) %>%
  ggplot(aes(x = time_day, y = protein, color = genotype)) +
  geom_smooth() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Time of day", y = "Protein levels", color = "Genotype") +
  facet_wrap(~prot+geno_col) + theme_minimal()

ggsave("plots/FS_rsid.png", p1, width = 7, height = 5)

plot_data %>% filter(prot == "SPINK5") %>%
  lm(res ~ genotype, data = .) %>% summary()

### plot all
data_all <- analysis_df %>%
  mutate(t = round(time_day, 0)) %>%
  pivot_longer(starts_with("rs"), names_to = "geno_col", values_to = "genotype") %>%
  filter(!is.na(genotype)) %>%

  pivot_longer(contains(target_prots2$prot_lower), names_to = "protein", values_to = "prot_values") %>%
  group_by(t, protein, genotype, geno_col) %>% summarise(m_p = mean(prot_values, na.rm = T),
                                                    sd_p = sd(prot_values, na.rm = T),
                                                    n = n()) %>%
  filter(!is.na(m_p)) %>% ungroup()


data_all %>% group_by(t, geno_col, genotype) %>% count()

p_all <- data_all %>%
  filter(!geno_col %in% pair_map$geno_col) %>%
  #filter(geno_col != "rs3787023_A") %>%
  filter(n > 10) %>%
  left_join(geno_map) %>%
  left_join(target_prots2) %>%
  mutate(protein = toupper(protein)) %>%
  ggplot(aes(x = t, y = m_p, color = as.factor(genotype))) +
  geom_point(size = 1) +
  facet_grid(rows = vars(protein), cols = vars(snp, prot)) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Time of day", y = "Protein levels", color = "Genotype") +
  theme_minimal()+
  theme(legend.position = "bottom")


ggsave("plots/FS_rsid_all.png", p_all, width = 8, height = 10)


## FORMAL TEST?

v <- setdiff(target_prots2$snp, pair_map$snp)

d <- analysis_df %>%
  select(eid, contains(v), any_of(target_prots2$prot_lower)) %>%
  pivot_longer(starts_with("rs"), names_to = "geno_col", values_to = "genotype") %>%
  pivot_longer(-c(eid, geno_col, genotype)) %>%
  group_by(geno_col, name) %>% nest() %>%
  mutate(mod = map(data, ~lm(genotype ~ value, data = .)),
         tid = map(mod, tidy)) %>% select(-data,-mod) %>% unnest(tid)

d2 <- d %>%
  filter(term != "(Intercept)" ) %>%
  left_join(geno_map) %>%
  left_join(target_prots2)

data.table::fwrite(d2, "data_share/associations_pqlts.txt", sep = "\t")

chrono <- sleep <- data.table::fread(project_path("chronotype2.tsv")) %>%
  select(eid,
         h_sleep = `1160-0.0`,
         chrono = `1180-0.0`,
         ever_insomnia = `1200-0.0`,
         wakeup = `1170-0.0`,
         snoring = `1210-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ "Definitely morning",
    chrono == 2 ~ "Rather morning",
    chrono == -1~ "Don't know",
    chrono == 3 ~ "Rather evening",
    chrono == 4 ~ "Definitely evening",
    TRUE ~ NA_character_),
    chrono = factor(chrono, levels = c("Definitely morning", "Rather morning", "Don't know", "Rather evening", "Definitely evening"))) %>% select(eid, chrono)


d <- res %>% select(-time_day, -pred_mean) %>%
  left_join(combined %>% select(eid = FID, starts_with("rs"))) %>%
  left_join(chrono)

summary(lm(res ~ ., data = d))

