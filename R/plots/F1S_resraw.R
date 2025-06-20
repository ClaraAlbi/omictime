

library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(lubridate)
library(ggplot2)
install.packages("ggtext")

time <- readRDS("/mnt/project/biomarkers/time.rds")

fields <- data.table::fread("field.tsv")

df_r2 <- bind_rows(readRDS("/mnt/project/biomarkers_3/covariate_res/aov_labs.rds") %>% mutate(type = "Biochemistry") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("/mnt/project/biomarkers_3/covariate_res/aov_counts.rds") %>% mutate(type = "Cell_counts") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id")),
                   readRDS("/mnt/project/biomarkers_3/covariate_res/aov_nmr.rds") %>% mutate(type = "Metabolomics-NMR") %>%
                     left_join(fields %>% select(field_id, title), by = c("phen" = "field_id"))) %>%
  mutate(phen = as.character(phen)) %>%
  bind_rows(readRDS("/mnt/project/biomarkers_3/covariate_res/aov_olink.rds") %>% mutate(type = "Proteomics-Olink") %>%
              mutate(title = phen)
  ) %>%
  mutate(color_var = case_when(type == "Proteomics-Olink" ~ "#76B041",
                               type == "Metabolomics-NMR" ~ "#2374AB",
                               type == "Cell_counts" ~ "#8F3985",
                               type == "Biochemistry" ~ "#E85F5C")) %>%
  mutate(title = case_when(title == "White blood cell (leukocyte) count" ~ "Leukocyte count",
                           title == "Phospholipids to Total Lipids in Small HDL percentage" ~ "Phosphlipid ratio SHDL",
                           title == "Cholesterol to Total Lipids in Very Large HDL percentage" ~ "Cholesterol ratio VLHDL",
                           title == "Phospholipids to Total Lipids in Very Large HDL percentage" ~ "Phospholipids ratio VLHDL",
                           title == "Cholesterol to Total Lipids in Small HDL percentage" ~ "Cholesterol ratio SHLD",
                           title == "Spectrometer-corrected alanine" ~ "Alanine",
                           TRUE ~ title)) %>%
  filter(term == "time_day")

df_top <- df_r2 %>%
  group_by(type, phen) %>%
  filter(any(p.value < 0.05)) %>%
  ungroup() %>%
  distinct(phen, .keep_all = TRUE) %>%
  arrange(desc(pr2)) %>%
  slice_head(n = 30) %>%
  select(phen) %>%
  inner_join(df_r2, by = "phen") %>%
  filter(term != "Residuals") %>%
  group_by(phen, color_var, type, title) %>%
  summarise(t_r2 = sum(pr2))

facet_levels <- df_top %>%
  arrange(desc(t_r2)) %>%
  # recreate the exact HTML string you’ll use below
  mutate(f_html = sprintf("<span style='color:%s'>%s</span>", color_var, title)) %>%
  pull(f_html)

prot <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_olink.rds") %>%
  select(eid, any_of(df_top$title))

cells <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_counts.rds") %>%
  select(eid, any_of(df_top$phen)) %>%
  filter(eid %in% prot$eid)

nmr <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_nmr.rds") %>%
  select(eid, any_of(df_top$phen)) %>%
  filter(eid %in% prot$eid)

bio <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_labs.rds") %>%
  select(eid, any_of(df_top$phen)) %>%
  filter(eid %in% prot$eid)

raw <- prot %>%
  left_join(cells) %>% left_join(nmr) %>% left_join(bio) %>%
  left_join(time %>% select(eid, time_day)) %>%
  pivot_longer(c(-eid, -time_day), names_to = "phen")

prot_res <- readRDS("/mnt/project/biomarkers_3/covariate_res/res_olink.rds") %>%
  select(eid, any_of(df_top$title))

cells_res <- readRDS("/mnt/project/biomarkers_3/covariate_res/res_counts.rds") %>%
  select(eid, any_of(df_top$phen)) %>%
  filter(eid %in% prot$eid)

nmr_res <- readRDS("/mnt/project/biomarkers_3/covariate_res/res_nmr.rds") %>%
  select(eid, any_of(df_top$phen)) %>%
  filter(eid %in% prot$eid)

bio_res <- readRDS("/mnt/project/biomarkers_3/covariate_res/res_labs.rds") %>%
  select(eid, any_of(df_top$phen)) %>%
  filter(eid %in% prot$eid)

res <- prot_res %>%
  left_join(cells_res) %>% left_join(nmr_res) %>% left_join(bio_res) %>%
  left_join(time %>% select(eid, time_day)) %>%
  pivot_longer(c(-eid, -time_day), names_to = "phen")


pl_res<- res %>%
  group_by(t = round(time_day, 0), phen) %>%
  summarise(
    n        = n(),               # sample size
    mean_val = mean(value, na.rm = T),       # replace y_var with your outcome
    sd_val   = sd(value, na.rm = T),         # standard deviation
    se_val   = sd_val / sqrt(n),  # standard error
    ci_lower = mean_val - 1.96 * se_val,  # lower 95% CI
    ci_upper = mean_val + 1.96 * se_val   # upper 95% CI
  ) %>%
  left_join(df_r2 %>%
              select(phen, title, type, color_var)) %>%
  mutate(f_html = sprintf("<span style='color:%s'>%s</span>", color_var, title),
         facet_html = factor(f_html, levels = facet_levels)) %>%
  ggplot(aes(x = t, y = mean_val, color = type)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(x = "Recorded time") +
  scale_x_continuous(breaks = c(8, 10, 12, 14, 16, 18, 20), limits = c(8, 21)) +
  scale_color_manual(
    name   = "Omic type",
    values = c(
      "Proteomics-Olink"  = "#76B041",
      "Metabolomics-NMR"  = "#2374AB",
      "Cell_counts"       = "#8F3985",
      "Biochemistry"      = "#E85F5C"
    )) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 15,
        size  = 6
      ), nrow = 2, byrow = TRUE
    )
  ) +
  theme_minimal() +
  ggtitle("Residualised values") +
  facet_wrap(~facet_html, scales = "free", ncol = 4) +
  theme(text = element_text(size = 14),
        strip.text = ggtext::element_markdown(size = 14, hjust = 0),
        axis.title.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        title = element_text(size = 18), legend.text = element_text(size = 16))

ggsave(plot = pl_raw, filename =  "plots/F1S_raw.png", height = 18, width = 15)
ggsave(plot = pl_res, filename =  "plots/F1S_res.png", height = 12, width = 15)
