
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

olink_cohort <- data.table::fread("/mnt/project/olink_instance_0.csv")

sleep <- data.table::fread("/mnt/project/chronotype2.tsv") %>%
  select(eid,
         chrono = `1180-0.0`) %>%
  mutate(chrono = case_when(
    chrono == 1 ~ 2,
    chrono == 2 ~ 1,
    chrono == -1~ 0,
    chrono == 3 ~ -1,
    chrono == 4 ~ -2))

df <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  rowwise() %>%
  mutate(pred_mean = mean(c(pred_lgb, pred_xgboost, pred_lasso, pred_lassox2))) %>%
  filter(!is.na(time_day)) %>%
  filter(i == 0) %>%
  select(eid, time_day, pred_mean)
df$res <- residuals(lm(pred_mean ~ time_day, data = df))

d <- olink_cohort %>%
  left_join(sleep) %>%
  left_join(df)

summary(lm(chrono ~ res, data = d))
summary(lm(chrono ~ spon2, data = d))
summary(lm(chrono ~ relt, data = d))
summary(lm(chrono ~ gdf15, data = d))


  ggplot(aes(x = chrono, y = spon2)) + geom_point()





olink_t <- readRDS("/mnt/project/biomarkers_res/res_olink_tech_14panels.rds")

cor_pmat <- function(mat, method = "pearson") {
    vars <- colnames(mat)

    expand.grid(prot1 = vars, prot2 = vars, stringsAsFactors = FALSE) %>%
      rowwise() %>%
      mutate(
        test = list(cor.test(
          mat[, prot1],
          mat[, prot2],
          method = method
        )),
        r = test$estimate,
        p = test$p.value
      ) %>%
      ungroup() %>%
      select(-test)
  }


target_prots <- c(
    "MYOC","HYAL1","EFNA1","TNR","FAS","RELT","CD276",
    "SPON2","SPINK5","GDF15","PGF","ANGPTL1",
    "KLK12","LGALS1","HS3ST3B1", "KLK13"
  )

mat <- olink_t %>%
  rename_with(tolower) %>%
  select(any_of(tolower(target_prots))) %>%
  as.matrix()

cs_long <- cor_pmat(mat) %>%
  mutate(pfdr = p.adjust(p))


ggplot(cs_long, aes(x = prot1, y = prot2, fill = r)) +
  geom_tile(color = "white") +
  #geom_text(aes(label = ifelse(pfdr < 0.05, "*", "")), size = 5) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "r"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  )


time <- readRDS("/mnt/project/biomarkers/time.rds")

olink_t %>%
  rename_with(tolower) %>%
  select(eid, any_of(tolower(target_prots))) %>%
  left_join(time) %>%
  filter(time_day > 9 & time_day < 20) %>%
  pivot_longer(2:17) %>%
  ggplot(aes(x = time_day, y = value, color = name)) +
  geom_smooth() +
  theme_classic(base_size = 16)
