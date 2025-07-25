
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)

a <- data.table::fread("/mnt/project/phenotypes.txt")

olink_raw_file <- data.table::fread("/mnt/project/olink_instance_0.csv")


phen <- readRDS("/mnt/project/olink_int_replication.rds") %>%
  filter(i == 0)
phen$res <- residuals(lm(pred_lasso ~ time_day, data = phen))

df <-
phen %>%
  left_join(olink_raw_file %>% select(eid, angptl1, spon2, gdf15)) %>%
  mutate(q_res = ntile(res, 5))

cor(df$res, df$angptl1, use = "complete.obs")
summary(lm(res ~ angptl1, data = df))

ggplot(df, aes(x = time_day, y = gdf15, color = factor(q_res))) +
  ggtitle("ANGPTL1") +
  geom_smooth() +
  theme_classic()


