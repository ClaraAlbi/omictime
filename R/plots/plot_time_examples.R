
library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(lubridate)
library(ggplot2)

time <- readRDS("/mnt/project/biomarkers/time.rds")

cells <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_counts.rds") %>%
  select(eid, `30120`)

pc <- cells %>%
  select(eid, `30120`) %>%
  left_join(time %>% select(eid, time_day)) %>%
  ggplot(aes(x = time_day, y = `30120`)) +
  geom_smooth(color = "#8F3985") +
  labs(x = "Recorded time") +
  scale_x_continuous(breaks = c(8, 10, 12, 14, 16, 18, 20)) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank())


ggsave(plot = pc, filename =  "lymphocites.png", height = 3, width = 4)

bio <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_labs.rds") %>%
  select(eid, `30740`)

pb <- bio %>%
  select(eid, `30740`) %>%
  left_join(time %>% select(eid, time_day)) %>%
  ggplot(aes(x = time_day, y = `30740`)) +
  geom_smooth(color = "#E85F5C") +
  labs(x = "Recorded time") +
  scale_x_continuous(breaks = c(8, 10, 12, 14, 16, 18, 20)) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank())

ggsave(plot = pb, filename =  "glucose.png", height = 3, width = 4)





prot <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_olink.rds") %>%
  select(eid, pomc)

pp <- prot %>%
  select(eid, pomc) %>%
  left_join(time %>% select(eid, time_day)) %>%
  ggplot(aes(x = time_day, y = pomc)) +
  geom_smooth(color = "#76B041") +
  labs(x = "Recorded time") +
  scale_x_continuous(breaks = c(8, 10, 12, 14, 16, 18, 20)) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank())

ggsave(plot = pp, filename =  "pomc.png", height = 3, width = 4)



nmr <- readRDS("/mnt/project/biomarkers_3/covariate_res/raw_nmr.rds") %>%
  select(eid, `23471`)

pn <- nmr %>%
  select(eid, `23471`) %>%
  left_join(time %>% select(eid, time_day)) %>%
  ggplot(aes(x = time_day, y = `23471`)) +
  geom_smooth(color = "#2374AB") +
  labs(x = "Recorded time") +
  scale_x_continuous(breaks = c(8, 10, 12, 14, 16, 18, 20)) +
  theme_minimal() +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title.y = element_blank())

ggsave(plot = pn, filename =  "lactate.png", height = 3, width = 4)

