library(stringr)
library(tidyr)
library(dplyr)



time_i0 <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  filter(time_day >= 9 & time_day <= 20)

time_i1 <- data.table::fread("/mnt/project/blood_sampling_instance1.tsv") %>%
  mutate(max_time = pmax(`3166-1.0`,`3166-1.1`,`3166-1.2`,`3166-1.3`,`3166-1.4`, `3166-1.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date) %>%
  filter(!is.na(time_day))

i2_meta <- data.table::fread("/mnt/project/blood_sampling_instance2.tsv") %>%
  filter(eid %in% data.table::fread("/mnt/project/OLINK_i2.tsv")$eid) %>%
  filter(!is.na(`3166-2.0`)) %>%
  mutate(max_time = pmax(`3166-2.0`,`3166-2.1`,`3166-2.2`,`3166-2.3`,`3166-2.4`, `3166-2.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date) %>%
  filter(!is.na(time_day))

i3_meta <- data.table::fread("/mnt/project/blood_sampling_instance3.tsv") %>%
  filter(eid %in% data.table::fread("/mnt/project/OLINK_i3.tsv")$eid) %>%
  filter(!is.na(`3166-3.0`)) %>%
  mutate(max_time = pmax(`3166-3.0`,`3166-3.1`,`3166-3.2`,`3166-3.3`,`3166-3.4`, `3166-3.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(date_bsampling = date) %>%
  filter(!is.na(time_day))


light_band <- data.frame(
  xmin = 5.4,
  xmax = 20.5,
  ymin = -Inf,
  ymax = Inf
)

night_band <- data.frame(
  xmin = c(0, 20.5),
  xmax = c(5.4, 24),
  ymin = -Inf,
  ymax = Inf
)


i0_hist <- time_i0 %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  ggtitle(label = "Initial assessment (2006-2010)", subtitle = paste0("n=", nrow(time_i0))) +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

i1_hist <- time_i1 %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  labs(x = "Time of day") +
  ggtitle(label = "First repeat assessment (2012-13)", subtitle = paste0("n=", nrow(time_i1))) +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())


i2_hist <- i2_meta %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  ggtitle(label ="Imaging (2014+)" , subtitle = paste0("n=", nrow(i2_meta))) +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())

i3_hist <- i3_meta %>%
  filter(time_day < 24 & time_day > 0) %>%
  ggplot(aes(x = time_day)) +
  geom_rect(data = light_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "lightyellow", alpha = 0.3, inherit.aes = FALSE) +
  geom_rect(data = night_band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "gray", alpha = 0.2, inherit.aes = FALSE) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  ggtitle(label = "First repeat imaging (2019+)", subtitle = paste0("n=", nrow(i3_meta))) +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        title = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())


### Histograms
install.packages("cowplot")
library(cowplot)
library(ggplot2)
plot_intval <- plot_grid(i0_hist, i1_hist, i2_hist, i3_hist, nrow = 2)

ggsave("plots/time_histograms.png", plot_intval, width = 12, height = 12)

