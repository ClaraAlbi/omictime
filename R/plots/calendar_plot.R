library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)


covs <- data.table::fread("/mnt/project/covariates.tsv") %>%
  mutate(across(c(2,3, 7, 9), as.factor),
         bmi = `21002-0.0` / (`50-0.0`/100)^2,
         smoking = case_when(`20116-0.0` == "-3" ~ NA_character_,
                             TRUE ~ `20116-0.0`)) %>% select(-"50-0.0", -"21002-0.0", -"20116-0.0")

colnames(covs) <- c("eid", "sex", "birth_year",  "age_recruitment",  "assessment_centre", "month_attending", "bmi", "smoking")

codings <- data.table::fread("coding10.tsv")

time  <- data.table::fread("/mnt/project/blood_sampling.tsv") %>%
  mutate(max_time = pmax(`3166-0.0`, `3166-0.1`, `3166-0.2`, `3166-0.3`, `3166-0.4`, `3166-0.5`, na.rm = T))
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  filter(date != "1900-01-01")

time <- readRDS("/mnt/project/biomarkers/time.rds") %>%
  filter(!y %in% c(1900, 2006))

time_i1 <- data.table::fread("/mnt/project/blood_sampling_instance1.tsv") %>%
  mutate(max_time = pmax(`3166-1.0`,`3166-1.1`,`3166-1.2`,`3166-1.3`,`3166-1.4`, `3166-1.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  filter(date != "1900-01-01")

time_i2 <- data.table::fread("/mnt/project/blood_sampling_instance2.tsv") %>%
  filter(!is.na(`3166-2.0`)) %>%
  mutate(max_time = pmax(`3166-2.0`,`3166-2.1`,`3166-2.2`,`3166-2.3`,`3166-2.4`, `3166-2.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  filter(date != "1900-01-01")

time_i3 <- data.table::fread("/mnt/project/blood_sampling_instance3.tsv") %>%
  filter(!is.na(`3166-3.0`)) %>%
  mutate(max_time = pmax(`3166-3.0`,`3166-3.1`,`3166-3.2`,`3166-3.3`,`3166-3.4`, `3166-3.5`, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  filter(date != "1900-01-01")

#########


cs <- sort(table(time$date_bsampling), decreasing = T)

cal <- time %>%
  # pull out year, week-of-year, weekday
  mutate(
    appointment_date = date_bsampling,
    year    = year(appointment_date),
    # lubridate::week gives 1–53 (ISO week might jump to prev/following year around Jan)
    week    = week(appointment_date),
    weekday = wday(appointment_date, label=TRUE, abbr=TRUE, week_start=1)
  ) %>%
  # count per calendar cell
  count(year, week, weekday, name="n") %>%
  ungroup() %>%
  # fill in any missing days (so empty days show white)
  complete(
    year,
    week    = full_seq(1:53, 1),
    weekday = wday(1:7, label=TRUE, abbr=TRUE, week_start=1),
    fill = list(n = 0)
  ) %>% filter(!is.na(year))

# ── 3) Plot ──────────────────────────────────────────────────────────────────
ggplot(cal_i1, aes(x = weekday, y = week, fill = n)) +
  geom_tile(color="white") +
  # week 1 at top
  scale_y_reverse(expand = c(0,0), breaks = seq(1,53,2)) +
  # nice perceptual palette
  scale_fill_viridis_c(option="magma", name="# appts") +
  # one panel per year
  facet_wrap(~year, nrow = 1) +
  # clean up labels
  labs(
    x = NULL,
    y = "Week of year",
    title = "Daily appointment counts (calendar view), instance 1"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 6),
    strip.background = element_rect(fill="#EEEEEE", color=NA),
    strip.text = element_text(face="bold")
  )

### histograms


time_i2 %>%
  inner_join(covs) %>%
  left_join(codings, by = c("assessment_centre" = "coding")) %>%
  ggplot(aes(x = time_day)) +
  geom_histogram(bins = 60) +
  coord_polar(start = 0) +
  labs(x = "Time of day") +
  scale_x_continuous(limits = c(0, 24), breaks = c(0, 3, 6, 9, 12, 15, 18, 21)) +
  facet_wrap(~meaning) +
  theme_minimal() +
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_blank(), panel.grid.minor = element_blank())
