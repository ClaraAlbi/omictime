install.packages("renv")
renv::restore()  # Reads renv.lock and installs everything


s <- data.table::fread("blood_sampling.tsv") %>% janitor::clean_names() %>%
  mutate(max_time = pmax(x3166_0_0, x3166_0_1, x3166_0_2, x3166_0_3, x3166_0_4, x3166_0_5, na.rm = T)) %>%
  separate(max_time, into = c("date", "time"), sep = " ") %>%
  separate(time, into = c("h", "min", "s"), sep = ":") %>%
  separate(date, into = c("y", "m", "d"), sep = "-", remove = F) %>%
  mutate(across(y:s, as.numeric),
         time_day = h + min/60) %>%
  rename(fasting = x74_0_0,
         num_bsamples = x68_0_0,
         date_bsampling = date)
