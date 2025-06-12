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

covs <- data.table::fread("covariates.tsv") %>%
  mutate(across(c(2,3, 6, 7, 9), as.factor),
         bmi = `21002-0.0` / (`50-0.0`/100)^2,
         smoking = case_when(`20116-0.0` == "-3" ~ NA_character_,
                             TRUE ~ `20116-0.0`)) %>% select(-"50-0.0", -"21002-0.0", -"20116-0.0")

colnames(covs) <- c("eid", "sex", "birth_year",  "age_recruitment",  "assessment_centre", "month_attending", "bmi", "smoking")

gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, "22009-0.1":"22009-0.20")
colnames(gen_covs) <- c("eid", paste0("PC", 1:20))
