
library(dplyr)
library(stringr)
library(writexl)
library(tidyverse)

### Labels

make_label <- function(title, type, max_len = 15) {

  label <- case_when(

    type == "proteins" ~ toupper(title),

    TRUE ~ title |>
      str_squish() |>
      str_replace_all("^Total\\s+", "") |>    # remove leading "Total "
      str_replace_all("^Concentration of\\s+", "") |>
      str_replace_all("\\bpercentage\\b", "pct") |>
      str_replace_all("Free Cholesterol", "FC") |>
      str_replace_all("Cholesteryl Esters", "CE") |>
      str_replace_all("Cholesterol", "Chol") |>
      str_replace_all("Phospholipids", "PL") |>
      str_replace_all("Triglycerides", "TG") |>
      str_replace_all("Monounsaturated Fatty Acids", "MUFA") |>
      str_replace_all("Very Large", "VL") |>
      str_replace_all("Very Small", "VS") |>
      str_replace_all("Large", "L") |>
      str_replace_all("Medium", "M") |>
      str_replace_all("Small", "S") |>
      str_replace_all("\\s+", "_") |>
      str_replace_all("^_+", "") |>           # remove leading _
      str_replace_all("_+$", "")              # remove trailing _
  )

  # enforce max length
  label <- str_sub(label, 1, max_len)

  label
}


data <- readRDS("data/combined_effects.rds")

d1 <- data %>%
  mutate(label = make_label(title = title, type = type_clean))

d1 <- d1 %>%
  mutate(label = make.unique(label, sep = "_"))

d1 %>%
  mutate(n = nchar(label)) %>%
  filter(n > 15) %>%
  select(phen, label, n) %>%
  arrange(desc(n))

# Duplicated
d1 %>%
  count(label) %>%
  filter(n > 1)

d_out <- d1 %>% select(title, label, type_clean, phen, amplitude_24hfreq, acrophase_24hfreq, pvalue_h, starts_with("estimate"), starts_with("p.value"), starts_with("std"))

colnames(d_out)[1:4] <- c("Biomarker", "Label", "Type", "FID")

writexl::write_xlsx(d_out, "tables/harmonic_models.xlsx")
saveRDS(d_out, "tables/harmonic_models.rds")

### data

data2 <- readRDS("data/combined_variance.rds") %>%
  pivot_wider(id_cols = c(phen, type_clean, title), names_from = term, values_from = c(t_r2, p.value)) %>%
  arrange(type_clean) %>%
  mutate(label = make_label(title = title, type = type_clean),
         label = make.unique(label, sep = "_"))  %>%
  select(title, label, type_clean, everything())

colnames(data2)[1:4] <- c("Biomarker", "Label", "Type", "FID")

writexl::write_xlsx(data2, "tables/variance_covariates.xlsx")
saveRDS(data2, "tables/variance_covariates.rds")



data2 %>%
  filter(t_r2_time_day >= 0.01) %>%
  filter(type_clean == "Proteins") %>%
  pull(title) %>% writeLines(., "data_share/top_rhythmic_prots.txt")
