library(tidyr)
library(dplyr)
library(glue)
library(stringr)
library(broom)

type <- "olink_tech"

rint <- function(x) {
  ranks <- rank(x, ties.method = "average")

  # Calculate the rank-inverse normal transformation
  n <- length(ranks)
  qnorm((ranks - 0.5) / n)
}

time <- readRDS("/mnt/project/biomarkers/time.rds")
covs <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight / (height/100)^2)
gen_covs <- data.table::fread("/mnt/project/covariates.txt") %>%
  select(eid = 1, contains("PC"))

olink_times <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_processing_start_date.dat")
batch <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_batch_number.dat")
panel <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_assay.dat") %>%
  mutate(prot = tolower(Assay), Panel = str_remove(Panel, " ")) %>% select(prot, Panel)

metadata <- data.table::fread("/mnt/project/olink_instance_0_meta.csv") %>%
  rename(
    num_prots = p30900_i0,
    PlateID   = p30901_i0,
    well_id   = p30902_i0,
    ppp_sel   = p30903_i0
  ) %>%
  filter(!is.na(PlateID)) %>%
  left_join(batch, by = "PlateID") %>%
  left_join(
    olink_times %>% pivot_wider(names_from = Panel, values_from = Processing_StartDate),
    by = "PlateID"
  ) %>%
  mutate(
    ppp_sel = if_else(ppp_sel == "Yes", 1L, 0L)
  ) %>%
  rename_with(~ str_remove(.x, " ")) %>%
  mutate(across(Cardiometabolic:OncologyII, factor), .keep = "all")


prots <- readRDS("/mnt/project/olink_instance_0_QC.rds")


for (b in seq(2, ncol(prots), 100)){
  l1 <- list()
  l2 <- list()
  l3 <- list()

  end_col <- min(b + 99, ncol(prots))

  for (i in b:end_col) {

    print(i)
    var <- colnames(prots)[i]
    print(var)
    p <- panel$Panel[which(var == panel$prot)]

    covs_c <- metadata %>% select(eid, Batch, ppp_sel, any_of(p))

    d <- tibble(eid = prots$eid, raw = prots[[i]]) %>%
      drop_na() %>%
      left_join(gen_covs %>% select(eid, any_of(paste0("PC", 1:20)))) %>%
      left_join(covs %>% select(eid, sex, age_recruitment, assessment_centre, month_attending, bmi, smoking)) %>%
      left_join(time %>% select(eid, fasting, time_day, date_bsampling, y, m)) %>%
      left_join(covs_c) %>%
      filter(time_day > 0) %>%
      filter(fasting < 24) %>%
      mutate(smoking = case_when(smoking == -3 ~ NA_integer_,
                                 TRUE ~ smoking),
             across(c(assessment_centre, Batch, ppp_sel, sex, month_attending, age_recruitment, smoking), as.factor),
             rint_b = rint(raw)) %>%
      select(eid,
             raw,
             rint_b,
             fasting,
             assessment_centre,
             any_of(paste0("PC", 1:20)),
             Batch,
             ppp_sel,
             all_of(p)) %>%
      filter(if_all(where(is.factor), ~ . %in% names(which(table(.) >= 50)))) %>%
      mutate(across(where(is.factor), droplevels))  %>%
      filter(time_day >= 9 & time_day <= 20)

    m <- lm(rint_b ~ ., data = d %>% select(-eid, -raw), na.action = "na.exclude")
    res_rint <- residuals(m)

    out <- tibble(eid = d$eid, raw = d$raw, res = res_rint)
    l2 <- c(l2, list(out))

  }

  output_res <- tibble(res = l2) %>%
    mutate(phen = colnames(prots)[b:end_col]) %>%
    unnest(res) %>%
    pivot_wider(id_cols = eid, names_from = phen, values_from = res)

  saveRDS(output_res, glue("res_{type}_{b}.rds"))

}
