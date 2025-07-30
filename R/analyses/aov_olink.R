library(tidyr)
library(dplyr)
library(glue)
library(stringr)
install.packages("broom")

type <- "olink"

rint <- function(x) {
  ranks <- rank(x, ties.method = "average")

  # Calculate the rank-inverse normal transformation
  n <- length(ranks)
  qnorm((ranks - 0.5) / n)
}

time <- readRDS("/mnt/project/clara/time.rds")
covs <- readRDS("/mnt/project/clara/covs.rds")
gen_covs <- readRDS("/mnt/project/clara/gencovs.rds")

olink_times <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_processing_start_date.dat")
batch <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_batch_number.dat")
panel <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_assay.dat") %>%
  mutate(prot = tolower(Assay), Panel = str_remove(Panel, " ")) %>% select(prot, Panel)

metadata <- data.table::fread("/mnt/project/clara/olink_i0_meta.csv") %>%
  rename(num_prots = p30900_i0,
         PlateID = p30901_i0,
         well_id = p30902_i0,
         ppp_sel = p30903_i0) %>%
  filter(!is.na(PlateID)) %>%
  left_join(batch) %>%
  left_join(olink_times %>% pivot_wider(names_from = Panel, values_from = Processing_StartDate)) %>%
  mutate(ppp_sel = case_when(ppp_sel == "Yes" ~ 1,
                             TRUE ~ 0)) %>%
  rename_all(~str_remove(.x, " ")) %>%
  mutate(across(Cardiometabolic:OncologyII , factor)) %>%
  filter(Batch != 0)


prots <- data.table::fread("/mnt/project/clara/olink_i0.csv") %>%
  select(-glipr1)


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
      mutate(across(c(assessment_centre, Batch, ppp_sel, sex, age_recruitment, smoking), as.factor),
             rint_b = rint(raw)) %>%
      select(eid,
             raw,
             rint_b,
             time_day,
             sex,
             age_recruitment,
             fasting,
             assessment_centre,
             month_attending,
             any_of(paste0("PC", 1:20)),
             Batch,
             ppp_sel,
             bmi,
             smoking,
             all_of(p)) %>%
      filter(if_all(where(is.factor), ~ . %in% names(which(table(.) >= 50)))) %>%
      mutate(across(where(is.factor), droplevels))  %>%
      filter(time_day >= 9 & time_day <= 20)

    mod <- aov(rint_b ~ . + sex * age_recruitment, data = d %>% select(-eid, -raw))
    out_aov <- broom::tidy(mod) %>% mutate(pr2 = sumsq / sum(sumsq),
                                           phen = var)
    l1 <- c(l1, list(out_aov))

    m <- lm(rint_b ~ . + sex * age_recruitment, data = d %>% select(-eid, -raw, -time_day, -bmi, -smoking), na.action = "na.exclude")
    res_rint <- residuals(m)

    m_time_only_cos <- lm(res_rint ~ cos(2 * pi * time_day / 24) + sin(2 * pi * time_day / 24), data = d)

    out <- tibble(eid = d$eid, raw = d$raw, res = res_rint)
    l2 <- c(l2, list(out))

    effects <- broom::tidy(m_time_only_cos , conf.int = TRUE) %>% mutate(phen = var)
    l3 <- c(l3, list(effects))
  }


  output_aov <- do.call(bind_rows, l1)
  saveRDS(output_aov, glue("aov_{type}_{b}.rds"))

  output_effects <- do.call(bind_rows, l3)
  saveRDS(output_effects, glue("effects_{type}_{b}.rds"))

  output_res <- tibble(res = l2) %>%
    mutate(phen = colnames(prots)[b:end_col]) %>%
    unnest(res) %>%
    pivot_wider(id_cols = eid, names_from = phen, values_from = res)

  saveRDS(output_res, glue("res_{type}_{b}.rds"))

  output_raw <- tibble(res = l2) %>%
    mutate(phen = colnames(prots)[b:end_col]) %>%
    unnest(res) %>%
    pivot_wider(id_cols = eid, names_from = phen, values_from = raw)

  saveRDS(output_raw, glue("raw_{type}_{b}.rds"))

}
