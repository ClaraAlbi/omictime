library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)
install.packages("broom")

type <- "olink_tech_14panels"

# Rank-based inverse normal transform
rint <- function(x) {
  r <- rank(x, ties.method = "average")
  qnorm((r - 0.5) / length(r))
}

# --- read once
time     <- readRDS("/mnt/project/biomarkers/time.rds")
covs_raw <- readRDS("/mnt/project/biomarkers/covs.rds") %>%
  mutate(bmi = weight / (height/100)^2)
gen_covs <- data.table::fread("/mnt/project/genetic_covs.tsv") %>%
  select(eid, `22009-0.1`:`22009-0.20`)
colnames(gen_covs) <- c("eid", paste0("PC", 1:20))

olink_times <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_processing_start_date.dat")
batch       <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_batch_number.dat")
panel       <- data.table::fread("/mnt/project/Bulk/Protein biomarkers/Olink/helper_files/olink_assay.dat") %>%
  mutate(
    prot  = tolower(Assay),
    Panel = str_remove(Panel, " ")
  ) %>%
  select(prot, Panel)

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

prots <- readRDS("olink_instance_0_QC_14_panels.rds")

stopifnot("eid" %in% names(prots))

# --- prep lookups
pcs <- paste0("PC", 1:20)

# Base covariates (joined ONCE)
base <- as_tibble(prots["eid"]) %>%
  left_join(gen_covs %>% select(eid, any_of(pcs)), by = "eid") %>%
  left_join(
    covs_raw %>% select(eid, sex, age_recruitment, assessment_centre, month_attending, bmi, smoking),
    by = "eid"
  ) %>%
  left_join(
    time %>% select(eid, fasting, time_day, date_bsampling, y, m),
    by = "eid"
  ) %>%
  left_join(metadata %>% select(eid, Batch, ppp_sel), by = "eid") %>%
  mutate(
    smoking = if_else(smoking == -3, NA_integer_, smoking),
    across(c(assessment_centre, Batch, ppp_sel, sex, month_attending), as.factor)
  )

# Valid-level filters (>= 50 per level) computed once for stable factors
stable_factor_cols <- c("assessment_centre","Batch","ppp_sel","sex","month_attending","smoking")
valid_levels <- lapply(base[stable_factor_cols], function(f) {
  if (!is.factor(f)) return(NULL)
  names(which(table(f) >= 50L))
})

drop_rare_levels <- function(df) {
  for (nm in names(valid_levels)) {
    if (!nm %in% names(df)) next
    keep <- valid_levels[[nm]]
    if (is.null(keep) || !is.factor(df[[nm]])) next
    x <- df[[nm]]
    x[!(x %in% keep) & !is.na(x)] <- NA
    df[[nm]] <- droplevels(x)
  }
  df
}

# Panel mapping: protein name -> panel column name
prot_to_panel <- setNames(panel$Panel, panel$prot)
panel_names   <- unique(panel$Panel)[1:4]
panel_names   <- panel_names[panel_names %in% names(metadata)]  # keep only those present in metadata

# Fast access to metadata panel columns
metadata_dt <- as.data.table(metadata)
setkey(metadata_dt, eid)

# Terms to exclude from regressions
exclude_terms <- c("bmi", "smoking")

# --- main pass (no batching)
current_prots <- setdiff(names(prots), "eid")
nP <- length(current_prots)

l1 <- vector("list", nP)  # aov tidies
l2 <- vector("list", nP)  # residuals (eid, raw, res)
l3 <- vector("list", nP)  # cos/sin effects

# Pre-convert for speed
prots_dt <- as.data.table(prots)

for (k in seq_along(current_prots)) {

  var <- current_prots[k]
  var_lc <- tolower(var)

  idx <- match(tolower(var), panel$prot)
  phen_panel <- if (is.na(idx)) NULL else panel$Panel[idx]

  # Build working frame: base + protein raw + panel_date (if mapped)
  d <- base %>% mutate(raw = prots_dt[[var]])

  if (!is.null(phen_panel) && phen_panel %in% names(metadata_dt)) {
    d$panel_date <- factor(metadata_dt[[phen_panel]][match(d$eid, metadata_dt$eid)])
  }

  # Essential filters
  keep <- !is.na(d$raw) &
    is.finite(d$time_day) & d$time_day > 0 &
    is.finite(d$fasting)  & d$fasting < 24 &
    d$time_day >= 9 & d$time_day <= 20

  if (!any(keep)) next
  d <- d[keep, , drop = FALSE]

  # Drop rare levels for stable factors; and for panel_date specifically (>=50)
  d <- drop_rare_levels(d)
  if ("panel_date" %in% names(d)) {
    lv <- names(which(table(d$panel_date) >= 30L))
    d$panel_date <- droplevels(d$panel_date[d$panel_date %in% lv])
  }

  # INT
  d$rint_b <- rint(d$raw)

  # ----- ANOVA (explicit terms) -- EXCLUDING sex, age_recruitment, bmi, smoking
  candidate_terms <- c("fasting", "sex", "age_recruitment",
                       "assessment_centre","month_attending", pcs,
                       "Batch","ppp_sel","panel_date")
  # intersect with what is actually present in d
  present_terms <- intersect(candidate_terms, names(d))

  f_aov <- reformulate(present_terms, response = "rint_b")

  mod_aov <- aov(f_aov, data = d)
  out_aov <- broom::tidy(mod_aov) %>%
    mutate(pr2 = sumsq / sum(sumsq), phen = var)
  l1[[k]] <- out_aov

  # ----- Residualize EXCLUDING sex, age_recruitment, bmi, smoking
  rhs_res <- intersect(c("fasting", "sex", "age_recruitment",
                         "assessment_centre","month_attending", pcs,
                         "Batch","ppp_sel","panel_date"), names(d))
  # remove excluded terms if any accidentally present (defensive)
  rhs_res <- setdiff(rhs_res, exclude_terms)

  f_res <- reformulate(rhs_res, response = "rint_b")

  m_res  <- lm(f_res, data = d, na.action = na.exclude)
  d$res_rint <- residuals(m_res)

  # ----- Time-only cos/sin on residuals
  m_time <- lm(res_rint ~ cos(2*pi*time_day/24) + sin(2*pi*time_day/24), data = d)
  effects <- broom::tidy(m_time, conf.int = TRUE) %>% mutate(phen = var)

  l2[[k]] <- tibble(eid = d$eid, raw = d$raw, res = d$res_rint)
  l3[[k]] <- effects

  if (k %% 1000 == 0) gc()
}

# --- bind & save
output_aov     <- bind_rows(l1)
output_effects <- bind_rows(l3)
saveRDS(output_aov,     glue("aov_{type}.rds"))
saveRDS(output_effects, glue("effects_{type}.rds"))

# residuals/raw wide by eid
res_long <- bind_rows(l2, .id = "phen_id")
if (nrow(res_long)) {
  phen_map <- tibble(phen_id = as.character(seq_along(current_prots)),
                     phen     = current_prots)
  res_long <- res_long %>% left_join(phen_map, by = "phen_id")

  output_res <- res_long %>%
    select(eid, phen, res) %>%
    pivot_wider(id_cols = eid, names_from = phen, values_from = res)
  saveRDS(output_res, glue("res_{type}.rds"))

  output_raw <- res_long %>%
    select(eid, phen, raw) %>%
    pivot_wider(id_cols = eid, names_from = phen, values_from = raw)
  saveRDS(output_raw, glue("raw_{type}.rds"))
} else {
  warning("No residuals produced; skipping res/raw saves.")
}
