# Omic time-of-day analyses

This directory contains the analysis code used to quantify time-of-day variation in blood biomarkers, predict sampling time from multi-omic measurements, test phenotype associations, and investigate genetic architecture.

## Repository layout

| Directory | Purpose | Status |
|---|---|---|
| `00_variable_extraction/` | Prepare sampling-time, demographic, and genetic covariates | 
| `01_variance_analysis/` | Estimate time-of-day and technical variance; fit harmonic models |
| `02_prediction_and_validation/` | Train the main prediction models and assess validation | 
| `02_prediction_and_validation/sensitivity_analyses/` | Alternative cohorts, feature sets, technical adjustments, and subgroup models | 
| `03_phenotype_associations/` | Test phenotype, disease, sleep, medication, chronotype, and ancestry associations | 
| `04_genetics/` | GWAS preparation, GWAS/GREML, COJO, PRS, LDSC, MR, and fine-mapping | 

Files under `Figures/` and `Figures_tables/` generate manuscript figures or formatted tables.

## Configuration

Run scripts from the repository root. 

| Variable | Content |
|---|---|
| `OMICTIME_PROJECT_DIR` | Root of the UK Biobank/RAP project data |
| `OMICTIME_BMRC_DIR` | Root of the BMRC genetics analysis workspace |
| `OMICTIME_DOWNLOADS_DIR` | Local directory for externally downloaded inputs |

## R dependencies

Package installation is intentionally not performed inside analysis scripts. Install the required packages before running an analysis.

Core packages:

- `tidyverse`, `data.table`, `broom`, `glue`, `lubridate`, `purrr`, `Matrix`
- `ggplot2`, `cowplot`, `patchwork`, `ggrepel`, `ggtext`, `paletteer`
- `glmnet`, `lightgbm`, `xgboost`, `survival`
- `readxl`, `writexl`, `table1`

Packages used by selected analyses:

- `ggalluvial`, `ggh4x`, `ggmisc`, `ggpmisc`, `ggpubr`, `tidytext`

## Recommended execution order

1. Run `00_variable_extraction/time_var.R` to prepare sampling-time and covariate objects.
2. Run the platform-specific scripts in `01_variance_analysis/` for blood counts, clinical biochemistry, NMR metabolomics, and Olink proteomics.
3. Run the primary `run_ML_*.R` scripts directly under `02_prediction_and_validation/`.
4. Run only the required external-validation and sensitivity scripts. These are not necessary for the primary models.
5. Run `03_phenotype_associations/phenotypes_associations.R`, followed by the required figure/table scripts.
6. Run genetics workflows only in their target environment: RAP submission scripts on UK Biobank RAP and BMRC scripts on the BMRC cluster.

## Primary and sensitivity analyses

### Primary variance analyses

- `aov_counts.R`: blood-cell counts
- `aov_labs.R`: clinical biochemistry
- `aov_nmr.R`: NMR metabolites
- `aov_olink.R`: Olink proteins

### Primary prediction analyses

- `run_ML_counts.R`: blood-cell counts
- `run_ML_labs.R`: clinical biochemistry
- `run_ML_NMR.R`: NMR metabolites
- `run_ML_olink.R`: Olink proteins
- `run_ML_all.R`: combined biomarker platforms

## Proteomic prediction model weights

Fitted models for the 14-panel Olink proteomic time-of-day prediction analysis are provided separately in `models/proteomic_prediction/`. The directory contains LASSO, squared-feature LASSO, LightGBM, and XGBoost models from cross-validation fold 1. See the README in that directory for file formats and loading instructions.
