# Population-scale molecular reconstruction of human circadian phase from blood biomarkers

This repository contains the analysis code and fitted proteomic prediction models accompanying the manuscript:

> Albiñana C et al. **Population-scale molecular reconstruction of human circadian phase from blood biomarkers.** medRxiv. 2026. [https://doi.org/10.64898/2026.07.08.26356418](https://www.medrxiv.org/content/10.64898/2026.07.08.26356418v1)

## Overview

Circadian timing influences human physiology and disease risk, but molecular circadian phase is difficult to measure at population scale. This project uses circulating blood biomarkers from UK Biobank to:

- characterize time-of-day variation across plasma biomarkers;
- reconstruct sampling time using machine-learning models trained on plasma biomarkers;
- define circadian acceleration as deviation from the population-average predicted time;
- examine associations with chronotype, shift work, phenotypes, and disease; and
- investigate the heritability and genetic correlates of circadian acceleration.

Please refer to the [medRxiv manuscript](https://www.medrxiv.org/content/10.64898/2026.07.08.26356418v1) for the study design, methods, results, limitations, author list, and supplementary information.

## Repository contents

```text
code/
├── 00_variable_extraction/
├── 01_variance_analysis/
├── 02_prediction_and_validation/
├── 03_phenotype_associations/
└── 04_genetics/
models/
└── proteomic_prediction/
```

- [`code/`](code/) contains the R and shell scripts used for variable preparation, variance and harmonic analyses, prediction and validation, phenotype associations, and genetic analyses.
- [`code/README.md`](code/README.md) describes the code structure, dependencies, and recommended execution order.
- [`models/proteomic_prediction/`](models/proteomic_prediction/) contains fitted models for the 14-panel Olink proteomic time-of-day prediction analysis.
- [`models/proteomic_prediction/README.md`](models/proteomic_prediction/README.md) documents the LASSO, squared-feature LASSO, LightGBM, and XGBoost model files.

## Requirements

The primary analyses are written in R. Selected genetic workflows additionally require platform-specific tools and environments, including the UK Biobank Research Analysis Platform, PLINK, GCTA, LDSC, SBayesRC, PolyFun/FINEMAP, and ANNOVAR.

Package and workflow details are provided in the [code documentation](code/README.md).

## Citation

If you use this code or the accompanying models, please cite:

```text
Albiñana C, et al. Population-scale molecular
reconstruction of human circadian phase from blood biomarkers.
medRxiv. 2026. doi:10.64898/2026.07.08.26356418
```

## Contact

Questions about the methods and data should be directed to the corresponding authors listed in the manuscript.
