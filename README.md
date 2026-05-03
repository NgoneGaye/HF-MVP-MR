# eGFR → Heart Failure Subtypes: Mendelian Randomisation in African Ancestry Populations

**Ngone Diaba Gaye** | Cardiology & Genomics Research | Dakar, Senegal

---

## Overview

This repository contains the full analysis pipeline for a two-sample Mendelian randomisation (MR) study examining the causal effect of kidney function (estimated glomerular filtration rate, eGFR) on heart failure and its subtypes — heart failure with preserved ejection fraction (HFpEF) and heart failure with reduced ejection fraction (HFrEF) — in populations of African ancestry.

The study uses genome-wide association study (GWAS) summary statistics from African ancestry cohorts and applies multiple MR methods to assess causality while accounting for potential pleiotropy.

---

## Repository Structure

```
HF-MVP-MR/
├── MVP_config.R                  # Global settings, library loading, file paths
├── MVP_functions.R               # Custom MR functions and helpers
├── MVP_forward.R                 # Main forward MR analysis (eGFR → HF subtypes)
├── MVP_combined_forest.R         # Combined forest plot across MR methods
├── MVP_figures.R                 # Main manuscript figures
├── MVP_supplementary_figures.R   # Supplementary figures
├── MVP_tables.R                  # Manuscript and supplementary tables
├── data_raw/                     # Raw GWAS summary data (not tracked by git)
├── results/                      # Output files (not tracked by git)
└── session_info.txt              # R session and package version information
```

> **Note:** Raw data files are not included in this repository. GWAS summary statistics used in this study are publicly available. Sources and download links will be added upon preprint submission.

---

## Methods

- **Study design:** Two-sample Mendelian randomisation
- **Exposure:** Kidney function (eGFR), instrumented by independent genome-wide significant SNPs from African ancestry GWAS
- **Outcomes:** Heart failure (overall), HFpEF, HFrEF
- **Ancestry:** African
- **MR methods:** Inverse-variance weighted (IVW), MR-Egger, weighted median, weighted mode, simple mode
- **Sensitivity analyses:** Heterogeneity testing (Cochran's Q), MR-Egger intercept, leave-one-out analysis
- **Software:** R (TwoSampleMR, MendelianRandomization, ggplot2)

---

## Requirements

R version ≥ 4.2.0

Key packages:
```r
TwoSampleMR
MendelianRandomization
ggplot2
dplyr
tidyr
readr
```

Install core packages:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr"))
remotes::install_github("MRCIEU/TwoSampleMR")
```

Full session information including all package versions is available in `session_info.txt`.

---

## How to Reproduce

1. Clone the repository:
```bash
git clone https://github.com/NgoneGaye/HF-MVP-MR.git
cd HF-MVP-MR
```

2. Open `HF MVP MR.Rproj` in RStudio

3. Run scripts in this order:
```r
source("MVP_config.R")         # Load settings and paths
source("MVP_functions.R")      # Load custom functions
source("MVP_forward.R")        # Run main MR analysis
source("MVP_tables.R")         # Generate tables
source("MVP_figures.R")        # Generate main figures
source("MVP_combined_forest.R")        # Generate forest plots
source("MVP_supplementary_figures.R")  # Generate supplementary figures
```

---

## Status

Analysis complete. Manuscript in preparation. Preprint submission planned to medRxiv.

---

## Contact

**Dr Ngone Diaba Gaye**  
Cardiologist & Researcher  
Dakar, Senegal  
GitHub: [@NgoneGaye](https://github.com/NgoneGaye)

---

## Licence

This code is shared for transparency and reproducibility. If you use or adapt this pipeline, please cite the associated preprint (link to be added upon submission).
