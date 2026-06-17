# eGFR → Heart Failure Subtypes: Mendelian Randomisation in African Ancestry Populations

**Ngone Diaba Gaye** | Cardiology & Genomics Research | Dakar, Senegal

---

## Overview

This repository contains the full analysis pipeline for a two-sample Mendelian randomisation (MR) study examining the causal effect of kidney function (estimated glomerular filtration rate, eGFR) on heart failure subtypes, heart failure with preserved ejection fraction (HFpEF) and heart failure with reduced ejection fraction (HFrEF), in populations of African ancestry.

The study uses genome-wide association study (GWAS) summary statistics from African ancestry cohorts and applies multiple MR methods to assess causality while accounting for potential pleiotropy.

---

## Repository Structure
HF-MVP-MR/
├── HF MVP MR.Rproj               # RStudio project file
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

> **Note:** Raw data files are not included in this repository. Exposure GWAS summary statistics (eGFR, African ancestry, N = 67,943) are publicly available via figshare. Outcome GWAS summary statistics are available from the GWAS Catalog under accession numbers GCST90477963 (HFpEF) and GCST90477961 (HFrEF). Direct links are provided in the Data Sources section below.

---

## Data Sources

- **Exposure (eGFR):** [figshare](https://figshare.com/ndownloader/files/47668984)
- **Outcome (HFpEF):** GWAS Catalog, accession [GCST90477963](https://www.ebi.ac.uk/gwas/studies/GCST90477963)
- **Outcome (HFrEF):** GWAS Catalog, accession [GCST90477961](https://www.ebi.ac.uk/gwas/studies/GCST90477961)

---

## Methods

- **Study design:** Two-sample Mendelian randomisation
- **Exposure:** Kidney function (eGFR), instrumented by independent genome-wide significant SNPs from an African ancestry GWAS
- **Outcomes:** HFpEF, HFrEF
- **Ancestry:** African
- **Primary MR method:** Inverse-variance weighted (IVW), random-effects
- **Sensitivity analyses:** Weighted median, MR-Egger, weighted mode
- **Additional checks:** Heterogeneity testing (Cochran's Q), MR-Egger intercept test, leave-one-out analysis, Steiger filtering
- **Software:** R 4.5.0 (TwoSampleMR, ieugwasr, data.table, ggplot2, flextable, officer)

---

## Requirements

R version 4.5.0

Key packages:
```r
TwoSampleMR
ieugwasr
data.table
ggplot2
flextable
officer
```

Install core packages:
```r
install.packages(c("ggplot2", "data.table", "flextable", "officer"))
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/ieugwasr")
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
source("MVP_config.R")                 # Load settings and paths
source("MVP_functions.R")              # Load custom functions
source("MVP_forward.R")                # Run main MR analysis
source("MVP_tables.R")                 # Generate tables
source("MVP_figures.R")                # Generate main figures
source("MVP_combined_forest.R")        # Generate forest plots
source("MVP_supplementary_figures.R")  # Generate supplementary figures
```

---

## Status

Analysis complete. Manuscript submitted to *BMC Nephrology*.

---

## Contact

**Dr Ngone Diaba Gaye**
Cardiologist & Researcher
Dakar, Senegal
GitHub: [@NgoneGaye](https://github.com/NgoneGaye)

---

## Licence

This code is shared for transparency and reproducibility. If you use or adapt this pipeline, please cite the associated manuscript once published in *BMC Nephrology* (citation details to be added upon acceptance).