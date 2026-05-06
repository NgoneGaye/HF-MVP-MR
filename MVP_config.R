# =============================================================================
# MVP_config.R
# Project: Forward MR — eGFR → HFpEF / HFrEF (MVP AFR)
# BMC Nephrology paper
# NOTE: No reverse MR — insufficient HF subtype instruments in AFR ancestry
# =============================================================================

# ── OpenGWAS token (update monthly) ──
# Get new token at: https://api.opengwas.io
Sys.setenv(OPENGWAS_JWT = "PLACE YOUR TOKEN HERE")


library(here)

project_root <- here::here()

# --- Exposure: eGFR ---
egfr_cfg <- list(
  file     = file.path(project_root, "data_raw",
             "eGFRcrea_GWAS_AFR_67K_Individuals_Summary_Statistics.txt.gz"),
  id       = "eGFRcrea_AFR_67K",
  label    = "eGFR (AFR, N=67,943)",
  ancestry = "AFR",
  col_map  = c(
    MarkerName = "SNP",
    CHR        = "chr",
    POS        = "pos",
    ALT        = "effect_allele",
    REF        = "other_allele",
    Freq_ALT   = "eaf",
    BETA       = "beta",
    SE         = "se",
    P.value    = "pval",
    N          = "samplesize"
  ),
  has_eaf  = TRUE,
  type     = "continuous"
)

# --- Outcomes: MVP HF subtypes ---
# IMPORTANT: col_map entries marked "TO CONFIRM" must be verified with
# fread(file, nrows = 5) before running. Update and commit before analysis.

outcomes_cfg <- list(

  # ── HFpEF (Diastolic HF) ──
  # GCST90477963
  # N_cases = 5,379 | N_controls = 113,041 | N_eff = 20,539
  MVP_HFpEF_AFR = list(
    file     = file.path(project_root, "data_raw",
               "GCST90477963.tsv.gz"),  # update filename
    id       = "HFpEF_MVP_AFR",
    label    = "HFpEF (MVP, AFR)",
    ancestry = "AFR",
    col_map = c(
      chromosome              = "chr",
      base_pair_location      = "pos",
      rsid                    = "SNP",
      effect_allele_frequency = "eaf",
      p_value                 = "pval",
      num_cases               = "ncase",
      num_controls            = "ncontrol",
      ci_upper                = "ci_upper",
      ci_lower                = "ci_lower",
      effect_allele           = "effect_allele",   # already named correctly
      other_allele            = "other_allele",    # already named correctly
      odds_ratio              = "odds_ratio"       # needed for derived beta
    ),
    derived = list(
      beta       = "log(as.numeric(odds_ratio))",
      se         = "(log(as.numeric(ci_upper)) - log(as.numeric(ci_lower))) / (2 * 1.96)",
      samplesize = "4 / (1/as.numeric(ncase) + 1/as.numeric(ncontrol))"
    ),
    has_eaf  = TRUE,
    type     = "binary"
  ),

  # ── HFrEF (Systolic HF) ──
  # GCST90477961
  # N_cases = 9,104 | N_controls = 109,632 | N_eff = 33,624
  MVP_HFrEF_AFR = list(
    file     = file.path(project_root, "data_raw",
               "GCST90477961.tsv.gz"),  # update filename
    id       = "HFrEF_MVP_AFR",
    label    = "HFrEF (MVP, AFR)",
    ancestry = "AFR",
    col_map = c(
      chromosome              = "chr",
      base_pair_location      = "pos",
      rsid                    = "SNP",
      effect_allele_frequency = "eaf",
      p_value                 = "pval",
      num_cases               = "ncase",
      num_controls            = "ncontrol",
      ci_upper                = "ci_upper",
      ci_lower                = "ci_lower",
      effect_allele           = "effect_allele",   # already named correctly
      other_allele            = "other_allele",    # already named correctly
      odds_ratio              = "odds_ratio"       # needed for derived beta
    ),
    derived = list(
      beta       = "log(as.numeric(odds_ratio))",
      se         = "(log(as.numeric(ci_upper)) - log(as.numeric(ci_lower))) / (2 * 1.96)",
      samplesize = "4 / (1/as.numeric(ncase) + 1/as.numeric(ncontrol))"
    ),
    has_eaf  = TRUE,
    type     = "binary"
  )
)

# --- Instrument selection ---
instrument_cfg <- list(
  p_strict  = 5e-8,
  p_lenient = 1e-6,   # kept for sensitivity / power reporting
  clump_kb  = 10000,
  clump_r2  = 0.001,
  F_min     = 10
)

# --- MR methods ---
mr_methods_binary <- c("mr_ivw_mre", "mr_weighted_median",
                       "mr_egger_regression", "mr_weighted_mode")

# --- Output directory ---
out_root <- file.path(project_root, "results")
