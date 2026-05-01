# =============================================================================
# MVP_forward.R
# Forward MR: eGFR → HFpEF and HFrEF (MVP AFR)
# BMC Nephrology paper
# =============================================================================

library(here)
source(here::here("MVP_config.R"))
source(here::here("MVP_functions.R"))

if (length(outcomes_cfg) == 0) stop("outcomes_cfg is empty — check MVP_config.R")

# Confirm col_maps have been verified before proceeding
cat("\n====== MVP FORWARD MR: eGFR → HFpEF / HFrEF ======\n\n")
cat("REMINDER: Confirm MVP col_map entries before first run.\n")
cat("Run: fread('your_mvp_file.gz', nrows=5) and check names(dt)\n\n")

# ── STEP 1: Load eGFR exposure ────────────────────────────────────────────────
cat("--- Loading eGFR exposure ---\n")
egfr_dt <- load_gwas(egfr_cfg)
print(qc_summary(egfr_dt, egfr_cfg$id))

# ── STEP 2: Clump eGFR instruments (once, reused for both outcomes) ───────────
cat("\n--- Clumping eGFR instruments ---\n")

clump_dir <- file.path(out_root, "forward", "eGFR_instruments")
make_out_dirs(clump_dir)

egfr_instruments <- format_and_clump(
  dt             = egfr_dt,
  cfg            = egfr_cfg,
  p_thresh       = instrument_cfg$p_strict,
  instrument_cfg = instrument_cfg,
  save_dir       = clump_dir
)

n_preclump  <- nrow(egfr_instruments$preclump)
n_clumped   <- nrow(egfr_instruments$clumped)
exp_clumped <- egfr_instruments$clumped

cat("\n  eGFR instruments retained:", n_clumped, "\n")
if (n_clumped == 0) stop("No instruments after clumping — cannot proceed.")

# ── STEP 3: Loop over HFpEF and HFrEF ────────────────────────────────────────
all_results <- list()

for (out_name in names(outcomes_cfg)) {

  out_cfg <- outcomes_cfg[[out_name]]
  cat("\n\n====== OUTCOME:", out_cfg$label, "======\n")

  out_dir <- file.path(out_root, "forward", out_name)
  make_out_dirs(out_dir)

  # Load outcome
  cat("--- Loading outcome GWAS ---\n")
  out_dt <- load_gwas(out_cfg)
  print(qc_summary(out_dt, out_cfg$id))

  # Format outcome
  cat("\n--- Formatting outcome ---\n")
  out_fmt <- format_outcome(
    dt       = out_dt,
    cfg      = out_cfg,
    snp_list = exp_clumped$SNP,
    save_dir = out_dir
  )

  if (is.null(out_fmt) || nrow(out_fmt) == 0) {
    cat("  No matching SNPs — skipping", out_name, "\n"); next
  }

  n_outcome_match <- nrow(out_fmt)

  # Harmonise + F-stats
  cat("\n--- Harmonising ---\n")
  harm_obj <- harmonise_and_fstat(
    exp_clumped    = exp_clumped,
    out_fmt        = out_fmt,
    instrument_cfg = instrument_cfg,
    exp_id         = egfr_cfg$id,
    exp_label      = egfr_cfg$label,
    out_id         = out_cfg$id,
    out_label      = out_cfg$label,
    save_dir       = out_dir
  )

  dat_strong <- harm_obj$strong

  # MDE
  R2_total  <- sum(dat_strong$F_stat /
                  (dat_strong$F_stat + dat_strong$samplesize.exposure - 2),
                   na.rm = TRUE)
  N_eff_out <- mean(dat_strong$samplesize.outcome, na.rm = TRUE)
  mde       <- calc_mde(R2_total, N_eff_out, outcome_type = "binary")

  cat("  R² total:", round(R2_total, 6), "\n")
  cat("  N_eff outcome:", round(N_eff_out, 0), "\n")
  cat("  MDE (OR):", round(mde$mde_OR, 3), "\n")

  # Save MDE to table
  mde_tbl <- data.frame(
    outcome    = out_cfg$label,
    R2_total   = round(R2_total, 6),
    N_eff      = round(N_eff_out, 0),
    MDE_logOR  = round(mde$mde_logOR, 4),
    MDE_OR     = round(mde$mde_OR, 3),
    power_alpha= 0.05,
    power_pct  = 80
  )
  data.table::fwrite(mde_tbl,
    file.path(out_dir, "tables", "00_power_MDE.csv"))

  # Run MR
  cat("\n--- Running MR ---\n")
  mr_obj <- run_mr_full(
    dat_strong,
    method_list  = mr_methods_binary,
    outcome_type = "binary",
    save_dir     = out_dir
  )

  # Plots
  save_mr_plots(mr_obj, dat_strong, save_dir = out_dir)

  # Attrition
  make_attrition(
    n_preclump      = n_preclump,
    n_clumped       = n_clumped,
    n_outcome_match = n_outcome_match,
    n_harmonised    = nrow(harm_obj$all),
    n_strong        = nrow(dat_strong),
    save_dir        = out_dir
  )

  all_results[[out_name]] <- mr_obj
  cat("\n  Done:", out_name, "\n")
}

# ── STEP 4: Combined results table (HFpEF + HFrEF side by side) ──────────────
cat("\n--- Compiling combined results ---\n")

combined <- dplyr::bind_rows(lapply(names(all_results), function(nm) {
  res <- all_results[[nm]]
  if (is.null(res) || is.null(res$results)) return(NULL)
  res$results %>% dplyr::mutate(analysis = nm)
}))

if (!is.null(combined) && nrow(combined) > 0) {
  data.table::fwrite(combined,
    file.path(out_root, "forward", "COMBINED_HFpEF_HFrEF_results.csv"))
  cat("Combined results saved.\n")
}

cat("\n====== MVP FORWARD ANALYSIS COMPLETE ======\n")
cat("Results saved to:", file.path(out_root, "forward"), "\n")

# ── Session info ──────────────────────────────────────────────────────────────
session_info <- sessionInfo()
print(session_info)

# Save as RDS
saveRDS(session_info,
        file.path(out_root, "session_info.rds"))

# Save as readable text file
sink(file.path(out_root, "session_info.txt"))
print(session_info)
sink()

cat("\nSession info saved to:", file.path(out_root, "session_info.txt"), "\n")