# =============================================================================
# MVP_functions.R
# Shared helper functions for MVP forward MR project
# Identical logic to GBMI_functions.R — kept separate for project independence
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(TwoSampleMR)
  library(ieugwasr)
})

# ── 1. Directory setup ────────────────────────────────────────────────────────

make_out_dirs <- function(base_dir) {
  for (sub in c("tables", "plots", "rds")) {
    dir.create(file.path(base_dir, sub), recursive = TRUE, showWarnings = FALSE)
  }
  invisible(base_dir)
}

# ── 2. GWAS loader ────────────────────────────────────────────────────────────

load_gwas <- function(cfg, verbose = TRUE) {
  dt <- data.table::fread(cfg$file, select = names(cfg$col_map))
  if (verbose) message("Loaded: ", cfg$id, " (", nrow(dt), " rows)")

  map <- cfg$col_map
  existing_old <- intersect(names(map), names(dt))
  data.table::setnames(dt, old = existing_old, new = map[existing_old])

  if (!is.null(cfg$derived)) {
    for (new_col in names(cfg$derived)) {
      dt[, (new_col) := eval(parse(text = cfg$derived[[new_col]]))]
    }
  }

  is_acgt <- function(x) toupper(x) %in% c("A", "C", "G", "T")
  dt <- dt[
    nchar(effect_allele) == 1 & nchar(other_allele) == 1 &
      is_acgt(effect_allele) & is_acgt(other_allele)
  ]

  required <- intersect(c("SNP", "beta", "se", "pval"), names(dt))
  dt <- dt[complete.cases(dt[, ..required])]

  if ("chr" %in% names(dt)) dt[, chr := suppressWarnings(as.integer(chr))]
  if ("pos" %in% names(dt)) dt[, pos := as.integer(pos)]

  if (verbose) message("  After QC: ", nrow(dt), " rows retained")
  dt
}

# ── 3. QC summary ─────────────────────────────────────────────────────────────

qc_summary <- function(dt, label) {
  data.frame(
    dataset          = label,
    n_rows           = nrow(dt),
    n_unique_snp     = length(unique(dt$SNP)),
    n_duplicated_snp = sum(duplicated(dt$SNP)),
    p_min            = min(dt$pval, na.rm = TRUE),
    p_max            = max(dt$pval, na.rm = TRUE),
    missing_beta     = sum(is.na(dt$beta)),
    missing_se       = sum(is.na(dt$se)),
    missing_p        = sum(is.na(dt$pval))
  )
}

# ── 4. Format + clump exposure ────────────────────────────────────────────────

format_and_clump <- function(dt, cfg, p_thresh, instrument_cfg, save_dir = NULL) {
  dt_df <- as.data.frame(dt)

  fmt_args <- list(
    dat               = dt_df,
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    pval_col          = "pval",
    chr_col           = if ("chr" %in% names(dt_df)) "chr" else NULL,
    pos_col           = if ("pos" %in% names(dt_df)) "pos" else NULL,
    samplesize_col    = if ("samplesize" %in% names(dt_df)) "samplesize" else NULL
  )
  if (cfg$has_eaf && "eaf" %in% names(dt_df)) fmt_args$eaf_col <- "eaf"

  exp_fmt  <- do.call(TwoSampleMR::format_data, fmt_args)
  exp_inst <- exp_fmt[exp_fmt$pval.exposure < p_thresh, ]
  message("  Pre-clump instruments (p<", p_thresh, "): ", nrow(exp_inst))

  if (!is.null(save_dir)) {
    data.table::fwrite(exp_inst,
      file.path(save_dir, "tables", "01_exposure_instruments_preclump.csv"))
    saveRDS(exp_inst,
      file.path(save_dir, "rds", "01_exposure_instruments_preclump.rds"))
  }

  if (nrow(exp_inst) == 0) {
    warning("No instruments found at p<", p_thresh, " for ", cfg$id)
    return(list(preclump = exp_inst, clumped = exp_inst))
  }

  exp_clumped <- TwoSampleMR::clump_data(
    exp_inst,
    clump_kb = instrument_cfg$clump_kb,
    clump_r2 = instrument_cfg$clump_r2,
    pop      = cfg$ancestry
  )
  message("  Post-clump instruments: ", nrow(exp_clumped))

  if (!is.null(save_dir)) {
    data.table::fwrite(exp_clumped,
      file.path(save_dir, "tables", "02_exposure_instruments_clumped.csv"))
    saveRDS(exp_clumped,
      file.path(save_dir, "rds", "02_exposure_instruments_clumped.rds"))
  }

  list(preclump = exp_inst, clumped = exp_clumped)
}

# ── 5. Format outcome ─────────────────────────────────────────────────────────

format_outcome <- function(dt, cfg, snp_list, save_dir = NULL, suffix = "") {
  sub <- as.data.frame(dt[dt$SNP %in% snp_list, ])
  message("  Outcome rows matching instruments: ", nrow(sub))

  if (nrow(sub) == 0) {
    warning("No outcome SNPs match instruments for ", cfg$id)
    return(NULL)
  }

  fmt_args <- list(
    dat               = sub,
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    pval_col          = "pval",
    chr_col           = if ("chr" %in% names(sub)) "chr" else NULL,
    pos_col           = if ("pos" %in% names(sub)) "pos" else NULL,
    samplesize_col    = if ("samplesize" %in% names(sub)) "samplesize" else NULL
  )
  if (cfg$has_eaf && "eaf" %in% names(sub)) fmt_args$eaf_col <- "eaf"

  out_fmt <- do.call(TwoSampleMR::format_data, fmt_args)

  if (!is.null(save_dir)) {
    data.table::fwrite(out_fmt,
      file.path(save_dir, "tables", paste0("03_outcome_formatted", suffix, ".csv")))
    saveRDS(out_fmt,
      file.path(save_dir, "rds", paste0("03_outcome_formatted", suffix, ".rds")))
  }

  out_fmt
}

# ── 6. Harmonise + F-stats ────────────────────────────────────────────────────

harmonise_and_fstat <- function(exp_clumped, out_fmt, instrument_cfg,
                                exp_id, exp_label, out_id, out_label,
                                save_dir = NULL, suffix = "") {

  dat_h <- TwoSampleMR::harmonise_data(exp_clumped, out_fmt, action = 2)
  message("  Harmonised SNPs: ", nrow(dat_h))

  dat_h <- dat_h %>%
    mutate(
      F_stat      = (beta.exposure / se.exposure)^2,
      id.exposure = exp_id,
      exposure    = exp_label,
      id.outcome  = out_id,
      outcome     = out_label
    )

  f_summary <- dat_h %>%
    summarise(
      nsnp         = n(),
      min_F        = round(min(F_stat, na.rm = TRUE), 2),
      median_F     = round(median(F_stat, na.rm = TRUE), 2),
      mean_F       = round(mean(F_stat, na.rm = TRUE), 2),
      prop_F_lt_10 = mean(F_stat < 10, na.rm = TRUE)
    )
  message("  F-stat summary:"); print(f_summary)

  dat_strong <- dat_h %>% filter(F_stat >= instrument_cfg$F_min)
  message("  After F>=", instrument_cfg$F_min, ": ", nrow(dat_strong))

  if (!is.null(save_dir)) {
    data.table::fwrite(dat_h,
      file.path(save_dir, "tables", paste0("04_harmonised_with_F", suffix, ".csv")))
    data.table::fwrite(dat_strong,
      file.path(save_dir, "tables", paste0("05_harmonised_strongF", suffix, ".csv")))
    saveRDS(dat_h,
      file.path(save_dir, "rds", paste0("04_harmonised_with_F", suffix, ".rds")))
    saveRDS(dat_strong,
      file.path(save_dir, "rds", paste0("05_harmonised_strongF", suffix, ".rds")))
  }

  list(all = dat_h, strong = dat_strong, f_summary = f_summary)
}

# ── 7. Run MR + diagnostics ───────────────────────────────────────────────────

run_mr_full <- function(dat_strong, method_list, outcome_type = "binary",
                        save_dir = NULL, suffix = "") {

  if (nrow(dat_strong) == 0) {
    warning("No SNPs after F-filtering — skipping MR.")
    return(NULL)
  }

  if (nrow(dat_strong) == 1) {
    message("  Single instrument — Wald ratio only.")
    mr_res <- TwoSampleMR::mr(dat_strong, method_list = "mr_wald_ratio")
    if (!is.null(save_dir))
      data.table::fwrite(mr_res,
        file.path(save_dir, "tables", paste0("06_mr_wald_ratio", suffix, ".csv")))
    return(list(results = mr_res, raw = mr_res,
                het = NULL, pleio = NULL, loo = NULL, ss = NULL, steiger = NULL))
  }

  mr_res <- TwoSampleMR::mr(dat_strong, method_list = method_list)

  mr_res_out <- if (outcome_type == "binary") {
    TwoSampleMR::generate_odds_ratios(mr_res)
  } else {
    mr_res
  }

  het   <- TwoSampleMR::mr_heterogeneity(dat_strong)
  pleio <- TwoSampleMR::mr_pleiotropy_test(dat_strong)
  loo   <- TwoSampleMR::mr_leaveoneout(dat_strong)
  ss    <- TwoSampleMR::mr_singlesnp(dat_strong)

  steiger <- tryCatch(
    TwoSampleMR::steiger_filtering(dat_strong),
    error = function(e) { message("Steiger not run: ", e$message); NULL }
  )

  if (!is.null(save_dir)) {
    data.table::fwrite(mr_res_out,
      file.path(save_dir, "tables", paste0("06_mr_results", suffix, ".csv")))
    data.table::fwrite(het,
      file.path(save_dir, "tables", paste0("07_heterogeneity", suffix, ".csv")))
    data.table::fwrite(pleio,
      file.path(save_dir, "tables", paste0("08_pleiotropy", suffix, ".csv")))
    data.table::fwrite(loo,
      file.path(save_dir, "tables", paste0("09_leaveoneout", suffix, ".csv")))
    data.table::fwrite(ss,
      file.path(save_dir, "tables", paste0("10_singlesnp", suffix, ".csv")))
    if (!is.null(steiger))
      data.table::fwrite(steiger,
        file.path(save_dir, "tables", paste0("11_steiger", suffix, ".csv")))
    saveRDS(mr_res_out,
      file.path(save_dir, "rds", paste0("06_mr_results", suffix, ".rds")))
  }

  list(results = mr_res_out, raw = mr_res, het = het, pleio = pleio,
       loo = loo, ss = ss, steiger = steiger)
}

# ── 8. Plots ──────────────────────────────────────────────────────────────────

save_mr_plots <- function(mr_obj, dat_strong, save_dir, suffix = "") {
  if (is.null(mr_obj) || is.null(mr_obj$raw)) return(invisible(NULL))

  tryCatch({
    p_scatter <- TwoSampleMR::mr_scatter_plot(mr_obj$raw, dat_strong)[[1]]
    ggplot2::ggsave(
      file.path(save_dir, "plots", paste0("scatter", suffix, ".png")),
      p_scatter, width = 7, height = 5)
  }, error = function(e) message("scatter plot failed: ", e$message))

  if (!is.null(mr_obj$ss)) {
    tryCatch({
      p_forest <- TwoSampleMR::mr_forest_plot(mr_obj$ss)[[1]]
      ggplot2::ggsave(
        file.path(save_dir, "plots", paste0("forest", suffix, ".png")),
        p_forest, width = 7, height = 6)

      p_funnel <- TwoSampleMR::mr_funnel_plot(mr_obj$ss)[[1]]
      ggplot2::ggsave(
        file.path(save_dir, "plots", paste0("funnel", suffix, ".png")),
        p_funnel, width = 7, height = 6)
    }, error = function(e) message("forest/funnel failed: ", e$message))
  }

  if (!is.null(mr_obj$loo)) {
    tryCatch({
      p_loo <- TwoSampleMR::mr_leaveoneout_plot(mr_obj$loo)[[1]]
      ggplot2::ggsave(
        file.path(save_dir, "plots", paste0("leaveoneout", suffix, ".png")),
        p_loo, width = 7, height = 6)
    }, error = function(e) message("LOO plot failed: ", e$message))
  }

  invisible(NULL)
}

# ── 9. Attrition table ────────────────────────────────────────────────────────

make_attrition <- function(n_preclump, n_clumped, n_outcome_match,
                           n_harmonised, n_strong,
                           save_dir = NULL, suffix = "") {
  attrition <- data.frame(
    step  = c(
      "GWS instruments (pre-clump)",
      "Post-clump",
      "Present in outcome",
      "Harmonised",
      "Strong (F>=10)"
    ),
    n_snp = c(n_preclump, n_clumped, n_outcome_match, n_harmonised, n_strong)
  )
  if (!is.null(save_dir))
    data.table::fwrite(attrition,
      file.path(save_dir, "tables", paste0("00_attrition", suffix, ".csv")))
  attrition
}

# ── 10. MDE calculation ───────────────────────────────────────────────────────

calc_mde <- function(R2_total, N_eff, alpha = 0.05, power = 0.80,
                     outcome_type = "binary") {
  z_alpha  <- qnorm(1 - alpha / 2)
  z_power  <- qnorm(power)
  mde_beta <- (z_alpha + z_power) / sqrt(N_eff * R2_total)
  if (outcome_type == "binary") {
    list(mde_logOR = mde_beta, mde_OR = exp(mde_beta))
  } else {
    list(mde_beta = mde_beta)
  }
}
