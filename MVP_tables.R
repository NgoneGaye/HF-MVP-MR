# =============================================================================
# MVP_tables.R
# Generates publication-ready tables for BMC Nephrology submission
# Output: results/tables/MVP_MR_Tables.docx
#
# Run after MVP_forward.R has completed successfully.
# =============================================================================

library(here)
library(data.table)
library(dplyr)
library(flextable)
library(officer)

# ── Paths ─────────────────────────────────────────────────────────────────────
res_root  <- here::here("results", "forward")
out_dir   <- here::here("results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load all result files ─────────────────────────────────────────────────────

# Instruments (same for both outcomes — use HFpEF folder)
instruments <- data.table::fread(
  file.path(res_root, "eGFR_instruments", "tables", "02_exposure_instruments_clumped.csv"))

# MR results
pef_mr  <- data.table::fread(file.path(res_root, "MVP_HFpEF_AFR", "tables", "06_mr_results.csv"))
ref_mr  <- data.table::fread(file.path(res_root, "MVP_HFrEF_AFR", "tables", "06_mr_results.csv"))

# Heterogeneity
pef_het <- data.table::fread(file.path(res_root, "MVP_HFpEF_AFR", "tables", "07_heterogeneity.csv"))
ref_het <- data.table::fread(file.path(res_root, "MVP_HFrEF_AFR", "tables", "07_heterogeneity.csv"))

# Pleiotropy
pef_pleo <- data.table::fread(file.path(res_root, "MVP_HFpEF_AFR", "tables", "08_pleiotropy.csv"))
ref_pleo <- data.table::fread(file.path(res_root, "MVP_HFrEF_AFR", "tables", "08_pleiotropy.csv"))

# Attrition
pef_att <- data.table::fread(file.path(res_root, "MVP_HFpEF_AFR", "tables", "00_attrition.csv"))
ref_att <- data.table::fread(file.path(res_root, "MVP_HFrEF_AFR", "tables", "00_attrition.csv"))

# Power/MDE
pef_mde <- data.table::fread(file.path(res_root, "MVP_HFpEF_AFR", "tables", "00_power_MDE.csv"))
ref_mde <- data.table::fread(file.path(res_root, "MVP_HFrEF_AFR", "tables", "00_power_MDE.csv"))

# ── Helper: format p-values ───────────────────────────────────────────────────
fmt_p <- function(p) {
  ifelse(p < 0.001, "<0.001",
  ifelse(p < 0.01,  formatC(p, digits = 3, format = "f"),
                    formatC(p, digits = 3, format = "f")))
}

fmt_n <- function(x) formatC(round(x, 3), digits = 3, format = "f")

# ── TABLE 1: Instrument selection and strength ────────────────────────────────
cat("Building Table 1 — Instrument selection and strength\n")

# Compute F-stat from harmonised file (use HFpEF as reference)
harm <- data.table::fread(
  file.path(res_root, "MVP_HFpEF_AFR", "tables", "05_harmonised_strongF.csv"))

t1_data <- harm %>%
  mutate(F_stat = round((beta.exposure / se.exposure)^2, 1)) %>%
  select(SNP, chr.exposure, pos.exposure,
         effect_allele.exposure, other_allele.exposure,
         eaf.exposure, beta.exposure, se.exposure,
         pval.exposure, F_stat) %>%
  arrange(chr.exposure, pos.exposure) %>%
  mutate(
    eaf.exposure   = round(eaf.exposure, 3),
    beta.exposure  = round(beta.exposure, 4),
    se.exposure    = round(se.exposure, 4),
    pval.exposure  = formatC(pval.exposure, format = "e", digits = 2)
  )

# Add summary row
r2_total <- round(sum(harm$F_stat / (harm$F_stat + harm$samplesize.exposure - 2),
                      na.rm = TRUE), 5)

t1_display <- t1_data %>%
  rename(
    "rsID"         = SNP,
    "Chr"          = chr.exposure,
    "Position"     = pos.exposure,
    "EA"           = effect_allele.exposure,
    "OA"           = other_allele.exposure,
    "EAF"          = eaf.exposure,
    "β"            = beta.exposure,
    "SE"           = se.exposure,
    "P-value"      = pval.exposure,
    "F-statistic"  = F_stat
  )

ft1 <- flextable(t1_display) %>%
  set_caption(paste0(
    "Table 1. Genetic instruments for estimated glomerular filtration rate (eGFRcrea) ",
    "in African ancestry populations. SNPs were selected at genome-wide significance ",
    "(p < 5×10\u207B\u2078) and clumped for independence (r\u00B2 < 0.001; 10 Mb window). ",
    "Total variance explained R\u00B2 = ", r2_total, ". ",
    "EA, effect allele; OA, other allele; EAF, effect allele frequency; ",
    "\u03B2, effect estimate (SD units of eGFR); SE, standard error.")) %>%
  bold(part = "header") %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit()

# ── TABLE 2: Main MR results ──────────────────────────────────────────────────
cat("Building Table 2 — MR results\n")

# Method labels — clean up TwoSampleMR verbose names
clean_method <- function(x) {
  x <- gsub("Inverse variance weighted \\(multiplicative random effects\\)",
             "IVW (random effects)", x)
  x <- gsub("Weighted median", "Weighted median", x)
  x <- gsub("MR Egger", "MR-Egger", x)
  x <- gsub("Weighted mode", "Weighted mode", x)
  x
}

build_mr_table <- function(mr_res, het_res, pleo_res, outcome_label) {

  # IVW Q p-value
  ivw_q_p <- het_res %>%
    filter(grepl("inverse variance", method, ignore.case = TRUE)) %>%
    pull(Q_pval) %>% round(3)

  # Egger intercept p
  egger_p <- pleo_res$pval %>% round(3)

  mr_res %>%
    mutate(
      Outcome  = outcome_label,
      Method   = clean_method(method),
      SNPs     = nsnp,
      OR       = round(or, 2),
      CI       = paste0(round(or_lci95, 2), "–", round(or_uci95, 2)),
      P        = fmt_p(pval),
      Het_Q_p  = ifelse(Method == "IVW (random effects)",
                        fmt_p(ivw_q_p), ""),
      Egger_p  = ifelse(Method == "MR-Egger",
                        fmt_p(egger_p), "")
    ) %>%
    select(Outcome, Method, SNPs, OR, CI, P, Het_Q_p, Egger_p)
}

t2_pef <- build_mr_table(pef_mr, pef_het, pef_pleo, "HFpEF")
t2_ref <- build_mr_table(ref_mr, ref_het, ref_pleo, "HFrEF")
t2_data <- bind_rows(t2_pef, t2_ref)

# Clean column names for display
t2_display <- t2_data %>%
  rename(
    "Outcome"            = Outcome,
    "Method"             = Method,
    "SNPs (n)"           = SNPs,
    "OR"                 = OR,
    "95% CI"             = CI,
    "P-value"            = P,
    "Heterogeneity Q\np" = Het_Q_p,
    "Egger intercept\np" = Egger_p
  )

ft2 <- flextable(t2_display) %>%
  set_caption(
    "Table 2. Mendelian randomisation estimates for the effect of genetically predicted eGFRcrea on heart failure subtypes in African ancestry populations. Primary method: IVW random-effects. OR, odds ratio per SD increase in genetically predicted eGFRcrea; CI, confidence interval; IVW, inverse-variance weighted; HFpEF, heart failure with preserved ejection fraction; HFrEF, heart failure with reduced ejection fraction. Heterogeneity Q p-value reported for IVW only; MR-Egger intercept p-value reported for MR-Egger only.") %>%
  bold(part = "header") %>%
  bold(i = ~ Method == "IVW (random effects)") %>%
  merge_v(j = "Outcome") %>%
  hline(i = nrow(t2_pef), border = fp_border(width = 0.75, style = "dashed")) %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit()

# ── TABLE S1: Attrition ───────────────────────────────────────────────────────
cat("Building Table S1 — Attrition\n")

ts1_data <- pef_att %>%
  rename("Step" = step, "HFpEF (n SNPs)" = n_snp) %>%
  left_join(ref_att %>% rename("Step" = step, "HFrEF (n SNPs)" = n_snp),
            by = "Step")

fts1 <- flextable(ts1_data) %>%
  set_caption(
    "Additional file 1: Table S1. SNP attrition from instrument selection to analysis. GWS, genome-wide significant (p < 5×10\u207B\u2078); AFR, African ancestry linkage disequilibrium reference panel.") %>%
  bold(part = "header") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit()

# ── TABLE S2: Power/MDE ───────────────────────────────────────────────────────
cat("Building Table S2 — Power/MDE\n")

ts2_data <- bind_rows(pef_mde, ref_mde) %>%
  mutate(
    R2_total  = round(R2_total, 5),
    N_eff     = formatC(N_eff, format = "d", big.mark = ","),
    MDE_logOR = round(MDE_logOR, 3),
    MDE_OR    = round(MDE_OR, 3)
  ) %>%
  select(outcome, R2_total, N_eff, MDE_logOR, MDE_OR) %>%
  rename(
    "Outcome"         = outcome,
    "R\u00B2 (total)" = R2_total,
    "N\u1D49\u1DA0\u1DA0"    = N_eff,
    "MDE (log-OR)"    = MDE_logOR,
    "MDE (OR)"        = MDE_OR
  )

fts2 <- flextable(ts2_data) %>%
  set_caption(
    "Additional file 2: Table S2. Statistical power and minimum detectable effects (MDE) at 80% power and \u03B1 = 0.05. R\u00B2, total variance in eGFRcrea explained by retained instruments; N\u1D49\u1DA0\u1DA0, effective sample size computed as 4/(1/N\u2090\u2090\u2090\u2090 + 1/N\u209C\u2092\u2099\u209C\u2B63\u2B63\u2B63\u2B63); MDE, minimum detectable effect.") %>%
  bold(part = "header") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit()

# ── TABLE S3: Heterogeneity and pleiotropy ────────────────────────────────────
cat("Building Table S3 — Diagnostics\n")

diag_pef <- bind_cols(
  pef_het %>% filter(grepl("inverse variance", method, ignore.case = TRUE)) %>%
    select(Q, Q_df, Q_pval) %>%
    mutate(across(where(is.numeric), ~round(.x, 3))),
  pef_pleo %>%
    select(egger_intercept, se, pval) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)))
) %>% mutate(Outcome = "HFpEF")

diag_ref <- bind_cols(
  ref_het %>% filter(grepl("inverse variance", method, ignore.case = TRUE)) %>%
    select(Q, Q_df, Q_pval) %>%
    mutate(across(where(is.numeric), ~round(.x, 3))),
  ref_pleo %>%
    select(egger_intercept, se, pval) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)))
) %>% mutate(Outcome = "HFrEF")

ts3_data <- bind_rows(diag_pef, diag_ref) %>%
  select(Outcome, Q, Q_df, Q_pval, egger_intercept, se, pval) %>%
  rename(
    "Outcome"              = Outcome,
    "Cochran Q"            = Q,
    "df"                   = Q_df,
    "Q p-value"            = Q_pval,
    "Egger intercept"      = egger_intercept,
    "SE"                   = se,
    "Intercept p-value"    = pval
  )

fts3 <- flextable(ts3_data) %>%
  set_caption(
    "Additional file 3: Table S3. Heterogeneity and pleiotropy diagnostics. Cochran Q statistic from IVW model; MR-Egger intercept test for directional horizontal pleiotropy. SE, standard error of Egger intercept.") %>%
  bold(part = "header") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1.5)) %>%
  hline_bottom(part = "header", border = fp_border(width = 1)) %>%
  hline_bottom(part = "body", border = fp_border(width = 1.5)) %>%
  fontsize(size = 10, part = "all") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  autofit()

# ── Write Word document ───────────────────────────────────────────────────────
cat("Writing Word document...\n")

doc <- read_docx() %>%

  # Title
  body_add_par("Tables — MVP MR: eGFR → HFpEF / HFrEF",
               style = "heading 1") %>%
  body_add_par("Main text tables", style = "heading 2") %>%

  # Table 1
  body_add_flextable(ft1) %>%
  body_add_par("") %>%
  body_add_par("") %>%

  # Table 2
  body_add_flextable(ft2) %>%
  body_add_par("") %>%

  # Page break before supplementary
  body_add_break() %>%
  body_add_par("Supplementary tables", style = "heading 2") %>%
  body_add_par("") %>%

  # Table S1
  body_add_flextable(fts1) %>%
  body_add_par("") %>%
  body_add_par("") %>%

  # Table S2
  body_add_flextable(fts2) %>%
  body_add_par("") %>%
  body_add_par("") %>%

  # Table S3
  body_add_flextable(fts3)

out_path <- file.path(out_dir, "MVP_MR_Tables.docx")
print(doc, target = out_path)
cat("Tables saved to:", out_path, "\n")
