# =============================================================================
# MVP_supplementary_figures.R
# Generates Additional file 2: Supplementary Figures
# Three pages, paired side-by-side (HFpEF left, HFrEF right)
# Page 1: Scatter plots
# Page 2: Leave-one-out plots
# Page 3: Funnel plots
# Output: results/supplementary/Additional_file2_Supplementary_Figures.pdf
# =============================================================================

library(here)
library(data.table)
library(ggplot2)
library(TwoSampleMR)
library(patchwork)

res_root <- here::here("results", "forward")
out_dir  <- here::here("results", "supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data from CSVs ───────────────────────────────────────────────────────
outcomes <- list(
  HFpEF = file.path(res_root, "MVP_HFpEF_AFR", "tables"),
  HFrEF = file.path(res_root, "MVP_HFrEF_AFR", "tables")
)

load_outcome_data <- function(tbl_dir) {
  list(
    dat_strong = as.data.frame(data.table::fread(
      file.path(tbl_dir, "05_harmonised_strongF.csv"))),
    mr_res     = as.data.frame(data.table::fread(
      file.path(tbl_dir, "06_mr_results.csv"))),
    ss         = as.data.frame(data.table::fread(
      file.path(tbl_dir, "10_singlesnp.csv"))),
    loo        = as.data.frame(data.table::fread(
      file.path(tbl_dir, "09_leaveoneout.csv")))
  )
}

pef_data <- load_outcome_data(outcomes$HFpEF)
ref_data <- load_outcome_data(outcomes$HFrEF)

# ── Theme for supplementary plots ─────────────────────────────────────────────
supp_theme <- theme_bw(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    legend.position = "none",
    plot.margin   = margin(8, 8, 8, 8)
  )

# ── Helper: clean up TwoSampleMR plots ───────────────────────────────────────
clean_plot <- function(p, title, subtitle = NULL) {
  p +
    guides(colour = "none", shape = "none", fill = "none") +
    labs(title = title, subtitle = subtitle) +
    supp_theme
}

# ── PAGE 1: Scatter plots ──────────────────────────────────────────────────────
cat("Building Page 1: Scatter plots\n")

p_scatter_pef <- clean_plot(
  TwoSampleMR::mr_scatter_plot(pef_data$mr_res, pef_data$dat_strong)[[1]] +
    scale_x_continuous(expand = expansion(mult = 0.3)),
  title    = "Figure S1a. eGFR \u2192 HFpEF",
  subtitle = "Scatter plot"
)

p_scatter_ref <- clean_plot(
  TwoSampleMR::mr_scatter_plot(ref_data$mr_res, ref_data$dat_strong)[[1]] +
    scale_x_continuous(expand = expansion(mult = 0.3)),
  title    = "Figure S1b. eGFR \u2192 HFrEF",
  subtitle = "Scatter plot"
)

page1 <- p_scatter_pef + p_scatter_ref +
  plot_layout(ncol = 2) +
  plot_annotation(
    title   = "Figure S1. Scatter plots of SNP effects on eGFR and heart failure subtypes",
    caption = paste0(
      "Each point represents a single SNP instrument. Lines show MR method estimates.\n",
      "IVW = inverse-variance weighted (primary); sensitivity methods shown for comparison.\n",
      "HFpEF, heart failure with preserved ejection fraction; ",
      "HFrEF, heart failure with reduced ejection fraction."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 8, colour = "grey40", hjust = 0)
    )
  )

# ── PAGE 2: Leave-one-out plots ────────────────────────────────────────────────
cat("Building Page 2: Leave-one-out plots\n")

p_loo_pef <- clean_plot(
  TwoSampleMR::mr_leaveoneout_plot(pef_data$loo)[[1]] +
    guides(colour = "none", shape = "none"),
  title    = "Figure S2a. eGFR \u2192 HFpEF",
  subtitle = "Leave-one-out analysis"
)

p_loo_ref <- clean_plot(
  TwoSampleMR::mr_leaveoneout_plot(ref_data$loo)[[1]] +
    guides(colour = "none", shape = "none"),
  title    = "Figure S2b. eGFR \u2192 HFrEF",
  subtitle = "Leave-one-out analysis"
)

page2 <- p_loo_pef + p_loo_ref +
  plot_layout(ncol = 2) +
  plot_annotation(
    title   = "Figure S2. Leave-one-out sensitivity analyses",
    caption = paste0(
      "Each row shows the IVW estimate after sequentially removing one SNP instrument.\n",
      "The bottom row (All) shows the overall IVW estimate using all instruments.\n",
      "Stability of estimates across rows indicates no single SNP drives the result.\n",
      "HFpEF, heart failure with preserved ejection fraction; ",
      "HFrEF, heart failure with reduced ejection fraction."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 8, colour = "grey40", hjust = 0)
    )
  )

# ── PAGE 3: Funnel plots ───────────────────────────────────────────────────────
cat("Building Page 3: Funnel plots\n")

p_funnel_pef <- clean_plot(
  TwoSampleMR::mr_funnel_plot(pef_data$ss)[[1]],
  title    = "Figure S3a. eGFR \u2192 HFpEF",
  subtitle = "Funnel plot"
)

p_funnel_ref <- clean_plot(
  TwoSampleMR::mr_funnel_plot(ref_data$ss)[[1]],
  title    = "Figure S3b. eGFR \u2192 HFrEF",
  subtitle = "Funnel plot"
)

page3 <- p_funnel_pef + p_funnel_ref +
  plot_layout(ncol = 2) +
  plot_annotation(
    title   = "Figure S3. Funnel plots of single-SNP MR estimates",
    caption = paste0(
      "Each point represents one SNP instrument. Vertical lines show IVW (light blue) ",
      "and MR-Egger (dark blue) pooled estimates.\n",
      "Symmetry around the pooled estimate suggests absence of directional horizontal pleiotropy.\n",
      "HFpEF, heart failure with preserved ejection fraction; ",
      "HFrEF, heart failure with reduced ejection fraction."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 8, colour = "grey40", hjust = 0)
    )
  )

# ── Save as multi-page PDF ─────────────────────────────────────────────────────
cat("Saving supplementary figures PDF...\n")

out_path <- file.path(out_dir, "Additional_file2_Supplementary_Figures.pdf")

pdf(out_path, width = 11, height = 7, paper = "a4r")   # landscape A4
print(page1)
print(page2)
print(page3)
dev.off()

cat("Saved to:", out_path, "\n")
cat("Pages: 3 (Scatter | Leave-one-out | Funnel)\n")
