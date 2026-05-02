# MVP_plots.R
# Regenerate all plots at 600 DPI from saved CSV files
# (RDS not required)

library(here)
library(TwoSampleMR)
library(ggplot2)
library(data.table)

res_root <- here::here("results", "forward")
plot_dir <- here::here("results", "figures")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

outcomes <- list(
  HFpEF = file.path(res_root, "MVP_HFpEF_AFR", "tables"),
  HFrEF = file.path(res_root, "MVP_HFrEF_AFR", "tables")
)

for (nm in names(outcomes)) {
  
  tbl_dir <- outcomes[[nm]]
  cat("Generating plots for", nm, "\n")
  
  # Load from CSVs
  dat_strong <- as.data.frame(data.table::fread(
    file.path(tbl_dir, "05_harmonised_strongF.csv")))
  mr_res     <- as.data.frame(data.table::fread(
    file.path(tbl_dir, "06_mr_results.csv")))
  ss         <- as.data.frame(data.table::fread(
    file.path(tbl_dir, "10_singlesnp.csv")))
  loo        <- as.data.frame(data.table::fread(
    file.path(tbl_dir, "09_leaveoneout.csv")))
  
  # Restore class attributes needed by TwoSampleMR plot functions
  class(mr_res) <- c("MRresult", "data.frame")
  
  # ── Forest plot ────────────────────────────────────────────────────────────
  p_forest <- TwoSampleMR::mr_forest_plot(ss)[[1]] +
    guides(colour = "none", shape = "none") +  # remove internal legend
    labs(
      title    = paste0("eGFR \u2192 ", nm),
      subtitle = "African ancestry (MVP)",
      caption  = "IVW = primary; sensitivity methods shown for comparison"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, colour = "grey40"),
      plot.caption  = element_text(size = 8, colour = "grey50"),
      legend.position = "none"   # belt and suspenders
    )
  
  ggsave(file.path(plot_dir, paste0("Figure_Forest_", nm, ".tiff")),
         p_forest, width = 7, height = 5, dpi = 600, compression = "lzw")
  ggsave(file.path(plot_dir, paste0("Figure_Forest_", nm, ".png")),
         p_forest, width = 7, height = 5, dpi = 600)
  
  # ── Scatter plot ───────────────────────────────────────────────────────────
  p_scatter <- TwoSampleMR::mr_scatter_plot(mr_res, dat_strong)[[1]] +
    labs(
      title   = paste0("eGFR \u2192 ", nm, " — Scatter"),
      caption = "African ancestry (MVP)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7)) +
    scale_x_continuous(expand = expansion(mult = 0.3))  # force x-axis to expand
  
  ggsave(file.path(plot_dir, paste0("FigureS_Scatter_", nm, ".tiff")),
         p_scatter, width = 7, height = 5, dpi = 600, compression = "lzw")
  ggsave(file.path(plot_dir, paste0("FigureS_Scatter_", nm, ".png")),
         p_scatter, width = 7, height = 5, dpi = 600)
  
  # ── Leave-one-out ──────────────────────────────────────────────────────────
  p_loo <- TwoSampleMR::mr_leaveoneout_plot(loo)[[1]] +
    guides(colour = "none", shape = "none") +
    labs(
      title   = paste0("eGFR \u2192 ", nm, " — Leave-one-out"),
      caption = "African ancestry (MVP)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  
  ggsave(file.path(plot_dir, paste0("FigureS_LOO_", nm, ".tiff")),
         p_loo, width = 7, height = 5, dpi = 600, compression = "lzw")
  ggsave(file.path(plot_dir, paste0("FigureS_LOO_", nm, ".png")),
         p_loo, width = 7, height = 5, dpi = 600)
  
  # ── Funnel plot ────────────────────────────────────────────────────────────
  p_funnel <- TwoSampleMR::mr_funnel_plot(ss)[[1]] +
    labs(title   = paste0("eGFR \u2192 ", nm, " — Funnel"),
         caption = "African ancestry (MVP)") +
    theme_bw(base_size = 12)
  
  ggsave(file.path(plot_dir, paste0("FigureS_Funnel_", nm, ".tiff")),
         p_funnel, width = 7, height = 5, dpi = 600, compression = "lzw")
  ggsave(file.path(plot_dir, paste0("FigureS_Funnel_", nm, ".png")),
         p_funnel, width = 7, height = 5, dpi = 600)
  
  cat("  Done:", nm, "\n")
}

cat("\nAll figures saved to:", plot_dir, "\n")