library(ggplot2)
library(dplyr)
library(data.table)
library(here)

# ── Load MR results ────────────────────────────────────────────────────────────
res_root <- here::here("results", "forward")

pef_mr <- data.table::fread(
  file.path(res_root, "MVP_HFpEF_AFR", "tables", "06_mr_results.csv"))
ref_mr <- data.table::fread(
  file.path(res_root, "MVP_HFrEF_AFR", "tables", "06_mr_results.csv"))

# ── Build combined data frame ──────────────────────────────────────────────────
clean_method <- function(x) {
  x <- gsub("Inverse variance weighted \\(multiplicative random effects\\)",
             "IVW (random effects)", x)
  x <- gsub("Weighted median", "Weighted median", x)
  x <- gsub("MR Egger", "MR-Egger", x)
  x <- gsub("Weighted mode", "Weighted mode", x)
  x
}

# Method order — IVW at top (primary), then sensitivity
method_order <- c("IVW (random effects)", "Weighted median",
                  "MR-Egger", "Weighted mode")

prepare_mr <- function(dt, outcome_label) {
  dt %>%
    mutate(
      method  = clean_method(method),
      outcome = outcome_label
    ) %>%
    filter(method %in% method_order) %>%
    mutate(method = factor(method, levels = rev(method_order)))
}

pef_data <- prepare_mr(pef_mr, "HFpEF")
ref_data <- prepare_mr(ref_mr, "HFrEF")

# Combine — HFpEF on top, HFrEF below, with gap between
plot_data <- bind_rows(pef_data, ref_data) %>%
  mutate(
    outcome = factor(outcome, levels = c("HFpEF", "HFrEF")),
    # Primary method flag
    is_primary = method == "IVW (random effects)"
  )

# ── Colour scheme ──────────────────────────────────────────────────────────────
col_primary   <- "#1B4F8A"   # dark blue — IVW
col_secondary <- "#6C757D"   # grey — sensitivity

# ── Build plot ─────────────────────────────────────────────────────────────────
p <- ggplot(plot_data,
            aes(x = or, xmin = or_lci95, xmax = or_uci95,
                y = method,
                colour = is_primary,
                shape  = is_primary,
                size   = is_primary)) +

  # Facet by outcome — HFpEF top, HFrEF bottom
  facet_grid(outcome ~ ., scales = "free_y", space = "free_y",
             switch = "y") +

  # Null line
  geom_vline(xintercept = 1, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +

  # Confidence intervals
  geom_errorbarh(height = 0.25, linewidth = 0.7) +

  # Point estimates
  geom_point() +

  # Colour / shape / size scales
  scale_colour_manual(values = c("TRUE"  = col_primary,
                                  "FALSE" = col_secondary),
                      guide = "none") +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 17),
                     guide = "none") +
  scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3),
                    guide = "none") +

  # X axis — log scale for OR
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5),
    labels = c("0.10", "0.25", "0.50", "0.75", "1.00",
               "1.50", "2.00", "3.00", "5.00"),
    limits = c(0.05, 6)
  ) +

  # Labels
  labs(
    x       = "Odds ratio (95% CI) per SD increase in eGFRcrea",
    y       = NULL,
    title   = NULL,
    subtitle= NULL,
    caption = paste0(
      "Diamond = IVW random-effects (primary); Triangle = sensitivity methods.\n",
      "Estimates on log scale. Dashed line indicates null (OR = 1)."
    )
  ) +

  # Add OR text labels to the right of each point
  geom_text(
    aes(x = or_uci95 * 1.05,
        label = sprintf("%.2f (%.2f–%.2f)", or, or_lci95, or_uci95)),
    hjust = 0, size = 2.8,
    colour = "grey20",
    nudge_x = 0.01
  ) +

  # Theme
  theme_bw(base_size = 12) +
  theme(
    strip.text.y.left    = element_text(face = "bold", size = 12,
                                         angle = 0, hjust = 1),
    strip.placement      = "outside",
    strip.background     = element_rect(fill = "grey93", colour = "grey70"),
    panel.grid.major.y   = element_blank(),
    panel.grid.minor     = element_blank(),
    panel.spacing        = unit(0.8, "lines"),
    axis.text.y          = element_text(size = 10),
    axis.text.x          = element_text(size = 9),
    axis.title.x         = element_text(size = 10, margin = margin(t = 8)),
    plot.title           = element_text(face = "bold", size = 13),
    plot.subtitle        = element_text(size = 9, colour = "grey40"),
    plot.caption         = element_text(size = 8, colour = "grey50",
                                         hjust = 0),
    plot.margin          = margin(10, 120, 10, 10)   # extra right margin for labels
  )

# ── Save ───────────────────────────────────────────────────────────────────────
fig_dir <- here::here("results", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# TIFF for submission
ggsave(file.path(fig_dir, "Figure2_Combined_Forest_HFpEF_HFrEF.tiff"),
       p, width = 10, height = 6, dpi = 600, compression = "lzw")

# PNG for sharing/preprint
ggsave(file.path(fig_dir, "Figure2_Combined_Forest_HFpEF_HFrEF.png"),
       p, width = 10, height = 6, dpi = 600)

cat("Combined forest plot saved to:", fig_dir, "\n")
