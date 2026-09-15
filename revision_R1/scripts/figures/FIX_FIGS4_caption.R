##############################################################################
# FIX_FIGS4_caption.R
# Re-renders ONLY Supplementary Figure S4 with a corrected footnote.
#
# The FigS4 block is copied verbatim from APPLY_REVIEWER_FIXES.R (lines 462-551);
# the only change is the caption, which wrongly said the p-tau and Abeta42
# panels came from "different GWAS discovery cohorts". Both come from the same
# study (Jansen et al. 2022) with different discovery sample sizes, matching the
# S4 legend corrected on 2026-08-18.
#
# Do NOT run APPLY_REVIEWER_FIXES.R / APPLY_MCH.R to regenerate S4: they also
# regenerate Figure 3 from a superseded code path (retired serum-iron p=0.817).
#
# Usage: Rscript -e ".libPaths(...); source('FIX_FIGS4_caption.R')" <out_dir>
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[1] else
  "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper/Final Figures and data/figures/supplementary"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FONT <- "Arial"
COL_PRIMARY <- "#1B4F8A"
COL_SECONDARY <- "#4A90D9"
COL_NULL <- "#AAAAAA"
COL_RISK <- "#C62828"
mm2in <- function(mm) mm / 25.4

save_fig <- function(path, plot, w, h, dpi = 300) {
  ggsave(paste0(path, ".png"), plot, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(path, ".pdf"), plot, width = w, height = h, device = cairo_pdf, bg = "white")
  cat("  Saved:", basename(path), "\n")
}

cat("FigS4: CSF subgroups...\n")

unified_xlim <- c(-0.25, 0.27)

ptau <- data.frame(
  label = c("Main (n=7,798)", "APOE4 non-carriers (n=3,141)",
            "APOE4 carriers (n=3,047)", "Normal amyloid (n=3,174)",
            "Abnormal amyloid (n=3,534)"),
  beta = c(-0.050, -0.116, 0.016, -0.075, -0.018),
  se = c(0.033, 0.059, 0.053, 0.057, 0.049),
  pt_color = c(COL_SECONDARY, COL_PRIMARY, COL_RISK, COL_SECONDARY, COL_SECONDARY),
  p = c(0.134, 0.048, 0.769, 0.185, 0.713), stringsAsFactors = FALSE
)
ptau$ci_lo <- ptau$beta - 1.96 * ptau$se
ptau$ci_hi <- ptau$beta + 1.96 * ptau$se
ptau$row <- rev(seq_len(nrow(ptau)))
ptau$p_lab <- ifelse(ptau$p < 0.05, paste0(sprintf("%.3f", ptau$p), "*"), sprintf("%.3f", ptau$p))

p4a <- ggplot(ptau, aes(x = beta, y = row)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 3.55, ymax = 4.45, fill = "#E8F5E9") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2.55, ymax = 3.45, fill = "#FFEBEE") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0,
    color = ptau$pt_color, linewidth = 0.6) +
  geom_point(color = ptau$pt_color, size = 3, shape = 16) +
  geom_text(aes(label = p_lab), x = 0.20, hjust = 0, size = 2.5, family = FONT,
    fontface = ifelse(ptau$p < 0.05, "bold", "plain"),
    color = ifelse(ptau$p < 0.05, COL_PRIMARY, "grey40")) +
  scale_y_continuous(breaks = ptau$row, labels = ptau$label,
    expand = expansion(add = c(0.5, 0.5))) +
  scale_x_continuous(limits = unified_xlim) +
  labs(title = "A  CSF phosphorylated tau (p-tau)",
    x = expression("Effect on CSF p-tau (IVW " * beta * ", 95% CI)"), y = NULL) +
  theme_minimal(base_size = 9, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3))

ab <- data.frame(
  label = c("Main (n=8,074)", "APOE4 non-carriers (n=3,201)",
            "APOE4 carriers (n=3,240)", "Normal amyloid (n=3,182)",
            "Abnormal amyloid (n=3,775)"),
  beta = c(0.012, 0.002, -0.002, 0.008, -0.009),
  se = c(0.034, 0.052, 0.051, 0.052, 0.055),
  p = c(0.719, 0.977, 0.967, 0.879, 0.875),
  stringsAsFactors = FALSE
)
ab$ci_lo <- ab$beta - 1.96 * ab$se
ab$ci_hi <- ab$beta + 1.96 * ab$se
ab$row <- rev(seq_len(nrow(ab)))

p4b <- ggplot(ab, aes(x = beta, y = row)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0,
    color = COL_NULL, linewidth = 0.6) +
  geom_point(color = COL_NULL, size = 3, shape = 18) +
  geom_text(aes(label = sprintf("%.3f", p)), x = 0.20, hjust = 0, size = 2.5,
    family = FONT, color = "grey40") +
  scale_y_continuous(breaks = ab$row, labels = ab$label,
    expand = expansion(add = c(0.5, 0.5))) +
  scale_x_continuous(limits = unified_xlim) +
  labs(title = expression("B  CSF amyloid-" * beta * "-42"),
    x = expression("Effect on CSF A" * beta * "42 (IVW " * beta * ", 95% CI)"), y = NULL) +
  theme_minimal(base_size = 9, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3))

figs4 <- p4a / p4b + plot_layout(heights = c(1.5, 1.3)) +
  plot_annotation(
    caption = paste0(
      "* p < 0.05 (nominal significance). ",
      ### FIX (2026-09-13): same study, different discovery sample sizes -- not different cohorts
      "Sample sizes differ between panels A and B because the p-tau and A\u03b242 GWAS\n",
      "(Jansen et al. 2022) have different discovery sample sizes."
    ),
    theme = theme(
      ### FIX (2026-09-13): 6 pt fell to ~4.7 pt once this page is scaled into the
      ### combined supplementary PDF; 7.5 pt keeps it ~5.9 pt at final size
      plot.caption = element_text(size = 7.5, color = "grey45", hjust = 0,
        family = FONT, margin = margin(t = 4))
    )
  )

save_fig(file.path(OUT, "FigS4_CSF_ptau_subgroups"), figs4, mm2in(165), mm2in(165))
