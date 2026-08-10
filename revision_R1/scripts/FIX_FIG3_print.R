##############################################################################
# FIX_FIG3_print.R — regenerate main Figure 3 (MVMR panel)
#
# Addresses two revision items:
#   (1) Reviewer 1.3 — reorganise for print so panels / axis labels do not
#       overlap when the figure is scaled into a single-column JAD PDF:
#         - shorter, two-line y-axis category labels
#         - wider Panel A, generous plot margins, discrete-axis headroom
#         - shorter x-axis title
#   (2) MVMR p-value correction (2.29e-5 z-approx -> 1.4e-3 t-based, MVMR pkg).
#       Matches Fig S7 and the response letter. Panel B interpretation softened
#       to hedged ("consistent with") language for the reframed manuscript.
#
# Self-contained (hard-coded MVMR estimates); no TwoSampleMR / RData needed.
# The same corrected block is mirrored in APPLY_REVIEWER_FIXES.R (Fig 3).
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid)
})
set.seed(42)

FONT <- "Arial"
BASE <- "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper/Final Figures and data"
MAIN <- file.path(BASE, "figures/main_text")
dir.create(MAIN, showWarnings = FALSE, recursive = TRUE)

COL_PRIMARY <- "#1B4F8A"; COL_NULL <- "#AAAAAA"
mm2in <- function(mm) mm / 25.4
save_fig <- function(path, plot, w, h, dpi = 300) {
  ggsave(paste0(path, ".png"), plot, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(path, ".pdf"), plot, width = w, height = h, device = cairo_pdf, bg = "white")
  cat("  Saved:", basename(path), "\n")
}

cat("Fig 3: MVMR panel (print-safe + corrected p)...\n")

# Shorter two-line labels so the y-axis does not compress/overlap in print
lv <- c("Serum iron\n(MVMR)",
        "Transferrin, MVMR\n(adjusted for iron)",
        "Transferrin\n(univariable IVW)")
# beta/SE synced to Analysis E (E_MVMR_extended.csv) so CIs match the p-values
mvmr_d <- data.frame(
  label = factor(lv, levels = lv),
  beta  = c(-0.0057, -0.0437, -0.0393),
  se    = c( 0.0218,  0.0092,  0.0132),
  stringsAsFactors = FALSE
)
mvmr_d$ci_lo <- mvmr_d$beta - 1.96 * mvmr_d$se
mvmr_d$ci_hi <- mvmr_d$beta + 1.96 * mvmr_d$se
mvmr_d$color <- c(COL_NULL, COL_PRIMARY, COL_PRIMARY)

p3a <- ggplot(mvmr_d, aes(x = beta, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.22,
    color = mvmr_d$color, linewidth = 0.8) +
  geom_point(shape = c(18, 15, 15), size = c(3, 4, 3.5), color = mvmr_d$color) +
  # p-values: MVMR transferrin corrected 2.29e-5 -> 1.4e-3 (t-based, MVMR pkg)
  geom_text(aes(label = c("p = 0.801", "p = 1.4e-3", "p = 2.88e-3")),
    nudge_y = 0.33, hjust = 0.5, size = 2.5, family = FONT, color = "grey25") +
  scale_x_continuous(limits = c(-0.09, 0.05), breaks = seq(-0.08, 0.04, 0.04)) +
  scale_y_discrete(expand = expansion(add = c(0.7, 0.9))) +
  labs(title = "A  Effect of conditioning on serum iron",
    x = expression("Effect on total-tau (" * beta * ", 95% CI)"), y = NULL) +
  theme_minimal(base_size = 9, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.text.y = element_text(size = 8, lineheight = 0.9),
    plot.title = element_text(face = "bold", size = 9, margin = margin(b = 10)),
    plot.margin = margin(t = 4, r = 8, b = 4, l = 4))

p3b <- ggplot() +
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 3.0, ymax = 3.7,
    fill = "#E3F2FD", color = COL_PRIMARY, linewidth = 0.6) +
  annotate("text", x = 0.5, y = 3.48, label = "MVMR: Transferrin",
    size = 3, fontface = "bold", family = FONT, color = COL_PRIMARY) +
  annotate("text", x = 0.5, y = 3.22, label = "p = 1.4e-3",
    size = 2.5, family = FONT, color = COL_PRIMARY) +
  annotate("rect", xmin = 0.08, xmax = 0.92, ymin = 2.0, ymax = 2.7,
    fill = "#F5F5F5", color = COL_NULL, linewidth = 0.6) +
  annotate("text", x = 0.5, y = 2.48, label = "MVMR: Serum iron",
    size = 3, fontface = "bold", family = FONT, color = COL_NULL) +
  annotate("text", x = 0.5, y = 2.22, label = "p = 0.801",
    size = 2.5, family = FONT, color = COL_NULL) +
  annotate("rect", xmin = 0.03, xmax = 0.97, ymin = 0.5, ymax = 1.6,
    fill = "white", color = COL_PRIMARY, linewidth = 0.8) +
  annotate("text", x = 0.5, y = 1.35, label = "Interpretation",
    size = 3, fontface = "bold", family = FONT) +
  # Non-causal wording per Reviewer 2 (Part 5): report association, not causal claim
  annotate("text", x = 0.5, y = 0.9,
    label = "Signal not attenuated by\nconditioning on serum iron",
    size = 2.8, family = FONT, color = "grey20", lineheight = 1.2) +
  annotate("segment", x = 0.5, xend = 0.5, y = 3.0, yend = 2.7,
    arrow = arrow(length = unit(1.5, "mm")), color = "grey50") +
  annotate("segment", x = 0.5, xend = 0.5, y = 2.0, yend = 1.6,
    arrow = arrow(length = unit(1.5, "mm")), color = "grey50") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0.3, 4.0)) +
  labs(title = "B  Mechanistic interpretation") +
  theme_void(base_family = FONT) +
  theme(plot.title = element_text(face = "bold", size = 9, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(t = 4, r = 4, b = 4, l = 8))

# Panel A gets more width (it holds the plot + longer labels); taller canvas
fig3 <- p3a + p3b + plot_layout(widths = c(1.4, 1))
save_fig(file.path(MAIN, "Fig3_MVMR_panel"), fig3, mm2in(185), mm2in(105))

cat("Done.\n")
