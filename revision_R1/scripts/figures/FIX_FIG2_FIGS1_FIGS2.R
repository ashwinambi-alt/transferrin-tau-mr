##############################################################################
# FIX_FIG2_FIGS1_FIGS2.R  (2026-09-13)
# Re-renders ONLY Figure 2, Supplementary Figure S1 and Supplementary Figure S2.
# Blocks copied from scripts/APPLY_REVIEWER_FIXES.R (Fig 2: 158-229; FigS1:
# 310-377; FigS2: 380-412). Do NOT run that script as a whole: it also rewrites
# Figure 3 from a superseded code path (retired serum-iron p = 0.817).
#
# Fixes applied here (audit, 2026-09-13) — no numeric value changes:
#   1. Fig 2 / FigS1 annotations printed a Latin "b" for the effect estimate and
#      an ASCII hyphen; now β and a true minus (U+2212).
#   2. Fig 2 / FigS2 printed "I2"; now I².
#   3. Fig 2 printed "p = 2.88e-3"; now 2.88×10⁻³ to match the tables and
#      legends. Arial has NO superscript minus (U+207B), so the exponent is set
#      with plotmath and a true minus, as in FIX_FIGS5_S7_vector.R. The Fig 2
#      statistics box is therefore drawn as a rect + one text layer per line
#      (a single multi-line label cannot mix plotmath and plain text).
#   4. FigS1 bolded "rs1495741 (NAT2)" instead of "All SNPs": the face vector was
#      indexed loo$is_all[order(loo$row)], but vectorised element_text() is applied
#      in the order the breaks are supplied (loo$row). Same bug as Fig 1.
#   5. Label overlaps: Fig 2 "TFRC" sat on the IVW line; FigS2 "HLA*"/"NAT2" sat
#      on their points / the dashed line. Nudged clear (ggrepel).
#
# Usage: Rscript -e ".libPaths(...); source('FIX_FIG2_FIGS1_FIGS2.R')" <project_dir> <out_dir>
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(grid); library(TwoSampleMR)
})
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- args[1]
OUT <- if (length(args) >= 2) args[2] else
  file.path(PROJECT_DIR, "Final Figures and data", "figures", "_staging_F2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FONT <- "Arial"
COL_PRIMARY <- "#1B4F8A"
COL_SECONDARY <- "#4A90D9"
COL_RISK <- "#C62828"
MINUS <- "−"
mm2in <- function(mm) mm / 25.4

save_fig <- function(path, plot, w, h, dpi = 300) {
  ggsave(paste0(path, ".png"), plot, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(path, ".pdf"), plot, width = w, height = h, device = cairo_pdf, bg = "white")
  cat("  Saved:", basename(path), "\n")
}

PLEIO_FOOTNOTE <- "* Instruments in known pleiotropic loci; excluded in 6-SNP sensitivity analysis."

load(file.path(PROJECT_DIR, "data", "FINAL_analysis_results.RData"))

gene_labels <- c(rs8177240 = "TF", rs1800562 = "HFE", rs1495741 = "NAT2",
  rs9990333 = "TFRC", rs174577 = "FADS*", rs744653 = "SLC40A1",
  rs7646473 = "TF (near)", rs9268633 = "HLA*")

# ============================================================
# FIGURE 2 — scatter
# ============================================================

cat("Fig 2: Scatter...\n")

dat_tt <- all_harmonised[["transferrin_circ_total_tau"]]
mr_tt <- all_results[["transferrin_circ_total_tau"]]
dat_tt$F_stat <- (dat_tt$beta.exposure / dat_tt$se.exposure)^2
dat_tt$gene <- gene_labels[dat_tt$SNP]
dat_tt$pleio <- dat_tt$SNP %in% c("rs174577", "rs9268633")

ivw_b <- mr_tt$b[mr_tt$method == "Inverse variance weighted"]
wm_b <- mr_tt$b[mr_tt$method == "Weighted median"]

fig2_base <- ggplot(dat_tt, aes(x = beta.exposure, y = beta.outcome)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_abline(intercept = 0, slope = ivw_b, color = COL_PRIMARY, linewidth = 0.7) +
  geom_abline(intercept = 0, slope = wm_b, color = COL_SECONDARY, linewidth = 0.5, linetype = "dotted") +
  geom_errorbar(aes(ymin = beta.outcome - se.outcome, ymax = beta.outcome + se.outcome),
    width = 0, color = "grey75", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = beta.exposure - se.exposure, xmax = beta.exposure + se.exposure),
    height = 0, color = "grey75", linewidth = 0.3) +
  geom_point(aes(size = F_stat, color = pleio, shape = pleio)) +
  scale_color_manual(values = c("FALSE" = COL_PRIMARY, "TRUE" = COL_RISK), guide = "none") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1), guide = "none") +
  scale_size_continuous(range = c(1.5, 5), name = "F-statistic",
    breaks = c(50, 700, 1500), guide = guide_legend(
      override.aes = list(color = COL_PRIMARY))) +
  ### FIX 5: nudge TFRC off the IVW regression line
  geom_text_repel(aes(label = gene), size = 2.3, family = FONT, fontface = "italic",
    color = ifelse(dat_tt$pleio, COL_RISK, "grey30"),
    nudge_y = ifelse(dat_tt$gene == "TFRC", 0.004, 0),
    nudge_x = ifelse(dat_tt$gene == "TFRC", -0.02, 0),
    max.overlaps = 15, segment.size = 0.2, seed = 42, box.padding = 0.6,
    min.segment.length = 0.3, force = 2) +
  annotate("segment", x = 0.20, xend = 0.28, y = 0.028, yend = 0.028,
    color = COL_PRIMARY, linewidth = 0.7) +
  annotate("text", x = 0.29, y = 0.028, label = "IVW", hjust = 0, size = 2.2,
    family = FONT, color = COL_PRIMARY) +
  annotate("segment", x = 0.20, xend = 0.28, y = 0.024, yend = 0.024,
    color = COL_SECONDARY, linewidth = 0.5, linetype = "dotted") +
  annotate("text", x = 0.29, y = 0.024, label = "Weighted median", hjust = 0, size = 2.2,
    family = FONT, color = COL_SECONDARY) +
  labs(x = "SNP effect on transferrin (SD units)",
       y = "SNP effect on circulating total-tau (SD units)",
       caption = PLEIO_FOOTNOTE) +
  theme_minimal(base_size = 9, base_family = FONT) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(t = 2),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7.5),
    legend.key.size = unit(10, "pt"),
    plot.caption = element_text(size = 6, color = "grey45", hjust = 0, margin = margin(t = 4))
  )

### FIX 1-3: statistics box, drawn line by line so the exponent can be plotmath
yr <- ggplot_build(fig2_base)$layout$panel_params[[1]]$y.range
dy <- 0.042 * diff(yr)          # line spacing in data units
x0 <- -0.60                     # box left edge (old label sat at x = -0.55)
y0 <- yr[1] + 0.055 * diff(yr)  # box bottom
lines <- list(
  list(txt = "n = 8 SNPs", parse = FALSE),
  list(txt = paste0('paste("IVW: ", beta, " = ', MINUS, '0.039, p = ", "2.88×10"^"', MINUS, '3")'), parse = TRUE),
  list(txt = "Egger intercept p = 0.202", parse = FALSE),
  list(txt = 'paste(I^2, " = 0%, Q p = 0.544")', parse = TRUE),
  list(txt = "MR-PRESSO p = 0.570", parse = FALSE)
)
fig2 <- fig2_base +
  annotate("rect", xmin = x0, xmax = x0 + 0.315,
           ymin = y0 - 0.55 * dy, ymax = y0 + (length(lines) - 0.35) * dy,
           fill = "white", color = "grey25", linewidth = 0.3) +
  lapply(seq_along(lines), function(i) {
    l <- lines[[length(lines) - i + 1]]
    annotate("text", x = x0 + 0.012, y = y0 + (i - 1) * dy, label = l$txt,
             parse = l$parse, hjust = 0, vjust = 0.5, size = 2.2,
             family = FONT, color = "grey25")
  })

save_fig(file.path(OUT, "Fig2_scatter_Tf_tau"), fig2, mm2in(130), mm2in(140))

# ============================================================
# FigS1 — leave-one-out
# ============================================================

cat("FigS1: LOO...\n")

loo <- mr_leaveoneout(dat_tt)
loo$label <- ifelse(loo$SNP == "All", "All SNPs",
  paste0(loo$SNP, " (", gene_labels[loo$SNP], ")"))
loo$ci_lo <- loo$b - 1.96 * loo$se
loo$ci_hi <- loo$b + 1.96 * loo$se
loo$is_driver <- loo$SNP == "rs8177240"
loo$is_all <- loo$SNP == "All"
loo$row <- rev(seq_len(nrow(loo)))

pt_color <- ifelse(loo$is_all, "black", ifelse(loo$is_driver, COL_RISK, COL_PRIMARY))
pt_shape <- ifelse(loo$is_all, 18, 16)
pt_size <- ifelse(loo$is_all, 3.5, 2)
ci_color <- pt_color

figs1 <- ggplot(loo, aes(x = b, y = row)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
    ymin = loo$row[loo$is_driver] - 0.45, ymax = loo$row[loo$is_driver] + 0.45,
    fill = "#FFFDE7") +
  geom_hline(yintercept = loo$row[loo$is_all] + 0.5, color = "grey40",
    linewidth = 0.5, linetype = "solid") +
  geom_vline(xintercept = loo$b[loo$is_all], linetype = "dashed",
    color = COL_PRIMARY, linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0,
    color = ci_color, linewidth = 0.5) +
  geom_point(color = pt_color, size = pt_size, shape = pt_shape) +
  ### FIX 1: β and true minus in the callout
  annotate("label", x = 0.005, y = loo$row[loo$is_driver] - 0.75,
    label = paste0("Without rs8177240 (TF):\nβ = ", MINUS,
                   "0.020, p = 0.358\n(direction preserved, significance lost)"),
    size = 2, family = FONT, color = COL_RISK, hjust = 0.5,
    fill = "#FFFDE7", linewidth = 0.3, label.padding = unit(3, "pt"),
    fontface = "bold") +
  scale_y_continuous(breaks = loo$row, labels = loo$label,
    expand = expansion(add = c(0.8, 0.5))) +
  labs(x = expression("IVW " * beta * " for effect on circulating total-tau (leave-one-out)"),
       y = NULL, caption = PLEIO_FOOTNOTE) +
  theme_minimal(base_size = 9, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    ### FIX 4: breaks are supplied in loo$row order, so the face vector must be too
    axis.text.y = element_text(
      face = ifelse(loo$is_all, "bold", "plain"),
      size = 7.5, family = FONT),
    plot.caption = element_text(size = 6, color = "grey45", hjust = 0, margin = margin(t = 4)))

save_fig(file.path(OUT, "FigS1_LOO"), figs1, mm2in(160), mm2in(135))
cat("  bold row:", loo$label[loo$is_all], "\n")

# ============================================================
# FigS2 — funnel
# ============================================================

cat("FigS2: Funnel...\n")

ss <- mr_singlesnp(dat_tt)
ss_pts <- ss[!grepl("^All", ss$SNP), ]
ss_pts$gene <- gene_labels[ss_pts$SNP]
ivw_est <- ss$b[grepl("Inverse variance", ss$SNP)]

figs2 <- ggplot(ss_pts, aes(x = b, y = 1/se)) +
  geom_vline(xintercept = ivw_est, linetype = "dashed", color = COL_PRIMARY, linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(color = COL_PRIMARY, size = 2.5) +
  ### FIX 5: push HLA*/NAT2 clear of their points and of the dashed line
  geom_text_repel(aes(label = gene), size = 2.3, family = FONT, fontface = "italic",
    color = "grey40", seed = 42, segment.size = 0.2, max.overlaps = 20,
    box.padding = 0.8, force = 4, min.segment.length = 0,
    nudge_x = ifelse(ss_pts$gene %in% c("HLA*", "NAT2"), -0.012, 0),
    nudge_y = ifelse(ss_pts$gene == "NAT2", 4, ifelse(ss_pts$gene == "HLA*", -3, 0))) +
  ### FIX 2: I²
  annotate("label", x = 0.20, y = 55,
    label = "Egger intercept: p = 0.202\nCochran Q: p = 0.544, I² = 0%",
    size = 2.2, family = FONT, fill = "#F5F5F5", linewidth = 0.3, hjust = 1) +
  labs(x = "Per-SNP MR estimate (Wald ratio)", y = expression("Precision (1/SE)"),
       caption = PLEIO_FOOTNOTE) +
  theme_minimal(base_size = 9, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    plot.caption = element_text(size = 6, color = "grey45", hjust = 0, margin = margin(t = 4)))

save_fig(file.path(OUT, "FigS2_funnel"), figs2, mm2in(130), mm2in(120))
