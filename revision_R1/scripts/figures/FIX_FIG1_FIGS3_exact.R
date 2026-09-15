##############################################################################
# FIX_FIG1_FIGS3.R  (2026-09-13)
# Re-renders ONLY Figure 1 and Supplementary Figure S3.
#
# Both blocks come from scripts/APPLY_REVIEWER_FIXES.R (Fig 1: lines 46-155;
# FigS3: lines 415-459), plus the two changes made later by the now-lost
# APPLY_MCH.R (control relabelled to MCH; hippocampal-volume row dropped from
# Fig 1; FigS3 given a bold title).
#
# Do NOT run APPLY_REVIEWER_FIXES.R to regenerate these: it also rewrites
# Figure 3 from a superseded code path (retired serum-iron p = 0.817).
#
# Fixes applied here (audit, 2026-09-13):
#   1. Depression and multiple sclerosis were plotted with standard errors ~30x
#      too large (0.029 / 0.028 instead of 0.00091 / 0.00039), so their CIs were
#      drawn at +/-0.055 instead of hairline width. Both figures were affected.
#   2. The two control rows did not match Table S2 (MCH 0.70 se 0.18; height
#      0.03 se 0.02) -> now +0.704 se 0.1929 and +0.031 se 0.02296.
#   3. "CSF p-tau - abnormal amyloid" was missing from Fig 1 although every other
#      CSF p-tau subgroup was plotted -> added (Exploratory).
#   4. Bold y-axis labels were mirrored: the code indexed the significance vector
#      with d$sig[order(d$row)], but vectorised element_text() is applied in the
#      order the breaks are supplied (d$row, descending). Bold therefore landed on
#      the reverse rows. Now plain d$sig, i.e. bold = nominally significant.
#   5. The clipped-CI arrow for the positive control stopped short of its own
#      point; the CI is now clipped at 0.72 and the arrow runs to the axis edge.
#   6. Notation: ASCII "->" -> "→", "-ve" -> "−ve", "Abeta42" -> "Aβ42", and true
#      minus signs on the x-axis tick labels.
#   All plotted values are taken from Supplementary Table S2 (CI -> se = (hi-lo)/3.92).
#
# Usage: Rscript -e ".libPaths(...); source('FIX_FIG1_FIGS3.R')" <out_dir>
##############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(grid)
})
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[1] else
  "C:/Users/ashwi/OneDrive/Documents/EB1A Docs/MR Paper/Final Figures and data/figures/_staging_F1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

FONT <- "Arial"
COL_PRIMARY <- "#1B4F8A"
COL_SECONDARY <- "#4A90D9"
COL_NULL <- "#AAAAAA"
COL_CONTROL <- "#2E7D32"
COL_SPEC <- "#6A1B9A"
MINUS <- "\u2212"   # true minus; Arial has no superscript minus, this one is fine
mm2in <- function(mm) mm / 25.4

save_fig <- function(path, plot, w, h, dpi = 300) {
  ggsave(paste0(path, ".png"), plot, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(path, ".pdf"), plot, width = w, height = h, device = cairo_pdf, bg = "white")
  cat("  Saved:", basename(path), "\n")
}

# ============================================================
# FIGURE 1 — gradient forest
# ============================================================

cat("Fig 1: Gradient forest...\n")

### FIX_FIG1_FIGS3_exact.R (2026-09-14): identical to FIX_FIG1_FIGS3.R except
###   (a) full-precision beta/SE/p from the saved analysis outputs replace the values
###       back-calculated from Table S2's rounded CIs (se = (hi-lo)/3.92), which had
###       left e.g. the primary row at SE 0.01276 instead of 0.01318. Sources:
###       data/FINAL_analysis_results.RData (primary_results, controls),
###       results/reviewer_response/E_MVMR_extended.csv (MVMR),
###       results/primary/CSF_all_Jansen_subanalyses_2026-03-29.csv (CSF rows),
###       results/primary/triangulation_AD_diagnosis_2026-03-29.csv (AD),
###       results/primary/causal_gradient_2026-03-29.csv (NfL, depression).
###       FTD-TDP, ALS (excl. HFE), PD, stroke and MS have no saved full-precision
###       estimate; their Table S2 values already reproduce their p and are kept.
###   (b) "Entorhinal thickness (bilateral)" was a plain average of the L and R
###       results (beta, SE and p averaged; p = 0.379 did not match its own beta/SE)
###       -> replaced by the separate L and R estimates, as in Table S2.
###   (c) the MVMR row uses a t-based 95% CI, qt(0.975, 8), matching its t-based p
###       (as in Fig 3 / Fig S7); every other row keeps 1.96 (normal-based p).
###   (d) FigS3: the same full-precision control values.
d <- tibble::tribble(
  ~cat,          ~outcome,                              ~beta,                 ~se,                  ~pval,               ~crit,
  "PRIMARY",     "Circulating total-tau",               -0.0392719,            0.0131775,            0.00288032,          1.96,
  "PRIMARY",     "MVMR (indep. of serum iron)",         -0.0436538820642196,   0.00915310475559959,  0.00140988031217651, qt(0.975, 8),
  "SECONDARY",   "CSF p-tau - APOE4 non-carriers",      -0.116409635519206,    0.0587440239457868,   0.04751930306598,    1.96,
  "SECONDARY",   "FTD-TDP",                             -0.306,                0.13980,              0.028,               1.96,
  "EXPLORATORY", "ALS (excl. HFE)",                     -0.060,                0.02602,              0.020,               1.96,
  "EXPLORATORY", "CSF p-tau - main analysis",           -0.0497812042592195,   0.0331934875230201,   0.133684861967294,   1.96,
  "EXPLORATORY", "NfL heavy chain",                     -0.0799336546721662,   0.0544261027832454,   0.141923962163707,   1.96,
  "EXPLORATORY", "CSF p-tau - normal amyloid",          -0.0752718117981438,   0.0568209819611203,   0.185264500948273,   1.96,
  "EXPLORATORY", "Entorhinal thickness (L)",            -0.0184258,            0.0151792,            0.22479181,          1.96,
  "EXPLORATORY", "Entorhinal thickness (R)",            -0.0203044,            0.0325758,            0.53308860,          1.96,
  "EXPLORATORY", "CSF p-tau - abnormal amyloid",        -0.0180682294968807,   0.0490627726982219,   0.712673709774801,   1.96,
  "NULL",        "PD diagnosis",                         0.050,                0.03520,              0.156,               1.96,
  "NULL",        "AD diagnosis",                         0.0181094472946676,   0.0216122755699681,   0.402073313315196,   1.96,
  "NULL",        "Stroke",                              -0.0002,               0.00051,              0.663,               1.96,
  "NULL",        "CSF Aβ42 - main",                 0.0121432139664284,   0.0337084000253439,   0.71866544722239,    1.96,
  "NULL",        "CSF p-tau - APOE4 carriers",           0.0155438648999537,   0.05292814961546,     0.769003508023431,   1.96,
  "NULL",        "Multiple sclerosis",                   0.0001,               0.00041,              0.779,               1.96,
  "NULL",        "Depression",                          -0.000128644224388586, 0.000910682232183447, 0.887663461794211,   1.96,
  "NULL",        "CSF Aβ42 - APOE4 carriers",      -0.00209200962265214,  0.0512259575875278,   0.967424360750609,   1.96,
  "NULL",        "CSF Aβ42 - APOE4 non-carriers",   0.00151099125038972,  0.0516257969148371,   0.976650733538908,   1.96,
  "CONTROL",     paste0("Iron → MCH (+ve control)"), 0.704406,            0.192681,             0.000256371,         1.96,
  "CONTROL",     paste0("Iron → height (", MINUS, "ve control)"), 0.0305968, 0.0233856,          0.19075,             1.96
)
d <- as.data.frame(d)
d$ci_lo <- d$beta - d$crit * d$se
d$ci_hi <- d$beta + d$crit * d$se
d$sig <- d$pval < 0.05
d$cat <- factor(d$cat, levels = c("PRIMARY", "SECONDARY", "EXPLORATORY", "NULL", "CONTROL"))
d <- d %>% arrange(cat, pval) %>% mutate(row = rev(seq_len(n())))

d$cat_label <- case_when(
  d$cat == "PRIMARY" ~ "Primary", d$cat == "SECONDARY" ~ "Secondary",
  d$cat == "EXPLORATORY" ~ "Exploratory", d$cat == "NULL" ~ "Null",
  d$cat == "CONTROL" ~ "Control")
CAT_ORDER <- c("Primary", "Secondary", "Exploratory", "Null", "Control")
d$cat_label <- factor(d$cat_label, levels = CAT_ORDER)

### FIX 5: clip the CI at 0.72 (right of the MCH point at 0.704) and run the
### arrow from there to the axis edge, so it no longer stops short of its point.
CLIP_R <- 0.72
CLIP_L <- -0.65
d$ci_lo_clip <- pmax(d$ci_lo, CLIP_L)
d$ci_hi_clip <- pmin(d$ci_hi, CLIP_R)
d$clipped_right <- d$ci_hi > CLIP_R
d$clipped_left <- d$ci_lo < CLIP_L

fig1 <- ggplot(d, aes(x = beta, y = row)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
    ymin = d$row[c(TRUE, FALSE)] - 0.45, ymax = d$row[c(TRUE, FALSE)] + 0.45,
    fill = "#F8F8F8") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_segment(aes(x = ci_lo_clip, xend = ci_hi_clip, yend = row),
    color = case_when(d$cat == "PRIMARY" ~ COL_PRIMARY, d$cat == "SECONDARY" ~ COL_SECONDARY,
      d$cat == "EXPLORATORY" ~ COL_SECONDARY, d$cat == "NULL" ~ COL_NULL,
      d$cat == "CONTROL" ~ COL_CONTROL),
    linewidth = ifelse(d$sig, 0.9, 0.5)) +
  geom_segment(data = d[d$clipped_right, ], aes(x = CLIP_R, xend = 0.75, yend = row),
    arrow = arrow(length = unit(1.5, "mm"), type = "open"),
    color = COL_CONTROL, linewidth = 0.5) +
  geom_segment(data = d[d$clipped_left, ], aes(x = CLIP_L, xend = -0.68, yend = row),
    arrow = arrow(length = unit(1.5, "mm"), type = "open"),
    color = COL_SECONDARY, linewidth = 0.5) +
  geom_point(aes(shape = cat_label, color = cat_label),
    size = ifelse(d$cat == "PRIMARY", 3, 2.5)) +
  scale_shape_manual(name = NULL, breaks = CAT_ORDER,
    values = c("Primary" = 15, "Secondary" = 16, "Exploratory" = 1,
               "Null" = 18, "Control" = 17)) +
  scale_color_manual(name = NULL, breaks = CAT_ORDER,
    values = c("Primary" = COL_PRIMARY, "Secondary" = COL_SECONDARY,
               "Exploratory" = COL_SECONDARY, "Null" = COL_NULL,
               "Control" = COL_CONTROL)) +
  scale_y_continuous(breaks = d$row, labels = d$outcome,
    expand = expansion(add = c(0.5, 1))) +
  scale_x_continuous(limits = c(-0.65, 0.78),
    breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6),
    labels = c(paste0(MINUS, "0.6"), paste0(MINUS, "0.4"), paste0(MINUS, "0.2"),
               "0", "0.2", "0.4", "0.6")) +
  labs(x = expression("Effect estimate (IVW " * beta * ", 95% CI)"), y = NULL,
    caption = paste0("Controls use serum iron (ieu-a-1049, 3 SNPs) as exposure, not transferrin. ",
                     "Arrows indicate CI extends beyond axis range.\n",
                     "Bold labels and thicker lines mark nominal p < 0.05.")) +
  theme_minimal(base_size = 8, base_family = FONT) +
  theme(
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.3),
    ### FIX 4: breaks are supplied in d$row order, so the face vector must be too
    axis.text.y = element_text(size = 6.5, hjust = 1,
      face = ifelse(d$sig, "bold", "plain")),
    legend.position = c(0.87, 0.42),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    legend.margin = margin(3, 6, 3, 6),
    legend.key.size = unit(10, "pt"),
    legend.text = element_text(size = 7, family = FONT),
    plot.caption = element_text(size = 6, color = "grey45", hjust = 0, margin = margin(t = 6)),
    plot.margin = margin(5, 5, 5, 5)
  )

save_fig(file.path(OUT, "Fig1_gradient_forest"), fig1, mm2in(170), mm2in(210))
cat("  rows top-to-bottom:\n")
print(d[order(-d$row), c("cat", "outcome", "beta", "se", "pval", "sig")], row.names = FALSE)

# ============================================================
# FigS3 — pipeline / specificity controls
# ============================================================

cat("FigS3: Controls...\n")

ctrl <- data.frame(
  label = c(paste0("Iron \u2192 MCH\n(positive control)"),
            paste0("Iron \u2192 height\n(negative control)"),
            paste0("Transferrin \u2192 depression\n(neurodegeneration-specificity control)")),
  beta = c(0.704406, 0.0305968, -0.000128644224388586),
  se   = c(0.192681, 0.0233856, 0.000910682232183447),
  color = c(COL_CONTROL, COL_CONTROL, COL_SPEC),
  shape = c(17, 2, 5),
  note = c("Significant (+ve)", "Null", "Null (specificity)"),
  stringsAsFactors = FALSE
)
ctrl$ci_lo <- ctrl$beta - 1.96 * ctrl$se
ctrl$ci_hi <- ctrl$beta + 1.96 * ctrl$se
ctrl$row <- rev(seq_len(nrow(ctrl)))

figs3 <- ggplot(ctrl, aes(x = beta, y = row)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15,
    color = ctrl$color, linewidth = 0.7) +
  geom_point(shape = ctrl$shape, size = 3.5, color = ctrl$color) +
  geom_text(aes(label = note), nudge_y = -0.25, size = 2.5, family = FONT,
    color = "grey40", fontface = "italic") +
  scale_y_continuous(breaks = ctrl$row, labels = ctrl$label,
    expand = expansion(add = c(0.8, 0.5))) +
  scale_x_continuous(limits = c(-0.15, 1.15), breaks = seq(0, 1.0, 0.25)) +
  labs(title = "Pipeline and specificity controls",
    x = expression("Effect estimate (IVW " * beta * ", 95% CI)"), y = NULL) +
  theme_minimal(base_size = 9, base_family = FONT) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 9, margin = margin(b = 6)))

save_fig(file.path(OUT, "FigS3_controls"), figs3, mm2in(170), mm2in(85))
cat("  S3 CIs:\n")
print(ctrl[, c("beta", "se", "ci_lo", "ci_hi")], row.names = FALSE)
