##############################################################################
# R1_FIGURES.R — New figures for reviewer response
#
# Matches the style established in PRODUCE_ALL_FIGURES_v2.R:
#   - Colors: COL_PRIMARY #1B4F8A, COL_SECONDARY #4A90D9, COL_NULL #AAAAAA,
#     COL_CONTROL #2E7D32, COL_RISK #C62828
#   - Font Arial, base 9pt
#   - save_fig() writes PNG (300 dpi) and cairo_pdf
#   - Sizing via mm2in(); typical 130–175 mm wide
#   - Zero-line dashed #888888, panel border black 0.5
#
# Produces:
#   FigS5_coloc_TF_locus       — regional plot of colocalization
#   FigS6_leave_locus_out      — forest plot of subsets
#   FigS7_MVMR_extended        — 2-way + 3-way MVMR with cond-F
#   FigS8_causal_DAG           — DAG (reviewer explicit ask)
#   FigS9_cis_TF_comparison    — full 8-SNP vs cis-only comparison
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(patchwork); library(dplyr)
  library(grid); library(ggforce); library(ggdag); library(dagitty)
})
set.seed(42)

FONT <- "Arial"
BASE <- file.path(PROJECT_DIR, "Final Figures and data")
SUPP <- file.path(BASE, "figures", "supplementary")
REV  <- file.path(PROJECT_DIR, "results", "reviewer_response")
dir.create(SUPP, showWarnings = FALSE, recursive = TRUE)

COL_PRIMARY   <- "#1B4F8A"
COL_SECONDARY <- "#4A90D9"
COL_NULL      <- "#AAAAAA"
COL_CONTROL   <- "#2E7D32"
COL_RISK      <- "#C62828"

mm2in <- function(mm) mm / 25.4

save_fig <- function(path, plot, w, h, dpi = 300) {
  ggsave(paste0(path, ".png"), plot, width = w, height = h, dpi = dpi, bg = "white")
  ggsave(paste0(path, ".pdf"), plot, width = w, height = h,
         device = cairo_pdf, bg = "white")
  cat("  Saved:", basename(path), "\n")
}

theme_forest <- function(base = 9) {
  theme_minimal(base_size = base, base_family = FONT) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "#EBEBEB", linewidth = 0.2),
      axis.ticks         = element_line(color = "black", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = base,
                                  margin = margin(b = 6))
    )
}

cat("=== FIGURES for reviewer response ===\n\n")

# ============================================================
# FigS5 — Colocalization regional plot at TF locus
# Two vertical tracks sharing x-axis:
#   Top:    transferrin -log10(p)
#   Bottom: total-tau  -log10(p)
# Highlight rs8177240; annotate PP.H4 = 0.447, PP.H3 = 0.017
# ============================================================

cat("FigS5: colocalization regional plot\n")
coloc <- read.csv(file.path(REV, "D_coloc_TF_snp_posteriors.csv"),
                   stringsAsFactors = FALSE)
coloc$p_tf  <- 2 * pnorm(-abs(coloc$z.df1))
coloc$p_tau <- 2 * pnorm(-abs(coloc$z.df2))
coloc$logp_tf  <- -log10(pmax(coloc$p_tf,  1e-300))
coloc$logp_tau <- -log10(pmax(coloc$p_tau, 1e-300))
coloc$mb <- coloc$position / 1e6

# TF gene body coordinates
TF_LO <- 133.464073
TF_HI <- 133.494388

# Annotate top variant
top_h4 <- coloc[order(-coloc$SNP.PP.H4)[1], ]
cat(sprintf("  Top PP.H4 SNP: %s at %.4f Mb (PP.H4=%.3f)\n",
            top_h4$snp, top_h4$mb, top_h4$SNP.PP.H4))

# All five coloc.abf posteriors (D_coloc_TF_summary.csv), 3 dp, sum = 1.000
pp_h0 <- 0.000; pp_h1 <- 0.536; pp_h2 <- 0.000; pp_h3 <- 0.017; pp_h4 <- 0.447

# Gradient color by SNP.PP.H4
coloc$is_top <- coloc$snp == top_h4$snp

# Extra top headroom so the italic "TF" gene label and the rs8177240 callout
# sit clearly ABOVE the tallest points instead of overlapping them.
y_lim_tf  <- c(0, max(coloc$logp_tf,  na.rm = TRUE) * 1.20)
y_lim_tau <- c(0, max(4.3, max(coloc$logp_tau, na.rm = TRUE) * 1.20))

make_track <- function(logp, top_row, title, y_lab, y_lim, annotation = NULL) {
  ggplot(coloc, aes(x = mb, y = logp)) +
    # TF gene body shading
    annotate("rect", xmin = TF_LO, xmax = TF_HI,
             ymin = -Inf, ymax = Inf,
             fill = "#E3F2FD", alpha = 0.6) +
    # Gene label parked in the headroom, above the peak points
    annotate("text", x = (TF_LO + TF_HI) / 2, y = y_lim[2] * 0.975,
             label = "TF", family = FONT, fontface = "italic",
             size = 2.5, color = COL_PRIMARY) +
    geom_point(color = ifelse(coloc$is_top, COL_RISK, COL_PRIMARY),
               alpha = ifelse(coloc$is_top, 1, 0.55),
               size  = ifelse(coloc$is_top, 2.5, 1)) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
               color = "#888888", linewidth = 0.3) +
    annotate("text", x = min(coloc$mb) + 0.05,
             y = -log10(5e-8) + y_lim[2] * 0.03,
             label = "p = 5×10⁻⁸", family = FONT, size = 2, color = "grey40") +
    # rs8177240 callout nudged up-and-left so it clears both the peak and the
    # "TF" gene label (which sits at band-centre in the headroom)
    geom_text_repel(data = top_row, aes(label = snp),
                    color = COL_RISK, size = 2.5, family = FONT,
                    fontface = "bold", segment.color = COL_RISK,
                    segment.size = 0.3,
                    nudge_y = y_lim[2] * 0.05, nudge_x = -0.055,
                    box.padding = 0.5, seed = 42) +
    { if (!is.null(annotation))
        annotate("label", x = max(coloc$mb) - 0.02, y = y_lim[2] * 0.95,
                 label = annotation, family = FONT, size = 2.3,
                 hjust = 1, vjust = 1, fill = "#F5F5F5",
                 label.size = 0.3, label.padding = unit(3, "pt"))
    } +
    scale_y_continuous(limits = y_lim, expand = expansion(mult = c(0, 0.02))) +
    labs(title = title, x = NULL, y = y_lab) +
    theme_forest(base = 8) +
    theme(axis.text.x = element_text(size = 7))
}

# Single joint (two-trait) quantity — describes both panels, so it lives once in A.
# All five posteriors shown (H0/H2 previously omitted).
annot_h4 <- sprintf(paste0("coloc.abf (transferrin vs. tau)\n",
                    "PP.H0 = %.3f   PP.H2 = %.3f\n",
                    "PP.H1 = %.3f  (Tf only)\n",
                    "PP.H3 = %.3f  (distinct variants)\n",
                    "PP.H4 = %.3f  (shared variant)\n",
                    "Conditional H4/(H3+H4) = %.1f%%"),
                    pp_h0, pp_h2, pp_h1, pp_h3, pp_h4, 100 * pp_h4 / (pp_h3 + pp_h4))

p_tf  <- make_track(coloc$logp_tf,  top_h4[, c("snp","mb","logp_tf")]  %>%
                       rename(logp = logp_tf),
                     "A  Transferrin (Benyamin 2014, n = 23,986)",
                     expression(-log[10] * "(p)"),
                     y_lim_tf, annot_h4)
p_tau <- make_track(coloc$logp_tau, top_h4[, c("snp","mb","logp_tau")] %>%
                       rename(logp = logp_tau),
                     "B  Circulating total-tau (Sarnowski 2022, n = 14,721)",
                     expression(-log[10] * "(p)"),
                     y_lim_tau) +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005))) +
  labs(x = "Chromosome 3 position (Mb, GRCh37)")

figs5 <- p_tf / p_tau + plot_layout(heights = c(1, 1)) +
  plot_annotation(caption = paste0(
    "PP.H4 = 0.447 is below the conventional 0.8 threshold; the conditional H4/(H3+H4) = 96.3% is reported because the\n",
    "residual mass sits on PP.H1 (transferrin-only), reflecting limited power in the tau GWAS at this locus (min p ≈ 10⁻³),\n",
    "not evidence against colocalization (PP.H3, distinct causal variants, ≈ 0.02)."),
    theme = theme(plot.caption = element_text(family = FONT, size = 6.5,
                  hjust = 0, color = "grey30", margin = margin(t = 4))))
save_fig(file.path(SUPP, "FigS5_coloc_TF_locus"),
         figs5, mm2in(175), mm2in(140))

# ============================================================
# FigS6 — Leave-locus-out forest
# ============================================================

cat("FigS6: leave-locus-out forest\n")
llo <- read.csv(file.path(REV, "B_leave_locus_out.csv"),
                 stringsAsFactors = FALSE)

# Human-readable labels & display order (top = primary, bottom = single-SNP)
label_map <- c(
  "all_8"                                                       = "All 8 SNPs (primary)",
  "excl_rs8177240"                                              = "Excl. rs8177240 (single-SNP LOO)",
  "excl_HLA"                                                    = "Excl. HLA (rs9268633)",
  "excl_FADS"                                                   = "Excl. FADS (rs174577)",
  "excl_HFE"                                                    = "Excl. HFE (rs1800562)",
  "excl_NAT2"                                                   = "Excl. NAT2 (rs1495741)",
  "excl_HLA_FADS"                                               = "Excl. HLA + FADS",
  "excl_HLA_FADS_HFE_NAT2 (joint)"                              = "Excl. HLA + FADS + HFE + NAT2 (joint)",
  "cis_TF_3snps (rs8177240,rs9990333,rs7646473)"                = "TF/TFRC-region subset\nof primary panel (3 SNPs)",
  "TF_only (rs8177240)"                                         = "rs8177240 alone (Wald)"
)
llo$display <- label_map[llo$subset]
llo <- llo[!is.na(llo$display), ]
llo$display <- factor(llo$display, levels = rev(unname(label_map)))

# Standardize on multiplicative random-effects IVW (matches TwoSampleMR, the
# convention behind the manuscript primary). Inflate the fixed-effect SE by the
# overdispersion factor sqrt(Q/(n-1)), floored at 1 (never shrink below FE).
# Reconciliation 1.2: resolves the Fig S6 (FE) vs manuscript/Fig S1 (RE) mismatch.
# Rows with I2 below the flooring threshold are unchanged; only the over-dispersed
# subsets (excl_HLA_FADS, joint 4-SNP, TF-region 3-SNP) widen.
llo$phi     <- ifelse(llo$n_snps > 1 & !is.na(llo$ivw_Q),
                      sqrt(pmax(1, llo$ivw_Q / (llo$n_snps - 1))), 1)
llo$re_se   <- llo$ivw_se * llo$phi
llo$re_pval <- 2 * pnorm(-abs(llo$ivw_beta / llo$re_se))
llo$ci_lo <- llo$ivw_beta - 1.96 * llo$re_se
llo$ci_hi <- llo$ivw_beta + 1.96 * llo$re_se

# Colour scheme
llo$color <- COL_PRIMARY
llo$color[grepl("joint", llo$subset, fixed = TRUE) |
           grepl("cis_TF_3snps", llo$subset, fixed = TRUE)] <- COL_CONTROL
llo$color[llo$subset == "excl_rs8177240"] <- COL_RISK
llo$color[llo$subset == "TF_only (rs8177240)"] <- COL_SECONDARY
llo$sig <- llo$re_pval < 0.05
llo$plab <- ifelse(is.na(llo$re_pval), "n/a",
             ifelse(llo$re_pval < 1e-3,
                    formatC(llo$re_pval, format = "e", digits = 2),
                    sprintf("%.3f", llo$re_pval)))

# For the single-SNP row, ivw_se is available (Wald SE)
figs6 <- ggplot(llo, aes(x = ivw_beta, y = display)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_vline(xintercept = llo$ivw_beta[llo$subset == "all_8"],
             linetype = "dotted", color = COL_PRIMARY, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                 color = llo$color, linewidth = 0.6) +
  geom_point(color = llo$color, size = 3,
             shape = ifelse(llo$sig, 15, 22),
             fill = "white") +
  # p-value column parked to the right of every CI (widest is the red
  # excl-rs8177240 row, whose CI reaches ~+0.023) so nothing overlaps
  geom_text(aes(x = 0.045, label = paste0("p = ", plab)),
            hjust = 0, size = 2.4, family = FONT,
            fontface = ifelse(llo$sig, "bold", "plain"),
            color = ifelse(llo$sig, llo$color, "grey40")) +
  scale_x_continuous(limits = c(-0.11, 0.105),
                     breaks = seq(-0.10, 0.02, 0.02),
                     expand = expansion(add = 0)) +
  labs(title = "Leave-locus-out sensitivity (transferrin → circulating total-tau)",
       x = expression("IVW " * beta * " (95% CI)"), y = NULL) +
  theme_forest() +
  theme(axis.text.y = element_text(size = 7.5,
    face = ifelse(rev(unname(label_map)) == "All 8 SNPs (primary)", "bold", "plain")))

save_fig(file.path(SUPP, "FigS6_leave_locus_out"),
         figs6, mm2in(175), mm2in(115))

# ============================================================
# FigS7 — Extended MVMR forest (2-way + 3-way + cond-F)
# ============================================================

cat("FigS7: extended MVMR forest\n")
mv2 <- read.csv(file.path(REV, "E_MVMR_extended.csv"),
                 stringsAsFactors = FALSE)
mv4 <- read.csv(file.path(REV, "F_MVMR_4way.csv"),
                 stringsAsFactors = FALSE)

# Univariable estimate from cached RData for reference row
if (file.exists(file.path(PROJECT_DIR, "data", "FINAL_analysis_results.RData"))) {
  e <- new.env()
  load(file.path(PROJECT_DIR, "data", "FINAL_analysis_results.RData"),
       envir = e)
  pr <- e$primary_results
  univar_row <- pr[pr$method == "Inverse variance weighted" &
                    pr$pair == "transferrin_circ_total_tau", ]
  univar_b   <- univar_row$b
  univar_se  <- univar_row$se
  univar_p   <- univar_row$pval
} else {
  univar_b <- -0.039; univar_se <- 0.013; univar_p <- 2.88e-3
}

# 3-way MVMR — MVMR package labels exposures positionally
mv3 <- mv4[mv4$model == "3-way (transferrin + ferritin + tsat)", ]
mv3$exposure_name <- c("transferrin", "ferritin", "tsat")

fdat <- rbind(
  data.frame(model = "Univariable IVW",
             exposure = "Transferrin",
             beta = univar_b, se = univar_se, pval = univar_p, cond_F = NA),
  data.frame(model = "2-way MVMR (transferrin + serum iron)",
             exposure = "Transferrin",
             beta = mv2$beta[mv2$exposure == "transferrin"],
             se   = mv2$se[mv2$exposure == "transferrin"],
             pval = mv2$pval[mv2$exposure == "transferrin"],
             cond_F = mv2$cond_F[mv2$exposure == "transferrin"]),
  data.frame(model = "2-way MVMR (transferrin + serum iron)",
             exposure = "Serum iron",
             beta = mv2$beta[mv2$exposure == "serum_iron"],
             se   = mv2$se[mv2$exposure == "serum_iron"],
             pval = mv2$pval[mv2$exposure == "serum_iron"],
             cond_F = mv2$cond_F[mv2$exposure == "serum_iron"]),
  data.frame(model = "3-way MVMR (transferrin + ferritin + TSAT)",
             exposure = "Transferrin",
             beta = mv3$beta[mv3$exposure_name == "transferrin"],
             se   = mv3$se[mv3$exposure_name == "transferrin"],
             pval = mv3$pval[mv3$exposure_name == "transferrin"],
             cond_F = mv3$cond_F[mv3$exposure_name == "transferrin"]),
  data.frame(model = "3-way MVMR (transferrin + ferritin + TSAT)",
             exposure = "Ferritin",
             beta = mv3$beta[mv3$exposure_name == "ferritin"],
             se   = mv3$se[mv3$exposure_name == "ferritin"],
             pval = mv3$pval[mv3$exposure_name == "ferritin"],
             cond_F = mv3$cond_F[mv3$exposure_name == "ferritin"]),
  data.frame(model = "3-way MVMR (transferrin + ferritin + TSAT)",
             exposure = "TSAT",
             beta = mv3$beta[mv3$exposure_name == "tsat"],
             se   = mv3$se[mv3$exposure_name == "tsat"],
             pval = mv3$pval[mv3$exposure_name == "tsat"],
             cond_F = mv3$cond_F[mv3$exposure_name == "tsat"])
)
fdat$row_lab <- paste(fdat$model, fdat$exposure, sep = " | ")
fdat$row_lab <- factor(fdat$row_lab, levels = rev(fdat$row_lab))
fdat$ci_lo <- fdat$beta - 1.96 * fdat$se
fdat$ci_hi <- fdat$beta + 1.96 * fdat$se
fdat$is_transferrin <- fdat$exposure == "Transferrin"
fdat$sig <- fdat$pval < 0.05
fdat$color <- ifelse(fdat$is_transferrin, COL_PRIMARY, COL_NULL)
# Uniform p-value labels across the panel: 2 significant figures throughout —
# scientific (m×10^e) below 0.01, plain decimal at/above. Keeps the 2-way MVMR
# transferrin estimate reading 1.4×10⁻³ (matches manuscript, Table S5, letter),
# not the ambiguous "0.001" that %.3f produced (Task 2).
.sup <- function(d) chartr("0123456789-",
  "⁰¹²³⁴⁵⁶⁷⁸⁹⁻", d)
fmt_p <- function(p) vapply(p, function(x){
  if (is.na(x)) return("")
  if (x < 0.01) { e <- floor(log10(x)); m <- x / 10^e
    paste0(formatC(m, format = "f", digits = 1), "×10", .sup(as.character(e))) }
  else if (x >= 0.1) sprintf("%.2f", x)
  else sprintf("%.3f", x)
}, character(1))
fdat$plab <- fmt_p(fdat$pval)
fdat$fann <- ifelse(is.na(fdat$cond_F), "",
                    sprintf("cond-F = %.1f", fdat$cond_F))

figs7 <- ggplot(fdat, aes(x = beta, y = row_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                 color = fdat$color, linewidth = 0.6) +
  geom_point(color = fdat$color, size = 3,
             shape = ifelse(fdat$is_transferrin, 15, 18)) +
  # Stat columns parked to the right of the widest CI (ferritin reaches ~+0.08)
  # so the p-values and cond-F never sit on top of an error bar
  geom_text(aes(x = 0.095, label = paste0("p = ", plab)),
            hjust = 0, size = 2.4, family = FONT,
            fontface = ifelse(fdat$sig, "bold", "plain"),
            color = ifelse(fdat$sig, fdat$color, "grey40")) +
  geom_text(aes(x = 0.15, label = fann),
            hjust = 0, size = 2.2, family = FONT, color = "grey40") +
  scale_x_continuous(limits = c(-0.13, 0.21),
                     breaks = seq(-0.10, 0.10, 0.05),
                     expand = expansion(add = 0)) +
  labs(title = "Extended MVMR reporting (transferrin → circulating total-tau)",
       x = expression("Effect estimate " * beta * " (95% CI)"), y = NULL) +
  theme_forest() +
  theme(axis.text.y = element_text(size = 7,
    face = ifelse(rev(fdat$is_transferrin), "bold", "plain")))

save_fig(file.path(SUPP, "FigS7_MVMR_extended"),
         figs7, mm2in(190), mm2in(105))

# ============================================================
# FigS8 — Causal DAG (reviewer explicit ask)
# ============================================================

cat("FigS8: causal DAG\n")

# Build DAG using dagitty syntax
# G_tf = genetic instruments for transferrin (cis TF variants)
# G_fe = genetic instruments for serum iron
# TF = transferrin protein
# Fe = serum iron
# NTBI = non-transferrin bound iron (mechanism)
# TAU = circulating total-tau
# APOE = APOE4 status (modifier)
# U = unmeasured confounders

# Layout tuned so labels sit inside their nodes and edges rarely cross:
#   instruments (left) -> exposures -> mechanism -> outcome (right);
#   confounder U parked top-left; downstream iron markers parked bottom-left.
dag <- dagify(
  TF   ~ G_tf + U,
  Fe   ~ G_fe + U,
  NTBI ~ TF + Fe,
  Tau  ~ NTBI + U,
  Ferritin ~ Fe,
  TSAT ~ Fe + TF,
  exposure = "TF",
  outcome  = "Tau",
  labels = c(
    G_tf = "TF-locus\nSNPs", G_fe = "Iron\nSNPs",
    TF = "Transferrin", Fe = "Serum\niron",
    NTBI = "NTBI",
    Tau = "Total\ntau",
    Ferritin = "Ferritin", TSAT = "TSAT",
    U = "Unmeasured\nconfounders"
  ),
  coords = list(
    x = c(G_tf = 0, G_fe = 0, TF = 2.0, Fe = 2.0, NTBI = 4.2,
          Tau = 6.0, Ferritin = 1.4, TSAT = 3.2, U = 1.2),
    y = c(G_tf = 4, G_fe = 1, TF = 4, Fe = 1, NTBI = 2.5,
          Tau = 2.5, Ferritin = -1.7, TSAT = -1.7, U = 5.8)
  )
)

dag_tidy <- tidy_dagitty(dag)

figs8 <- ggplot(dag_tidy, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_colour = "grey45",
                 arrow_directed = grid::arrow(length = grid::unit(2.2, "mm"),
                                              type = "closed")) +
  geom_dag_point(aes(fill = name), color = "black", shape = 21,
                 stroke = 0.5, size = 24, show.legend = FALSE) +
  geom_dag_text(aes(label = label), color = "grey10",
                family = FONT, fontface = "bold",
                size = 2.4, lineheight = 0.82) +
  scale_fill_manual(values = c(
    G_tf = "#F5F5F5", G_fe = "#F5F5F5",
    TF   = "#BBDEFB", Fe   = "#E0E0E0",
    NTBI = "#FFF9C4", Tau  = "#BBDEFB",
    Ferritin = "#EEEEEE", TSAT = "#EEEEEE",
    U = "#FFCDD2"
  )) +
  # column header for the instrument tier
  annotate("text", x = 0, y = 5.2, label = "Genetic\ninstruments",
           family = FONT, fontface = "italic", size = 2.3,
           color = "grey45", lineheight = 0.9) +
  # single interpretation box tucked into the empty lower-right region
  annotate("label", x = 4.25, y = 0.7,
    label = paste0(
      "Blue = exposure / outcome\n",
      "Pink = unmeasured confounders (U)\n",
      "Yellow = mechanism · grey = iron markers\n",
      "\n",
      "Serum iron is a parallel exposure that\n",
      "MVMR conditions on, not a mediator of\n",
      "transferrin's effect. NTBI = proposed\n",
      "mechanism. Ferritin & TSAT = downstream\n",
      "iron-status markers (added in 3-way MVMR)."),
    family = FONT, size = 2.2, hjust = 0, vjust = 1, lineheight = 1.15,
    fill = "#FFFDE7", label.size = 0.3, label.padding = unit(4, "pt")) +
  # Generous padding on all sides so the fixed-size circles (esp. the top U
  # node and the left-hand SNP nodes) are never cropped by the panel edge
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.13))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.16))) +
  coord_cartesian(clip = "off") +
  labs(title = "Causal DAG: transferrin → tau conditional on iron status") +
  theme_dag(base_family = FONT) +
  theme(plot.title = element_text(face = "bold", size = 9,
                                    margin = margin(b = 6)),
        plot.margin = margin(6, 6, 6, 6))

save_fig(file.path(SUPP, "FigS8_causal_DAG"),
         figs8, mm2in(200), mm2in(150))

# ============================================================
# FigS9 — Full 8-SNP panel vs cis-only 2-SNP comparison
# ============================================================

cat("FigS9: cis-only vs full-panel comparison\n")

cis <- read.csv(file.path(REV, "C_cis_TF_mr_500kb.csv"),
                 stringsAsFactors = FALSE)
cis_ivw <- cis[cis$method == "Inverse variance weighted", ]

comp <- data.frame(
  label = c(
    "Full 8-SNP panel (primary)\nrs8177240 + 7 iron-metabolism SNPs",
    "Cis-only TF locus (2 SNPs, ± 500 kb of TF)\nrs3811658 + rs17376530",
    "TF locus, single instrument (Wald)\nrs8177240 alone"),
  beta = c(univar_b, cis_ivw$b[1], -0.05031),
  # rs8177240 Wald SE = 0.01653 -> displays 0.017 (was 0.016); reconciliation 1.4
  se   = c(univar_se, cis_ivw$se[1], 0.01653),
  pval = c(univar_p, cis_ivw$pval[1], 2.343e-3),
  color = c(COL_PRIMARY, COL_CONTROL, COL_SECONDARY),
  stringsAsFactors = FALSE
)
comp$ci_lo <- comp$beta - 1.96 * comp$se
comp$ci_hi <- comp$beta + 1.96 * comp$se
comp$plab  <- formatC(comp$pval, format = "e", digits = 2)
comp$row   <- rev(seq_len(nrow(comp)))

BONF <- 0.05 / 23   # 0.00217

figs9 <- ggplot(comp, aes(x = beta, y = row)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                 color = comp$color, linewidth = 0.7) +
  geom_point(color = comp$color,
             shape = c(15, 17, 18), size = c(3.5, 4, 3.5)) +
  geom_text(aes(label = paste0("p = ", plab)),
            nudge_y = 0.28, size = 2.5, family = FONT,
            fontface = "bold",
            color = comp$color) +
  geom_text(aes(label = paste0("β = ", sprintf("%.3f", beta),
                                "  (SE ", sprintf("%.3f", se), ")")),
            nudge_y = -0.28, size = 2.2, family = FONT, color = "grey40") +
  scale_y_continuous(breaks = comp$row, labels = comp$label,
                     expand = expansion(add = c(0.7, 0.7))) +
  scale_x_continuous(limits = c(-0.11, 0.02)) +
  labs(title = "Comparison of instrument sets: transferrin → circulating total-tau",
       x = expression("IVW " * beta * " (95% CI)"), y = NULL,
       caption = paste0("Bonferroni threshold for 23 outcome pairs: p < ",
                        sprintf("%.4f", BONF), ".\n",
                        "The cis-only estimate clears Bonferroni; the full 8-SNP panel does not.")) +
  theme_forest() +
  theme(axis.text.y = element_text(size = 7.5, lineheight = 1.0),
        plot.caption = element_text(size = 7, hjust = 0,
                                      margin = margin(t = 6),
                                      color = "grey30"))

save_fig(file.path(SUPP, "FigS9_cis_TF_comparison"),
         figs9, mm2in(170), mm2in(100))

cat("\n=== DONE — 5 figures written to", SUPP, "===\n")
