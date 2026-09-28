##############################################################################
# FIX_FIGS5_S7_vector.R
# Re-renders ONLY Supplementary Figures S5 and S7 so their vector PDFs are correct.
#
# Both blocks are copied from R1_FIGURES.R. The only problem being fixed: the
# labels used the Unicode superscript minus (U+207B), which Arial does not
# contain, so the cairo PDFs printed "1.4x10[box]3". The PNGs were fine because
# the PNG device falls back to another font; cairo_pdf does not.
#
# Changes (content otherwise identical):
#   S5  dashed-line label "p = 5x10^-8"   -> plotmath superscript with U+2212 minus
#   S5  caption "(min p ~ 10^-3)"         -> "(min p ~ 0.0014)"  [exact min p = 1.44e-3,
#                                             matches the legend's "minimum p ~ 1.4x10^-3"]
#   S7  p-value labels "p = 1.4x10^-3"    -> plotmath superscript with U+2212 minus
#
# Usage: Rscript -e ".libPaths(...); source('FIX_FIGS5_S7_vector.R')" <project_dir> <out_dir>
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- args[1]
OUT <- if (length(args) >= 2) args[2] else
  file.path(PROJECT_DIR, "Final Figures and data", "figures", "supplementary")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(patchwork); library(dplyr); library(grid)
})
set.seed(42)

FONT <- "Arial"
REV  <- file.path(PROJECT_DIR, "results", "reviewer_response")

COL_PRIMARY   <- "#1B4F8A"
COL_SECONDARY <- "#4A90D9"
COL_NULL      <- "#AAAAAA"
COL_CONTROL   <- "#2E7D32"
COL_RISK      <- "#C62828"
MINUS <- "−"   # true minus sign: present in Arial, unlike U+207B

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

# plotmath label "p = m×10^e" with a real minus in the exponent
sci_label <- function(mantissa, exponent) {
  sprintf('p == "%s×10"^"%s"', mantissa,
          sub("^-", MINUS, as.character(exponent)))
}

# ============================================================
# FigS5 — Colocalization regional plot at TF locus (from R1_FIGURES.R)
# ============================================================

cat("FigS5: colocalization regional plot\n")
coloc <- read.csv(file.path(REV, "D_coloc_TF_snp_posteriors.csv"),
                   stringsAsFactors = FALSE)
coloc$p_tf  <- 2 * pnorm(-abs(coloc$z.df1))
coloc$p_tau <- 2 * pnorm(-abs(coloc$z.df2))
coloc$logp_tf  <- -log10(pmax(coloc$p_tf,  1e-300))
coloc$logp_tau <- -log10(pmax(coloc$p_tau, 1e-300))
coloc$mb <- coloc$position / 1e6

TF_LO <- 133.464073
TF_HI <- 133.494388

top_h4 <- coloc[order(-coloc$SNP.PP.H4)[1], ]
cat(sprintf("  Top PP.H4 SNP: %s at %.4f Mb (PP.H4=%.3f)\n",
            top_h4$snp, top_h4$mb, top_h4$SNP.PP.H4))
cat(sprintf("  Min tau p in window: %.3g\n", min(coloc$p_tau, na.rm = TRUE)))

pp_h0 <- 0.000; pp_h1 <- 0.536; pp_h2 <- 0.000; pp_h3 <- 0.017; pp_h4 <- 0.447

coloc$is_top <- coloc$snp == top_h4$snp

y_lim_tf  <- c(0, max(coloc$logp_tf,  na.rm = TRUE) * 1.20)
y_lim_tau <- c(0, max(4.3, max(coloc$logp_tau, na.rm = TRUE) * 1.20))

make_track <- function(logp, top_row, title, y_lab, y_lim, annotation = NULL) {
  ggplot(coloc, aes(x = mb, y = logp)) +
    annotate("rect", xmin = TF_LO, xmax = TF_HI,
             ymin = -Inf, ymax = Inf,
             fill = "#E3F2FD", alpha = 0.6) +
    annotate("text", x = (TF_LO + TF_HI) / 2, y = y_lim[2] * 0.975,
             label = "TF", family = FONT, fontface = "italic",
             size = 2.5, color = COL_PRIMARY) +
    geom_point(color = ifelse(coloc$is_top, COL_RISK, COL_PRIMARY),
               alpha = ifelse(coloc$is_top, 1, 0.55),
               size  = ifelse(coloc$is_top, 2.5, 1)) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed",
               color = "#888888", linewidth = 0.3) +
    ### FIX: plotmath exponent with U+2212 (was "p = 5×10⁻⁸", U+207B not in Arial)
    annotate("text", x = min(coloc$mb) + 0.05,
             y = -log10(5e-8) + y_lim[2] * 0.03,
             label = sci_label("5", -8), parse = TRUE,
             ### FIX (2026-09-13): size 2 -> 2.4 (this page is scaled to ~91% in the
             ### combined supplementary PDF, which left the exponent under 4 pt)
             family = FONT, size = 2.4, color = "grey40") +
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
### R2 FIX: the regional lead variant in the tau GWAS is rs3811658, not
### rs8177240 (which ranks 3rd). The revised Results say so explicitly, so
### panel B must label it. Added as extra layers; make_track is unchanged.
tau_lead <- coloc[which.min(coloc$p_tau), c("snp", "mb", "logp_tau")]
names(tau_lead)[3] <- "logp"
cat(sprintf("  Tau regional lead: %s at %.4f Mb (p=%.3g)
",
            tau_lead$snp, tau_lead$mb, min(coloc$p_tau, na.rm = TRUE)))
stopifnot(tau_lead$snp == "rs3811658")

p_tau <- make_track(coloc$logp_tau, top_h4[, c("snp","mb","logp_tau")] %>%
                       rename(logp = logp_tau),
                     "B  Circulating total-tau (Sarnowski 2022, n = 14,721)",
                     expression(-log[10] * "(p)"),
                     y_lim_tau) +
  geom_point(data = tau_lead, aes(x = mb, y = logp),
             color = COL_CONTROL, size = 2.5, inherit.aes = FALSE) +
  geom_text_repel(data = tau_lead, aes(x = mb, y = logp, label = snp),
                  inherit.aes = FALSE, color = COL_CONTROL, size = 2.5,
                  family = FONT, fontface = "bold",
                  segment.color = COL_CONTROL, segment.size = 0.3,
                  nudge_y = y_lim_tau[2] * 0.10, nudge_x = 0.075,
                  box.padding = 0.5, seed = 42) +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005))) +
  labs(x = "Chromosome 3 position (Mb, GRCh37)")

figs5 <- p_tf / p_tau + plot_layout(heights = c(1, 1)) +
  plot_annotation(caption = paste0(
    "PP.H4 = 0.447 is below the conventional 0.8 threshold; the conditional H4/(H3+H4) = 96.3% is reported because the\n",
    ### FIX: "(min p ≈ 10⁻³)" -> exact value, no superscript minus needed
    "residual mass sits on PP.H1 (transferrin-only), reflecting limited power in the tau GWAS at this locus (min p ≈ 0.0014),\n",
    "not evidence against colocalization (PP.H3, distinct causal variants, ≈ 0.02).
",
    "In panel B the regional lead variant is rs3811658 (green); rs8177240 (red) is the top PP.H4 variant and ranks third in the
",
    "tau GWAS. The two are 849 bp apart and in near-complete linkage disequilibrium (r² = 0.974)."),
    theme = theme(plot.caption = element_text(family = FONT, size = 7.2,
                  hjust = 0, color = "grey30", margin = margin(t = 4))))
save_fig(file.path(OUT, "FigS5_coloc_TF_locus"), figs5, mm2in(175), mm2in(140))

# ============================================================
# FigS7 — Extended MVMR forest (from R1_FIGURES.R)
# ============================================================

cat("FigS7: extended MVMR forest\n")
mv2 <- read.csv(file.path(REV, "E_MVMR_extended.csv"), stringsAsFactors = FALSE)
mv4 <- read.csv(file.path(REV, "F_MVMR_4way.csv"), stringsAsFactors = FALSE)

if (file.exists(file.path(PROJECT_DIR, "data", "FINAL_analysis_results.RData"))) {
  e <- new.env()
  load(file.path(PROJECT_DIR, "data", "FINAL_analysis_results.RData"), envir = e)
  pr <- e$primary_results
  univar_row <- pr[pr$method == "Inverse variance weighted" &
                    pr$pair == "transferrin_circ_total_tau", ]
  univar_b   <- univar_row$b
  univar_se  <- univar_row$se
  univar_p   <- univar_row$pval
} else {
  univar_b <- -0.039; univar_se <- 0.013; univar_p <- 2.88e-3
}

mv3 <- mv4[mv4$model == "3-way (transferrin + ferritin + tsat)", ]
mv3$exposure_name <- c("transferrin", "ferritin", "tsat")

fdat <- rbind(
  data.frame(model = "Univariable IVW", exposure = "Transferrin",
             beta = univar_b, se = univar_se, pval = univar_p, cond_F = NA),
  data.frame(model = "2-way MVMR (transferrin + serum iron)", exposure = "Transferrin",
             beta = mv2$beta[mv2$exposure == "transferrin"],
             se   = mv2$se[mv2$exposure == "transferrin"],
             pval = mv2$pval[mv2$exposure == "transferrin"],
             cond_F = mv2$cond_F[mv2$exposure == "transferrin"]),
  data.frame(model = "2-way MVMR (transferrin + serum iron)", exposure = "Serum iron",
             beta = mv2$beta[mv2$exposure == "serum_iron"],
             se   = mv2$se[mv2$exposure == "serum_iron"],
             pval = mv2$pval[mv2$exposure == "serum_iron"],
             cond_F = mv2$cond_F[mv2$exposure == "serum_iron"]),
  data.frame(model = "3-way MVMR (transferrin + ferritin + TSAT)", exposure = "Transferrin",
             beta = mv3$beta[mv3$exposure_name == "transferrin"],
             se   = mv3$se[mv3$exposure_name == "transferrin"],
             pval = mv3$pval[mv3$exposure_name == "transferrin"],
             cond_F = mv3$cond_F[mv3$exposure_name == "transferrin"]),
  data.frame(model = "3-way MVMR (transferrin + ferritin + TSAT)", exposure = "Ferritin",
             beta = mv3$beta[mv3$exposure_name == "ferritin"],
             se   = mv3$se[mv3$exposure_name == "ferritin"],
             pval = mv3$pval[mv3$exposure_name == "ferritin"],
             cond_F = mv3$cond_F[mv3$exposure_name == "ferritin"]),
  data.frame(model = "3-way MVMR (transferrin + ferritin + TSAT)", exposure = "TSAT",
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

### FIX: plotmath labels. Same rule as before (2 significant figures; scientific
### below 0.01, plain decimal at/above), but the exponent uses U+2212 so the
### glyph exists in Arial. Plain decimals are quoted so "0.80" keeps its zero.
fmt_p <- function(p) vapply(p, function(x) {
  if (is.na(x)) return('""')
  if (x < 0.01) { e <- floor(log10(x)); m <- x / 10^e
    sci_label(formatC(m, format = "f", digits = 1), e) }
  else if (x >= 0.1) sprintf('p == "%.2f"', x)
  else sprintf('p == "%.3f"', x)
}, character(1))
fdat$plab <- fmt_p(fdat$pval)
# fontface is ignored for parsed labels, so significant rows get plotmath bold()
fdat$plab <- ifelse(fdat$sig, paste0("bold(", fdat$plab, ")"), fdat$plab)
cat("  p labels:", paste(fdat$plab, collapse = " | "), "\n")
fdat$fann <- ifelse(is.na(fdat$cond_F), "", sprintf("cond-F = %.1f", fdat$cond_F))

figs7 <- ggplot(fdat, aes(x = beta, y = row_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25,
                 color = fdat$color, linewidth = 0.6) +
  geom_point(color = fdat$color, size = 3,
             shape = ifelse(fdat$is_transferrin, 15, 18)) +
  geom_text(aes(x = 0.095, label = plab), parse = TRUE,
            ### FIX (2026-09-13): 2.4 -> 2.7 so the exponent clears 5 pt at final size
            hjust = 0, size = 2.7, family = FONT,
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

save_fig(file.path(OUT, "FigS7_MVMR_extended"), figs7, mm2in(190), mm2in(105))
