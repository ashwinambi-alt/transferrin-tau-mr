# FIX_FIGS9_R2_caption.R (R2): copy of FIX_FIGS9_notation.R with ONLY the
# caption changed, to match the revised Results. Original header follows.
# FIX_FIGS9_notation.R (2026-09-13): lines 20-61, 261-274 and 460-518 of R1_FIGURES.R,
# extracted verbatim, then edited ONLY for label notation (see ### FIX comments).
# Output goes to a staging directory (2nd argument), not the supplementary folder.
args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)

suppressPackageStartupMessages({
  library(ggplot2); library(ggrepel); library(patchwork); library(dplyr)
  library(grid); library(ggforce); library(ggdag); library(dagitty)
})
set.seed(42)

FONT <- "Arial"
BASE <- file.path(PROJECT_DIR, "Final Figures and data")
### FIX: write to a staging directory given as the 2nd argument
SUPP <- if (length(args) >= 2) args[2] else stop("give the staging directory as the 2nd argument")
MINUS <- "−"   # true minus (U+2212); Arial has no superscript minus
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
### FIX: "2.88e-03" -> plotmath bold "2.88×10⁻³" (same 3 significant figures)
comp$plab  <- vapply(comp$pval, function(x) {
  e <- floor(log10(x)); m <- x / 10^e
  sprintf('bold(p == "%s×10"^"%s%d")', formatC(m, format = "f", digits = 2), MINUS, -e)
}, character(1))
comp$row   <- rev(seq_len(nrow(comp)))

BONF <- 0.05 / 23   # 0.00217

figs9 <- ggplot(comp, aes(x = beta, y = row)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#888888", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2,
                 color = comp$color, linewidth = 0.7) +
  geom_point(color = comp$color,
             shape = c(15, 17, 18), size = c(3.5, 4, 3.5)) +
  geom_text(aes(label = plab), parse = TRUE,
            nudge_y = 0.28, size = 2.5, family = FONT,
            fontface = "bold",
            color = comp$color) +
  geom_text(aes(label = paste0("β = ", sub("^-", MINUS, sprintf("%.3f", beta)),
                                "  (SE ", sprintf("%.3f", se), ")")),
            nudge_y = -0.28, size = 2.2, family = FONT, color = "grey40") +
  scale_y_continuous(breaks = comp$row, labels = comp$label,
                     expand = expansion(add = c(0.7, 0.7))) +
  scale_x_continuous(limits = c(-0.11, 0.02)) +
  labs(title = "Comparison of instrument sets: transferrin → circulating total-tau",
       x = expression("IVW " * beta * " (95% CI)"), y = NULL,
       caption = paste0("Bonferroni threshold for 23 outcome pairs: p < ",
                        sprintf("%.5f", BONF), ".\n",   ### FIX: 0.0022 -> 0.00217
                        "Neither estimate is presented as confirmatory: the cis-only estimate falls below
",
                        "this threshold but was added post hoc and re-estimates the same TF-locus
",
                        "signal, with rs3811658 carrying 89.5% of its inverse-variance weight.")) +
  theme_forest() +
  theme(axis.text.y = element_text(size = 7.5, lineheight = 1.0),
        plot.caption = element_text(size = 7, hjust = 0,
                                      margin = margin(t = 6),
                                      color = "grey30"))

save_fig(file.path(SUPP, "FigS9_cis_TF_comparison"),
         figs9, mm2in(170), mm2in(100))
