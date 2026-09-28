# FIX_FIGS11_R2_style.R
#
# Figure S11 restyled to match S1-S10. The default coloc::sensitivity() output
# is base graphics with "trait 1"/"trait 2" labels and no panel letters, which
# reads as a screenshot beside the other supplementary figures.
#
# Header constants (FONT, colours, theme_forest, save_fig, mm2in) are copied
# verbatim from FIX_FIGS9_notation.R so the house style is identical.
#
# Inputs: results/R2/coloc_sensitivity_curve.csv, results/R2/coloc_prior_grid.csv
# Usage:  Rscript ... <PROJECT_DIR> <OUTPUT_DIR>

args <- commandArgs(trailingOnly = TRUE)
PROJECT_DIR <- if (length(args) >= 1) args[1] else normalizePath("..", mustWork = TRUE)
OUT <- if (length(args) >= 2) args[2] else stop("give the output directory as arg 2")

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})
set.seed(42)

FONT  <- "Arial"
MINUS <- "−"
R2DIR <- file.path(PROJECT_DIR, "results", "R2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

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
      plot.title = element_text(face = "bold", size = base, margin = margin(b = 6))
    )
}

# ---------------------------------------------------------------- data
curve <- read.csv(file.path(R2DIR, "coloc_sensitivity_curve.csv"))
grid  <- read.csv(file.path(R2DIR, "coloc_prior_grid.csv"))
cat(sprintf("  curve rows: %d (p12 %.2e to %.2e)\n", nrow(curve),
            min(curve$p12), max(curve$p12)))

P12_DEFAULT <- 1e-5

# Thresholds are INTERPOLATED in log10(p12). Taking the first grid point at
# which a condition already holds overstates it -- that is what produced the
# earlier "5e-5" and "1.29e-5" claims.
log_root <- function(y) {
  i <- which(diff(sign(y)) != 0)[1]
  if (is.na(i)) return(NA_real_)
  x1 <- log10(curve$p12[i]); x2 <- log10(curve$p12[i + 1])
  10 ^ (x1 + (0 - y[i]) * (x2 - x1) / (y[i + 1] - y[i]))
}
P12_CROSS <- log_root(curve$PP.H4.abf - curve$PP.H1.abf)   # H4 overtakes H1
RULE_MIN  <- log_root(curve$PP.H4.abf - 0.5)               # H4 crosses 0.5
cat(sprintf("  H4 overtakes H1 at p12 = %.4e  (%.1f%% above default)\n",
            P12_CROSS, 100 * (P12_CROSS / 1e-5 - 1)))
cat(sprintf("  H4 crosses 0.5  at p12 = %.4e  (%.1f%% above default)\n",
            RULE_MIN, 100 * (RULE_MIN / 1e-5 - 1)))
dflt <- curve[which.min(abs(curve$p12 - P12_DEFAULT)), ]
cat(sprintf("  at default p12=1e-5: H1=%.3f H4=%.3f\n",
            dflt$PP.H1.abf, dflt$PP.H4.abf))

sup <- function(v) paste0('10^"', MINUS, abs(round(log10(v))), '"')

# facet label: 5e-05 -> "p12 = 5x10^-5", marking the default
fmt_p12 <- function(v) {
  e <- floor(log10(v)); m <- v / 10 ^ e
  mant <- ifelse(abs(m - 1) < 1e-9, "", paste0(round(m), "\u00d7"))
  paste0("p12 = ", mant, "10", MINUS, chartr("0123456789", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079", sub("-", "", as.character(e))),
         ifelse(abs(v - P12_DEFAULT) < 1e-12, "  (default)", ""))
}
xbreaks <- 10 ^ seq(-8, -4)

# ---------------------------------------------------------------- panel A
long <- curve %>%
  select(p12, H1 = PP.H1.abf, H3 = PP.H3.abf, H4 = PP.H4.abf) %>%
  pivot_longer(-p12, names_to = "hyp", values_to = "pp") %>%
  mutate(hyp = factor(hyp, levels = c("H4", "H1", "H3"),
                      labels = c("PP.H4  shared causal variant",
                                 "PP.H1  transferrin only",
                                 "PP.H3  distinct causal variants")))

pA <- ggplot(long, aes(p12, pp, colour = hyp)) +
  annotate("rect", xmin = RULE_MIN, xmax = max(curve$p12),
           ymin = -Inf, ymax = Inf, fill = "#E8F5E9", alpha = 0.7) +
  annotate("text", x = sqrt(RULE_MIN * max(curve$p12)), y = 0.97,
           label = "PP.H4 > 0.5", family = FONT, size = 2.3,
           colour = COL_CONTROL, fontface = "bold") +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             colour = "#888888", linewidth = 0.3) +
  geom_vline(xintercept = P12_DEFAULT, linetype = "dashed",
             colour = COL_RISK, linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  annotate("text", x = P12_DEFAULT * 0.82, y = 0.30,
           label = "default prior", family = FONT, size = 2.3,
           colour = COL_RISK, angle = 90, hjust = 0.5) +
  geom_point(data = data.frame(p12 = rep(P12_DEFAULT, 2),
                               pp  = c(dflt$PP.H4.abf, dflt$PP.H1.abf),
                               hyp = factor(c("PP.H4  shared causal variant",
                                              "PP.H1  transferrin only"),
                                            levels = levels(long$hyp))),
             size = 1.9) +
  scale_colour_manual(values = c(COL_CONTROL, COL_PRIMARY, COL_NULL), name = NULL) +
  scale_x_log10(breaks = xbreaks, labels = parse(text = sup(xbreaks)),
                expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "A  Posterior probabilities across the shared-variant prior",
       x = expression("Prior probability of a shared causal variant, " * p[12]),
       y = "Posterior probability") +
  theme_forest(base = 8) +
  theme(legend.position = c(0.015, 0.99), legend.justification = c(0, 1),
        legend.key.height = unit(9, "pt"), legend.key.width = unit(14, "pt"),
        legend.text = element_text(size = 6.6),
        legend.background = element_rect(fill = "#FFFFFFCC", colour = NA),
        axis.text = element_text(size = 7))

# ---------------------------------------------------------------- panel B
sel <- grid %>% filter(axis == "p12") %>%
  filter(p12 %in% c(1e-6, 1e-5, 5e-5, 1e-4))
bars <- sel %>%
  select(p12, PP.H0, PP.H1, PP.H2, PP.H3, PP.H4) %>%
  pivot_longer(-p12, names_to = "hyp", values_to = "pp") %>%
  mutate(hyp = factor(hyp, levels = paste0("PP.H", 0:4)),
         lab = factor(fmt_p12(p12), levels = fmt_p12(sel$p12)))

pB <- ggplot(bars, aes(hyp, pp, fill = hyp)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = ifelse(pp >= 0.001, sprintf("%.2f", pp), "")),
            vjust = -0.45, size = 2.1, family = FONT, colour = "grey25") +
  facet_wrap(~ lab, nrow = 1) +
  scale_fill_manual(values = c(COL_NULL, COL_PRIMARY, COL_NULL,
                               COL_SECONDARY, COL_CONTROL), guide = "none") +
  scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "B  All five unconditional posteriors at selected priors",
       x = NULL, y = "Posterior probability") +
  theme_forest(base = 8) +
  theme(axis.text.x = element_text(size = 6.2, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 7, family = FONT,
                                  margin = margin(b = 3, t = 1)),
        panel.spacing.x = unit(7, "pt"))

# ---------------------------------------------------------------- assemble
cap <- paste0(
  "Sensitivity of the TF-locus colocalization to the shared-variant prior p12. ",
  "Panel A: coloc.abf re-run across p12, with p1 = p2 = 1\u00d710", MINUS, "\u2074 held fixed.\n",
  "At the default p12 = 1\u00d710", MINUS, "\u2075 (dashed line) the transferrin-only hypothesis ",
  "remains more probable than a shared causal variant (0.54 vs 0.45),\n",
  "and the conventional rule PP.H4 > 0.5 is not met. The ordering is finely ",
  "balanced: PP.H4 overtakes PP.H1 at p12 \u2248 1.2\u00d710", MINUS, "\u2075, only about\n",
  "20% above the default, and crosses 0.5 shortly after, at p12 \u2248 1.24\u00d710",
  MINUS, "\u2075 (shaded). Both values are interpolated, not grid points. ",
  "Panel B: all five unconditional posteriors at four priors;\n",
  "PP.H0 and PP.H2 are below 0.001 throughout. Numerical values are given in ",
  "Supplementary Table S8.")

figs11 <- pA / pB + plot_layout(heights = c(1.15, 1)) +
  plot_annotation(caption = cap,
    theme = theme(plot.caption = element_text(family = FONT, size = 7.2,
                  hjust = 0, colour = "grey30", margin = margin(t = 4))))

save_fig(file.path(OUT, "FigS11_coloc_sensitivity"), figs11,
         mm2in(175), mm2in(135))
