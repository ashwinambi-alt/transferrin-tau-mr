# FIX_FIGS8_keysize.R (2026-09-13): lines 20-61 and 366-458 of R1_FIGURES.R extracted
# verbatim, then edited ONLY to enlarge the colour-key text (see ### FIX).
# Output goes to a staging directory (2nd argument).
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
           family = FONT, fontface = "italic", size = 2.6,
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
    ### FIX (2026-09-13): 2.2 -> 2.7; the page is scaled to 80% in the combined
    ### supplementary PDF, which left this key at ~5 pt
    family = FONT, size = 2.7, hjust = 0, vjust = 1, lineheight = 1.15,
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
