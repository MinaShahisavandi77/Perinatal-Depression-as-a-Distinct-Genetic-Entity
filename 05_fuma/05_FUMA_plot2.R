# =============================================================================
# MAGMA Gene-Set Enrichment — Publication-Quality Figures (Panel A & B)
# Traits: MDD_shared  |  PPD_unique
# Output: single 600 dpi TIFF with panels A & B side-by-side
#
# Paths come from config/config.R — edit them there, not here.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────
repos <- "https://cloud.r-project.org"
for (p in c("data.table", "dplyr", "ggplot2", "patchwork", "ggtext", "here")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggtext)   # rich-text axis / title annotations

# ── Config ────────────────────────────────────────────────────────────────────
source(here::here("GSEM","scripts", "config.R"))

# ── Paths ────────────────────────────────────────────────────────────────────
mdd_gsa <- file.path(DIR_FUMA_MDD, "magma.gsa.out")
ppd_gsa <- file.path(DIR_FUMA_PPD, "magma.gsa.out")
outdir  <- file.path(DIR_FUMA, "PATHWAY_PLOTS")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Global settings ───────────────────────────────────────────────────────────
trait_mdd  <- "MDD-shared"
trait_ppd  <- "PPD-unique"
top_n      <- 30          # top gene sets per panel
dpi        <- 600
FDR_THRESH <- 0.05
NOM_THRESH <- 0.05        # nominal p threshold for colouring

# ── Colour palette (colour-blind safe, publication ready) ────────────────────
COL_SIG  <- "#B2182B"   # FDR-significant  (deep red)
COL_NOM  <- "#4393C3"   # nominally significant (steel blue)
COL_NS   <- "#BABABA"   # not significant (grey)
COL_ANNO <- "#252525"   # annotation text

# ── Publication ggplot2 theme ─────────────────────────────────────────────────
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      # Axes
      axis.line        = element_line(colour = "#333333", linewidth = 0.45),
      axis.ticks       = element_line(colour = "#333333", linewidth = 0.35),
      axis.text.x      = element_text(colour = COL_ANNO, size = base_size - 1),
      axis.text.y      = element_text(colour = COL_ANNO, size = base_size - 1,
                                       hjust = 1),
      axis.title       = element_text(colour = COL_ANNO, size = base_size,
                                       face = "bold"),
      # Panel
      panel.grid.major.x = element_line(colour = "#E8E8E8", linewidth = 0.3,
                                          linetype = "dashed"),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      # Legend
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.key.size   = unit(0.45, "cm"),
      legend.text       = element_text(size = base_size - 1, colour = COL_ANNO),
      legend.title      = element_text(size = base_size - 1, colour = COL_ANNO,
                                        face = "bold"),
      legend.background = element_blank(),
      # Titles
      plot.title    = element_text(face = "bold", size = base_size + 1,
                                    colour = COL_ANNO, margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 1, colour = "#555555",
                                    margin = margin(b = 6)),
      plot.margin   = margin(8, 10, 6, 4)
    )
}

# ── Helper: read & annotate MAGMA GSA output ─────────────────────────────────
read_magma_gsa <- function(path, trait) {
  dt      <- fread(path, data.table = FALSE, comment.char = "#")
  stopifnot("P" %in% names(dt))

  set_col <- if ("VARIABLE" %in% names(dt)) "VARIABLE" else names(dt)[1]

  out <- dt %>%
    transmute(
      trait  = trait,
      SET    = as.character(.data[[set_col]]),
      NGENES = if ("NGENES" %in% names(dt)) as.integer(NGENES) else NA_integer_,
      BETA   = if ("BETA"   %in% names(dt)) as.numeric(BETA)   else NA_real_,
      SE     = if ("SE"     %in% names(dt)) as.numeric(SE)      else NA_real_,
      P      = as.numeric(P)
    ) %>%
    filter(is.finite(P), P > 0, P <= 1) %>%
    mutate(FDR = p.adjust(P, "BH")) %>%
    arrange(P)

  out
}

# ── Helper: clean gene-set labels for publication ────────────────────────────
clean_label <- function(x) {
  x <- gsub("_", " ", x)                 # underscores → spaces
  x <- gsub("GOBP|GOCC|GOMF|KEGG|REACTOME|HP", "", x)  # strip DB prefix
  x <- trimws(x)
  x <- tolower(x)
  # Title-case first letter only
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

# ── Core plot function ────────────────────────────────────────────────────────
make_panel <- function(gsa_df, trait, panel_label, top_n = 30) {

  d <- gsa_df %>%
    arrange(P) %>%
    slice_head(n = min(top_n, nrow(gsa_df))) %>%
    mutate(
      mlog10p  = -log10(P),
      label    = clean_label(SET),
      sig_cat  = case_when(
        FDR < FDR_THRESH ~ "FDR < 0.05",
        P   < NOM_THRESH ~ "Nominal p < 0.05",
        TRUE             ~ "Not significant"
      ),
      sig_cat = factor(sig_cat,
                       levels = c("FDR < 0.05",
                                  "Nominal p < 0.05",
                                  "Not significant")),
      label = factor(label, levels = rev(label))   # bottom-up ordering
    )

  # FDR threshold line on x-axis
  fdr_line <- -log10(max(d$P[d$FDR < FDR_THRESH], na.rm = TRUE) + 1e-30)

  # Number of significant sets for subtitle
  n_fdr <- sum(d$FDR < FDR_THRESH)
  n_nom <- sum(d$P   < NOM_THRESH & d$FDR >= FDR_THRESH)

  ggplot(d, aes(x = mlog10p, y = label, fill = sig_cat)) +
    # FDR threshold reference line
    geom_vline(xintercept = fdr_line, linetype = "dashed",
               colour = COL_SIG, linewidth = 0.45, alpha = 0.7) +
    # Bars
    geom_col(width = 0.72, colour = NA) +
    # NGENES annotation at bar end (if available)
    {if (!all(is.na(d$NGENES)))
      geom_text(aes(label = ifelse(!is.na(NGENES),
                                   paste0("n=", NGENES), "")),
                hjust = -0.15, size = 2.5, colour = "#555555")
    } +
    # Colour scale
    scale_fill_manual(
      name   = "Significance",
      values = c("FDR < 0.05"         = COL_SIG,
                 "Nominal p < 0.05"   = COL_NOM,
                 "Not significant"    = COL_NS),
      drop   = FALSE
    ) +
    # Expand x so gene-count labels aren't clipped
    scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
    # Labels
    labs(
      title    = paste0("**", panel_label, "** ", trait),
      subtitle = paste0("Top ", nrow(d), " gene sets  |  ",
                        n_fdr, " FDR-significant,  ", n_nom, " nominally significant"),
      x        = expression(-log[10](italic(p))),
      y        = NULL
    ) +
    theme_pub(base_size = 10) +
    theme(
      plot.title = element_markdown(face = "bold", size = 11),
      # Push legend to bottom-right so it doesn't crowd bars
      legend.position  = "bottom",
      legend.justification = "left"
    )
}

# ── Read data ─────────────────────────────────────────────────────────────────
cat("Reading MAGMA GSA results...\n")
mdd <- read_magma_gsa(mdd_gsa, trait_mdd)
ppd <- read_magma_gsa(ppd_gsa, trait_ppd)

# Save ranked tables
fwrite(mdd, file.path(outdir, paste0(gsub("-","_",trait_mdd), "_magma_gsa_ranked.tsv")), sep = "\t")
fwrite(ppd, file.path(outdir, paste0(gsub("-","_",trait_ppd), "_magma_gsa_ranked.tsv")), sep = "\t")

# ── Build panels ──────────────────────────────────────────────────────────────
p_mdd <- make_panel(mdd, trait_mdd, "A", top_n)
p_ppd <- make_panel(ppd, trait_ppd, "B", top_n)

# ── Combine into publication-ready A|B grid ───────────────────────────────────
# Shared legend is extracted from panel A via patchwork's guide_area()
combined <- (p_mdd | p_ppd) +
  plot_layout(
    ncol    = 2,
    guides  = "collect",   # single shared legend
    widths  = c(1, 1)
  ) +
  plot_annotation(
    title   = "MAGMA gene-set enrichment analysis",
    caption = paste0("FDR correction: Benjamini–Hochberg. ",
                     "Dashed line: FDR = 0.05 threshold. ",
                     "Bar labels: number of genes per set."),
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 13,
                                   colour = COL_ANNO, hjust = 0,
                                   margin = margin(b = 4)),
      plot.caption = element_text(size = 8, colour = "#777777",
                                   hjust = 0, margin = margin(t = 6)),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) &
  theme(legend.position = "bottom")

# ── Dynamic figure height (scales with number of gene sets) ──────────────────
n_rows    <- min(top_n, max(nrow(mdd), nrow(ppd)))
height_in <- max(10, 0.28 * n_rows + 4)
width_in  <- 18   # two panels side-by-side

# ── Save combined figure ──────────────────────────────────────────────────────
out_tiff <- file.path(outdir, paste0("Fig_MAGMA_pathway_enrichment_AB_top", top_n, ".tiff"))

ggsave(
  filename    = out_tiff,
  plot        = combined,
  device      = "tiff",
  dpi         = dpi,
  compression = "lzw",
  width       = width_in,
  height      = height_in,
  units       = "in"
)

cat("\n✓ Combined A/B figure saved →", out_tiff, "\n")

# ── Also save individual panels (useful for supplementary) ───────────────────
ggsave(
  file.path(outdir, paste0(gsub("-","_",trait_mdd), "_panelA_top", top_n, ".tiff")),
  p_mdd, device = "tiff", dpi = dpi, compression = "lzw",
  width = 9, height = height_in, units = "in"
)
ggsave(
  file.path(outdir, paste0(gsub("-","_",trait_ppd), "_panelB_top", top_n, ".tiff")),
  p_ppd, device = "tiff", dpi = dpi, compression = "lzw",
  width = 9, height = height_in, units = "in"
)

cat("✓ Individual panels saved to", outdir, "\n")

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n── Significant pathways (FDR < 0.05) ──────────────────────────\n")
cat(sprintf("  %-20s : %d\n", trait_mdd, sum(mdd$FDR < FDR_THRESH)))
cat(sprintf("  %-20s : %d\n", trait_ppd, sum(ppd$FDR < FDR_THRESH)))
cat("────────────────────────────────────────────────────────────────\n")