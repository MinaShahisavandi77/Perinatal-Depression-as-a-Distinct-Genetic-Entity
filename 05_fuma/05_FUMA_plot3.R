# =============================================================================
# Gene-based MAGMA Manhattan (per trait) + top-gene labels
# Author : M. Shahisavandi
#
# Reads the MAGMA gene-level results for each trait, draws a Manhattan with a
# Bonferroni line, labels the top genes, and writes a TIFF + a labelled-gene
# table per trait. Paths come from config/config.R — edit them there, not here.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
repos <- "https://cloud.r-project.org"
for (p in c("data.table", "dplyr", "ggplot2", "ggrepel", "here")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ---- Config -----------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R"))

# ---- Inputs / outputs -------------------------------------------------------
genes_file <- "magma.genes.out"
mdd_genes  <- file.path(DIR_FUMA_MDD, genes_file)
ppd_genes  <- file.path(DIR_FUMA_PPD, genes_file)

trait_mdd <- "MDD_Shared"
trait_ppd <- "PPD_Unique"

outdir <- file.path(DIR_FUMA, "GENE_MANHATTAN")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Parameters -------------------------------------------------------------
dpi         <- 600
label_top_n <- 20

# =============================================================================
# Read a MAGMA gene-level result table
# =============================================================================
read_magma_genes <- function(path, trait) {
  dt <- fread(path, data.table = FALSE, comment.char = "#")
  stopifnot(all(c("GENE", "CHR", "START", "STOP", "P") %in% names(dt)))

  dt %>%
    transmute(
      trait  = trait,
      SYMBOL = as.character(SYMBOL),
      CHR    = as.integer(CHR),
      START  = as.integer(START),
      STOP   = as.integer(STOP),
      P      = as.numeric(P)
    ) %>%
    filter(is.finite(P), P > 0, P <= 1,
           !is.na(CHR), CHR >= 1, CHR <= 22,
           is.finite(START), is.finite(STOP)) %>%
    mutate(
      MID  = floor((START + STOP) / 2),
      logp = -log10(P),
      FDR  = p.adjust(P, "BH")
    ) %>%
    arrange(CHR, MID)
}

# =============================================================================
# Gene-based Manhattan with Bonferroni line and top-gene labels
# =============================================================================
make_gene_manhattan <- function(df, trait, out_tiff, label_top_n = 20,
                                col_dark = "#08306B", col_light = "#6BAED6",
                                width_in = 12, height_in = 5.5, dpi = 600) {

  # cumulative positions
  chr_sizes <- df %>%
    group_by(CHR) %>%
    summarise(chr_len = max(MID, na.rm = TRUE), .groups = "drop") %>%
    arrange(CHR) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0))

  d <- df %>%
    inner_join(chr_sizes, by = "CHR") %>%
    mutate(
      x = MID + chr_start,
      chr_col = factor(CHR %% 2)
    )

  axis_df <- chr_sizes %>% mutate(center = chr_start + chr_len / 2)

  # top genes to label (by p-value)
  labs_df <- d %>% arrange(P) %>% slice_head(n = min(label_top_n, nrow(d)))

  # Bonferroni line (per trait, based on genes in this file)
  bonf <- 0.05 / nrow(d)

  g <- ggplot(d, aes(x = x, y = logp, color = chr_col)) +
    geom_point(size = 0.6, alpha = 0.85) +
    scale_color_manual(values = c("0" = col_dark, "1" = col_light), guide = "none") +
    geom_hline(yintercept = -log10(bonf), linetype = 2, linewidth = 0.55, color = "grey10") +
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$CHR,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    theme_classic(base_size = 12) +
    labs(
      title    = paste0("MAGMA gene-based Manhattan: ", trait),
      subtitle = paste0("Dashed line = Bonferroni (0.05 / ", nrow(d), " genes) = ", signif(bonf, 3)),
      x = "Chromosome",
      y = expression(-log[10](p))
    ) +
    ggrepel::geom_text_repel(
      data = labs_df,
      aes(label = SYMBOL),
      size = 3.0,
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15,
      segment.size = 0.25,
      segment.alpha = 0.6,
      max.overlaps = Inf
    )

  ggsave(out_tiff, g, device = "tiff",
         width = width_in, height = height_in, units = "in",
         dpi = dpi, compression = "lzw")

  list(plot = g, bonf = bonf, labels = labs_df)
}

# =============================================================================
# Run
# =============================================================================
mdd <- read_magma_genes(mdd_genes, trait_mdd)
ppd <- read_magma_genes(ppd_genes, trait_ppd)

res_mdd <- make_gene_manhattan(
  mdd, trait_mdd,
  out_tiff = file.path(outdir, paste0(trait_mdd, "_magma_gene_manhattan_top", label_top_n, ".tiff")),
  label_top_n = label_top_n,
  dpi = dpi
)

res_ppd <- make_gene_manhattan(
  ppd, trait_ppd,
  out_tiff = file.path(outdir, paste0(trait_ppd, "_magma_gene_manhattan_top", label_top_n, ".tiff")),
  label_top_n = label_top_n,
  dpi = dpi
)

# Save which genes were labelled
fwrite(res_mdd$labels %>% select(SYMBOL, CHR, START, STOP, P, FDR),
       file.path(outdir, paste0(trait_mdd, "_top", label_top_n, "_labeled_genes.tsv")),
       sep = "\t")
fwrite(res_ppd$labels %>% select(SYMBOL, CHR, START, STOP, P, FDR),
       file.path(outdir, paste0(trait_ppd, "_top", label_top_n, "_labeled_genes.tsv")),
       sep = "\t")

cat("Gene-based Manhattans + labelled-gene tables written to:", outdir, "\n")