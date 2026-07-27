# =============================================================================
# Mirror Manhattan — GSEM MDD_shared (top) vs PPD_resid (bottom)
# Author : M. Shahisavandi
#
# Reads the two GWAS-by-subtraction result tables (stage 04), draws a mirrored
# Manhattan with genome-wide-significant points highlighted and the top hit per
# chromosome labelled, and writes a publication TIFF.
#
# Paths come from config/config.R. 
# =============================================================================

# ---- Packages ---------------------------------------------------------------
repos <- "https://cloud.r-project.org"
cran_pkgs <- c("data.table", "dplyr", "ggplot2", "ggrepel", "here")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ---- Config -----------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R"))
if (!dir.exists(DIR_PLOTS)) dir.create(DIR_PLOTS, recursive = TRUE)

# ---- Inputs (stage 04 result tables) ----------------------------------------
file_trait1 <- file.path(DIR_GSEM, "MDD_F_shared_results.txt")
file_trait2 <- file.path(DIR_GSEM, "PPD_resid_F_results.txt")

trait1_name <- "MDD_F_shared"
trait2_name <- "PPD_resid_F"

# Column names in the result tables
snpcol <- "SNP"
chrcol <- "CHR"
bpcol  <- "BP"
pcol   <- "Pval_Estimate"

# ---- Plot parameters --------------------------------------------------------
dpi         <- 600
highlight_p <- 5e-8
ymax        <- 12          # set NULL to auto-scale

col_dark_blue  <- "#08306B"
col_light_blue <- "#6BAED6"
col_sig_red    <- "#D7191C"

out_tiff <- file.path(
  DIR_PLOTS,
  paste0("Mirror_Manhattan_", trait1_name, "_vs_", trait2_name, "_top1perchr.tiff")
)

# =============================================================================
# Read a GSEM result table into a tidy GWAS frame
# =============================================================================
read_gwas <- function(path,
                      trait_label,
                      snpcol = "SNP", chrcol = "CHR",
                      bpcol  = "BP",  pcol   = "Pval_Estimate") {

  # keep only real records (>=10 tab fields); drop the 4-field warning
  # continuation lines; cut to the 20 columns before `warning`
  cmd <- sprintf("awk -F'\\t' 'NF>=10' %s | cut -f1-20", shQuote(path))
  dt  <- fread(cmd = cmd, header = TRUE, data.table = FALSE, quote = "")

  need <- c(snpcol, chrcol, bpcol, pcol)
  missing <- setdiff(need, names(dt))
  if (length(missing))
    stop("Missing column(s): ", paste(missing, collapse = ", "),
         "\nFound: ", paste(names(dt), collapse = ", "))

  out <- dt[, match(need, names(dt))]
  names(out) <- c("SNP", "CHR", "BP", "P")

  out %>%
    transmute(
      trait = trait_label,
      SNP   = as.character(SNP),
      CHR   = suppressWarnings(as.integer(
                gsub("^chr", "", as.character(CHR), ignore.case = TRUE))),
      BP    = suppressWarnings(as.numeric(BP)),
      P     = suppressWarnings(as.numeric(P))
    ) %>%
    filter(!is.na(SNP), !is.na(CHR), is.finite(BP), is.finite(P),
           CHR >= 1, CHR <= 22, P > 0, P <= 1) %>%
    arrange(CHR, BP)
}

# =============================================================================
# Mirror Manhattan: top trait up, bottom trait down, top hit per chr labelled
# =============================================================================
make_mirror_manhattan <- function(d1, d2, out_tiff, title,
                                   highlight_p = 5e-8, ymax = NULL,
                                   width_in = 12, height_in = 6, dpi = 600,
                                   col_dark = "#08306B", col_light = "#6BAED6",
                                   col_sig = "#D7191C") {

  build_cumpos <- function(d) {
    chr_sizes <- d %>%
      group_by(CHR) %>%
      summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
      arrange(CHR) %>%
      mutate(chr_start = lag(cumsum(chr_len), default = 0))

    d2 <- d %>%
      inner_join(chr_sizes, by = "CHR") %>%
      mutate(BPcum = BP + chr_start)

    axis_df <- chr_sizes %>% mutate(center = chr_start + chr_len / 2)
    list(data = d2, axis = axis_df)
  }

  # Shared x-axis across both traits
  comb <- bind_rows(d1, d2) %>% arrange(CHR, BP)
  pc   <- build_cumpos(comb)
  axis_df <- pc$axis
  key  <- pc$data %>% select(trait, SNP, CHR, BP, BPcum)

  d1p <- d1 %>% inner_join(key, by = c("trait", "SNP", "CHR", "BP")) %>%
    mutate(logp = -log10(P), y =  logp)
  d2p <- d2 %>% inner_join(key, by = c("trait", "SNP", "CHR", "BP")) %>%
    mutate(logp = -log10(P), y = -logp)

  dd <- bind_rows(d1p, d2p) %>% mutate(chr_col = factor(CHR %% 2))
  if (is.null(ymax)) ymax <- ceiling(max(dd$logp, na.rm = TRUE))

  # Labels: strongest genome-wide-significant SNP per chromosome, per trait
  top1_per_chr <- function(dp) {
    dp %>%
      filter(P < highlight_p) %>%
      group_by(CHR) %>%
      slice_min(order_by = P, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(label = SNP)
  }

  labs_df <- bind_rows(top1_per_chr(d1p), top1_per_chr(d2p))

  g <- ggplot(dd, aes(x = BPcum, y = y, color = chr_col)) +
    geom_point(size = 0.55, alpha = 0.80) +
    scale_color_manual(values = c("0" = col_dark, "1" = col_light), guide = "none") +

    # threshold lines (suggestive + genome-wide, both directions)
    geom_hline(yintercept = c(-log10(1e-5),  log10(1e-5)), linetype = 2, linewidth = 0.40, color = "grey45") +
    geom_hline(yintercept = c(-log10(5e-8),  log10(5e-8)), linetype = 2, linewidth = 0.55, color = "grey10") +

    # highlight significant points
    geom_point(
      data = dd %>% filter(P < highlight_p),
      aes(x = BPcum, y = y),
      inherit.aes = FALSE,
      color = col_sig, size = 0.85, alpha = 0.95
    ) +

    # label only the top hit per chromosome (significant only)
    ggrepel::geom_text_repel(
      data = labs_df,
      aes(x = BPcum, y = y, label = label),
      inherit.aes = FALSE,
      size = 3.2, min.segment.length = 0,
      box.padding = 0.25, point.padding = 0.15,
      segment.size = 0.25, segment.alpha = 0.6,
      max.overlaps = Inf
    ) +

    scale_x_continuous(
      breaks = axis_df$center, labels = axis_df$CHR,
      expand = expansion(mult = c(0.005, 0.005))
    ) +
    coord_cartesian(ylim = c(-ymax, ymax)) +
    theme_classic(base_size = 12) +
    labs(
      title = title,
      x = "Chromosome",
      y = expression(paste("Mirror ", -log[10], "(p)"))
    ) +
    annotate("text", x = Inf, y =  ymax * 0.93, label = unique(d1$trait), hjust = 1.05, vjust = 1, size = 4.3) +
    annotate("text", x = Inf, y = -ymax * 0.93, label = unique(d2$trait), hjust = 1.05, vjust = 0, size = 4.3) +
    theme(
      plot.title  = element_text(face = "bold"),
      axis.text.x = element_text(size = 10)
    )

  ggsave(out_tiff, g, device = "tiff",
         width = width_in, height = height_in, units = "in",
         dpi = dpi, compression = "lzw")

  list(plot = g, labels = labs_df)
}

# =============================================================================
# Run
# =============================================================================
t1 <- read_gwas(file_trait1, trait1_name)
t2 <- read_gwas(file_trait2, trait2_name)

res <- make_mirror_manhattan(
  d1 = t1, d2 = t2,
  out_tiff = out_tiff,
  title = paste0("Mirror Manhattan: ", trait1_name, " (top) vs ", trait2_name, " (bottom)"),
  highlight_p = highlight_p,
  ymax = ymax,
  width_in = 12, height_in = 6, dpi = dpi,
  col_dark = col_dark_blue, col_light = col_light_blue, col_sig = col_sig_red
)

# Which SNPs were labelled
print(res$labels %>% select(trait, CHR, BP, SNP, P) %>% arrange(trait, CHR))

cat("Mirror Manhattan written to:", out_tiff, "\n")