# =============================================================================
# FUMA / MAGMA — GTEx v8 tissue expression enrichment: MDD_shared vs PPD_unique
# Author : M. Shahisavandi
#
# Reads the MAGMA gene-property (tissue) results for both traits, builds a
# side-by-side comparison table, and writes two figures (BETA bars + -log10P
# dot plot). Paths come from config/config.R — edit them there, not here.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
repos <- "https://cloud.r-project.org"
for (p in c("data.table", "dplyr", "ggplot2", "here")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}
library(data.table)
library(dplyr)
library(ggplot2)

# ---- Config -----------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R"))

# ---- Inputs / outputs -------------------------------------------------------
gtex_file <- "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"
mdd_gtex  <- file.path(DIR_FUMA_MDD, gtex_file)
ppd_gtex  <- file.path(DIR_FUMA_PPD, gtex_file)

trait_mdd <- "MDD_shared"
trait_ppd <- "PPD_unique"

outdir <- file.path(DIR_FUMA, "TISSUE_COMPARE")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Plot parameters --------------------------------------------------------
dpi   <- 600
top_n <- 30        # tissues shown, ranked by min(P) across traits

# =============================================================================
# Read a MAGMA gene-property (tissue) result table
# =============================================================================
read_gtex_magma <- function(path, trait) {
  dt <- fread(path, data.table = FALSE, comment.char = "#")
  # expected cols: VARIABLE TYPE NGENES BETA BETA_STD SE P FULL_NAME
  stopifnot(all(c("VARIABLE", "BETA", "SE", "P") %in% names(dt)))

  dt %>%
    transmute(
      trait  = trait,
      tissue = as.character(VARIABLE),
      BETA   = as.numeric(BETA),
      SE     = as.numeric(SE),
      P      = as.numeric(P)
    ) %>%
    filter(is.finite(P), P > 0, P <= 1, is.finite(BETA)) %>%
    mutate(FDR = p.adjust(P, "BH")) %>%
    arrange(P)
}

mdd <- read_gtex_magma(mdd_gtex, trait_mdd)
ppd <- read_gtex_magma(ppd_gtex, trait_ppd)

# =============================================================================
# Comparison table: one row per tissue (side-by-side stats)
# =============================================================================
comp <- full_join(
  mdd %>% select(tissue, BETA_MDD = BETA, SE_MDD = SE, P_MDD = P, FDR_MDD = FDR),
  ppd %>% select(tissue, BETA_PPD = BETA, SE_PPD = SE, P_PPD = P, FDR_PPD = FDR),
  by = "tissue"
) %>%
  mutate(
    minP       = pmin(P_MDD, P_PPD, na.rm = TRUE),
    delta_BETA = BETA_MDD - BETA_PPD
  ) %>%
  arrange(minP)

fwrite(comp, file.path(outdir, "GTEx_v8_tissue_property_compare_MDD_vs_PPD.tsv"),
       sep = "\t", quote = FALSE, na = "NA")

# =============================================================================
# Long format for the top-N tissues (by min P across traits)
# =============================================================================
comp_top <- comp %>%
  arrange(minP) %>%
  slice_head(n = top_n) %>%
  mutate(tissue = factor(tissue, levels = rev(tissue)))

long <- bind_rows(
  comp_top %>% transmute(tissue, trait = trait_mdd, BETA = BETA_MDD, SE = SE_MDD, P = P_MDD, FDR = FDR_MDD),
  comp_top %>% transmute(tissue, trait = trait_ppd, BETA = BETA_PPD, SE = SE_PPD, P = P_PPD, FDR = FDR_PPD)
) %>%
  filter(is.finite(BETA))

# =============================================================================
# Plot A — BETA bars (side-by-side facets)
# =============================================================================
p_beta <- ggplot(long, aes(x = BETA, y = tissue, fill = FDR < 0.05)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
  geom_col(width = 0.72) +
  facet_wrap(~trait, ncol = 2) +
  scale_fill_manual(values = c("TRUE" = "#D7191C", "FALSE" = "#2C7FB8")) +
  theme_classic(base_size = 12) +
  labs(
    title    = paste0("GTEx v8 tissue expression enrichment (MAGMA gene-property): ", trait_mdd, " vs ", trait_ppd),
    subtitle = paste0("Top ", top_n, " tissues by min(P) across traits; bars are BETA; red = FDR<0.05"),
    x = "BETA (positive = enrichment)", y = "Tissue", fill = "FDR<0.05"
  ) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "GTEx_compare_BETA_top30.tiff"),
       p_beta, device = "tiff", dpi = dpi, compression = "lzw",
       width = 14, height = 10, units = "in")

# =============================================================================
# Plot B — -log10(P) dot plot
# =============================================================================
p_p <- ggplot(long, aes(x = -log10(P), y = tissue, color = FDR < 0.05)) +
  geom_point(aes(size = abs(BETA)), alpha = 0.9) +
  facet_wrap(~trait, ncol = 2) +
  scale_color_manual(values = c("TRUE" = "#D7191C", "FALSE" = "grey35")) +
  theme_classic(base_size = 12) +
  labs(
    title    = paste0("GTEx v8 tissue expression enrichment: significance view (", trait_mdd, " vs ", trait_ppd, ")"),
    subtitle = "Point size = |BETA|; red = FDR<0.05",
    x = expression(-log[10](p)), y = "Tissue", color = "FDR<0.05", size = "|BETA|"
  ) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "GTEx_compare_log10P_top30.tiff"),
       p_p, device = "tiff", dpi = dpi, compression = "lzw",
       width = 14, height = 10, units = "in")

# =============================================================================
# Console summary: significant tissues (FDR < 0.05)
# =============================================================================
cat("\nSignificant tissues (FDR<0.05):\n")
cat(trait_mdd, ":", sum(mdd$FDR < 0.05), "\n")
cat(trait_ppd, ":", sum(ppd$FDR < 0.05), "\n")

if (sum(mdd$FDR < 0.05) > 0) {
  cat("\nTop significant tissues for ", trait_mdd, ":\n", sep = "")
  print(mdd %>% filter(FDR < 0.05) %>% select(tissue, BETA, SE, P, FDR) %>% head(20))
}
if (sum(ppd$FDR < 0.05) > 0) {
  cat("\nTop significant tissues for ", trait_ppd, ":\n", sep = "")
  print(ppd %>% filter(FDR < 0.05) %>% select(tissue, BETA, SE, P, FDR) %>% head(20))
}