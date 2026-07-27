# =============================================================================
# QC + harmonisation — PPD GWAS
# Author : M. Shahisavandi
# Updated: 2026-02-18
#
# Filters the raw PPD GWAS (MAF, finite beta/SE, heterogeneity, effective N),
# sanity-checks Z vs reported P, harmonises column names for GSEM / LAVA, and
# writes a cleaned sumstats file. PPD betas are already on a linear 0/1 scale.
# =============================================================================

library(dplyr)
library(data.table)

# ---- Paths ------------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R")); 
raw_file <- file.path(GWAS_RAW,"ppd.sumstats.gz" )
outfile  <- file.path(DIR_QC,"ppd.qc.sumstats.gz")

# ---- QC thresholds ----------------------------------------------------------
maf_min      <- 0.01    # drop rare variants
hetisq_max   <- 80      # drop extreme between-cohort heterogeneity
min_n_frac   <- 0.60    # keep SNPs with N >= 60% of max sample size

# =============================================================================
# 1. Load raw GWAS
# =============================================================================
ppd <- fread(raw_file)

# =============================================================================
# 2. Inspect raw distributions
# =============================================================================
summary(ppd$BETA)
summary(ppd$SE)
summary(ppd$Freq1)
summary(ppd$TotalSampleSize)

# Recorded from a prior run (reference):
#   BETA :  min -47.20  median  0.0001  mean -0.0005  max 45.45   <- implausible extremes
#   SE   :  min 0.0126  median  0.0248  mean  0.0997  max 30.27
#   Freq1:  min 0.000   median  0.409   mean  0.460   max 1.000
#   N    :  min 42      median 44334    mean 36202    max 54475   <- some very low N
#
# NOTE: very raw GWAS — biologically impossible min/max betas and some near-zero
#       N. Both the MAF filter and the effective-N filter below are needed.

# =============================================================================
# 3. Compute MAF
# =============================================================================
ppd <- ppd %>%
  mutate(
    MAF = ifelse(Freq1 > 0.5, 1 - Freq1, Freq1)
  )

summary(ppd$MAF)   # MAF looks fine -> proceed to filtering

# =============================================================================
# 4. QC filters
#    - rare variants (MAF outside [0.01, 0.5])
#    - missing / non-finite beta or SE
#    - extreme heterogeneity (HetISq > 80)
#    - small effective sample size (< 60% of max N)
# =============================================================================
clean_ppd <- ppd %>%
  filter(
    is.finite(BETA),
    is.finite(SE),
    MAF >= maf_min & MAF <= 0.5,
    HetISq <= hetisq_max,
    TotalSampleSize >= min_n_frac * max(TotalSampleSize, na.rm = TRUE)
  )

# Optional alternative: also screen on Direction string (disabled)
# clean_ppd <- ppd %>%
#   filter(
#     is.finite(BETA), is.finite(SE),
#     MAF >= maf_min & MAF <= 0.5,
#     HetISq <= hetisq_max,
#     !grepl("\\?", Direction),
#     grepl("^[+\\-]+$", gsub("-", "", Direction))
#   )

# =============================================================================
# 5. Post-filter distributions (extreme betas now gone)
# =============================================================================
summary(clean_ppd$BETA)
summary(clean_ppd$SE)
summary(clean_ppd$Freq1)
summary(clean_ppd$TotalSampleSize)

# Recorded from a prior run (reference):
#   BETA :  min -0.276  median 0.0001  mean 0.0002  max 0.282     <- linear 0/1 scale, no |BETA|>1
#   SE   :  min 0.0126  median 0.0191  mean 0.0246  max 0.100
#   Freq1:  min 0.010   median 0.424   mean 0.464   max 0.990
#   N    :  min 32685   median 45697   mean 46741   max 54475

# =============================================================================
# 6. Z-score vs reported P sanity check
# =============================================================================
clean_ppd <- clean_ppd %>%
  mutate(
    Z      = BETA / SE,
    Pcheck = 2 * (1 - pnorm(abs(Z)))
  )

cor_test <- cor.test(-log10(clean_ppd$P), -log10(clean_ppd$Pcheck))
cat("Correlation between reported P and computed Pcheck:", cor_test$estimate, "\n")
# Prior run: 0.9999919  -> consistent.

# =============================================================================
# 7. Harmonise column names for GSEM / LAVA
#    (CHR / BP / A1 / A2 / P / SNP already match; only N needs renaming)
# =============================================================================
clean_ppd <- clean_ppd %>%
  rename(N = TotalSampleSize) %>%
  select(SNP, CHR, BP, A1, A2, BETA, SE, P, N, MAF)

# =============================================================================
# 8. Final sanity check (expect median ~0, |Z| < ~8 for most SNPs)
# =============================================================================
clean_ppd <- clean_ppd %>%
  mutate(Z = BETA / SE)

summary(clean_ppd$Z)
# Prior run: min -5.20  median 0.005  mean 0.005  max 5.18

# =============================================================================
# 9. Write cleaned sumstats
# =============================================================================
fwrite(clean_ppd, outfile, sep = "\t")
cat("PPD QC completed. File ready for GSEM / MIXER.\n")


