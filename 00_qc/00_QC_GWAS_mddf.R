# =============================================================================
# QC + harmonisation — MDD female GWAS
# Author : M. Shahisavandi
# Updated: 2026-02-18
#
# Filters the raw MDD-female GWAS (MAF, finite beta/SE, heterogeneity, effective
# N), sanity-checks Z vs reported P, harmonises column names for GSEM / MiXeR /
# PRS and writes a cleaned sumstats file.
# =============================================================================

library(dplyr)
library(data.table)

# ---- Paths ------------------------------------------------------------------

source(here::here("GSEM","scripts", "config.R")); 
raw_file <- file.path(GWAS_RAW,"MDD_female_AllCohorts.tsv")
outfile  <- file.path(DIR_QC,"mdd_f.qc.sumstats.gz")

# ---- QC thresholds ----------------------------------------------------------
maf_min      <- 0.01    # drop rare variants
hetisq_max   <- 80      # drop extreme between-cohort heterogeneity
min_n_frac   <- 0.60    # keep SNPs with N >= 60% of max effective N

# =============================================================================
# 1. Load raw GWAS
# =============================================================================
mdd_f <- fread(raw_file)

# =============================================================================
# 2. Inspect raw distributions
# =============================================================================
summary(mdd_f$beta)
summary(mdd_f$standard_error)
summary(mdd_f$effect_allele_frequency)
summary(mdd_f$n)

# Recorded from a prior run (reference):
#   beta :  min -0.208  median  0.000  mean  4.3e-06  max  0.176-01
#   SE   :  min 0.0056  median  0.0088 mean  0.0118   max  0.082
#   EAF  :  min 0.010   median  0.425  mean  0.468    max  0.990
#   N    :  min 45079   median 266841  mean 258079    max 288830
#


# =============================================================================
# 3. Compute MAF
# =============================================================================
mdd_f <- mdd_f %>%
  mutate(
    MAF = ifelse(effect_allele_frequency > 0.5,
                 1 - effect_allele_frequency,
                 effect_allele_frequency)
  )

# =============================================================================
# 4. QC filters
#    - rare variants (MAF outside [0.01, 0.5])
#    - missing / non-finite beta or SE
#    - extreme heterogeneity (HetISq > 80)
#    - small effective sample size (< 60% of max N)
# =============================================================================
clean_mdd_f <- mdd_f %>%
  filter(
    is.finite(beta),
    is.finite(standard_error),
    MAF >= maf_min & MAF <= 0.5,
    HetISq <= hetisq_max,
    n >= min_n_frac * max(n, na.rm = TRUE)
  )

# =============================================================================
# 5. Z-score vs reported P sanity check
# =============================================================================
clean_mdd_f <- clean_mdd_f %>%
  mutate(
    Z      = beta / standard_error,
    Pcheck = 2 * (1 - pnorm(abs(Z)))
  )

cor_test <- cor.test(-log10(clean_mdd_f$p_value), -log10(clean_mdd_f$Pcheck))
cat("Correlation between reported P and computed Pcheck:", cor_test$estimate, "\n")
# Prior run: 0.9999582  -> consistent.

# =============================================================================
# 6. Rescale logistic betas to linear 0/1 scale
#    Currently disabled  
# =============================================================================
# P_mdd <- 130471 / (130471 + 159521)   # female case proportion in the meta-analysis
# clean_mdd_f <- clean_mdd_f %>%
#   mutate(
#     BETA_linear = beta            / sqrt(P_mdd * (1 - P_mdd)),
#     SE_linear   = standard_error  / sqrt(P_mdd * (1 - P_mdd))
#   )
# Reference (rescaled beta): min -0.268  median 0.000  mean 1.8e-05  max 0.282

# =============================================================================
# 7. Harmonise column names for GSEM / MiXeR 
# =============================================================================
clean_mdd_f <- clean_mdd_f %>%
  rename(
    CHR  = chromosome,
    BP   = base_pair_location,
    A1   = effect_allele,
    A2   = other_allele,
    P    = p_value,
    N    = n,
    SNP  = rs_id,
    SE   = standard_error,
    BETA = beta
  ) %>%
  select(SNP, CHR, BP, A1, A2, BETA, SE, P, N, MAF)

# =============================================================================
# 8. Final sanity check (expect median ~0, |Z| < ~8-10 for most SNPs)
# =============================================================================
clean_mdd_f <- clean_mdd_f %>%
  mutate(Z = BETA / SE)

summary(clean_mdd_f$Z)
# Prior run: min -6.92  median 0.000  mean 0.0006  max 6.70

# =============================================================================
# 9. Write cleaned sumstats
# =============================================================================
fwrite(clean_mdd_f, outfile, sep = "\t")
cat("MDD female QC completed. File ready for GSEM / MiXeR.\n")