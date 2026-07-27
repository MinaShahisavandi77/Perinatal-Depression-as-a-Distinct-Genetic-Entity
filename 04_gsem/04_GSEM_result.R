# =============================================================================
# GSEM model 1 — extract results + effective sample size + component h2
# Author : M. Shahisavandi
#
# From the GWAS-by-subtraction output (stage 04 model script):
#   1. split the MDD_shared and PPD_resid result tables
#   2. write genome-wide-significant SNP subsets
#   3. compute the effective N (Cholesky-adjusted) for each component
#   4. write the final per-component result tables (with Neff)
#   5. univariate LDSC on each derived GWAS -> SNP-h2, saved to the h2 folder
#
#
# =============================================================================

# ---- Packages ---------------------------------------------------------------
repos <- "https://cloud.r-project.org"
for (p in c("dplyr", "data.table", "here")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}
if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = repos)
  devtools::install_github("GenomicSEM/GenomicSEM")
}
library(dplyr)
library(data.table)
library(GenomicSEM)

# ---- Config -----------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R"))
if (!dir.exists(DIR_H2)) dir.create(DIR_H2, recursive = TRUE)

# Output files (all under DIR_GSEM)
f_mdd_shared <- file.path(DIR_GSEM, "MDD_F_shared_results.txt")
f_ppd_resid  <- file.path(DIR_GSEM, "PPD_resid_F_results.txt")
f_sig_mdd    <- file.path(DIR_GSEM, "significant_MDD_F_PPD_shared_SNPs.txt")
f_sig_ppd    <- file.path(DIR_GSEM, "significant_PPD_resid_F_SNPs.txt")

# ---- Analysis constants -----------------------------------------------------
gws_p  <- 5e-8        # genome-wide significance
maf_lo <- 0.10        # MAF window for the effective-N estimate
maf_hi <- 0.40

# Squared path loadings (residual heritabilities) used in the Neff adjustment.
# NB: these are component-specific; update if the preliminary GSEM figures change.
h2_ppd_resid  <- 0.001249188   # = 0.03534386^2
h2_mdd_shared <- 0.09347482    # MDD_shared path loading, squared
#  keeping the value that was actually used. Confirm which is intended.)

# =============================================================================
# 0. Load GWAS-by-subtraction output
# =============================================================================
load(file.path(DIR_GSEM, "outputGWAS.RData"))   # -> outputGWAS

# =============================================================================
# 1. Split the two latent-factor result tables
# =============================================================================
mdd_shared <- subset(outputGWAS[[1]], lhs == "MDD_shared")
ppd_resid  <- subset(outputGWAS[[2]], lhs == "PPD_resid")

head(ppd_resid)

# =============================================================================
# 2. Genome-wide-significant SNPs (before adding Neff)
# =============================================================================
write.table(subset(mdd_shared, Pval_Estimate < gws_p), f_sig_mdd,
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(subset(ppd_resid,  Pval_Estimate < gws_p), f_sig_ppd,
            sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
# 3. Effective sample size (Cholesky-adjusted)
#    Neff per SNP, then mean within a MAF window, floored.
# =============================================================================
compute_neff <- function(df, h2_sq) {
  ((df$Z_Estimate / (df$est * sqrt(h2_sq)))^2) / (2 * df$MAF * (1 - df$MAF))
}

effective_n <- function(df, maf_lo, maf_hi) {
  in_window <- df[df$MAF >= maf_lo & df$MAF <= maf_hi, ]
  floor(mean(in_window$Neff, na.rm = TRUE))
}

ppd_resid$Neff  <- compute_neff(ppd_resid,  h2_ppd_resid)
mdd_shared$Neff <- compute_neff(mdd_shared, h2_mdd_shared)

summary(ppd_resid$Neff)
summary(mdd_shared$Neff)

N_ppd_resid  <- effective_n(ppd_resid,  maf_lo, maf_hi)
N_mdd_shared <- effective_n(mdd_shared, maf_lo, maf_hi)

cat("Effective N  | PPD_resid :", N_ppd_resid,  "\n")  # prior run: 35172
cat("Effective N  | MDD_shared:", N_mdd_shared, "\n")  # prior run: 193060

# =============================================================================
# 4. Write final result tables (with Neff)
# =============================================================================
write.table(mdd_shared, f_mdd_shared, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(ppd_resid,  f_ppd_resid,  sep = "\t", row.names = FALSE, quote = FALSE)

# =============================================================================
# 5. Univariate LDSC on the derived GWAS -> SNP-h2 (observed scale)
#    The shared/residual factors have no population prevalence, so no
#    liability conversion is applied (sample.prev = population.prev = NA).
# =============================================================================

# 5a. Prepare munge-ready files (SNP, A1, A2, signed Z, P, N=Neff)
prep_for_munge <- function(df, path) {
  df %>%
    transmute(SNP = SNP, A1 = A1, A2 = A2,
              Z = Z_Estimate, P = Pval_Estimate, N = Neff) %>%
    fwrite(path, sep = "\t")
  path
}

pre_mdd <- prep_for_munge(mdd_shared, file.path(DIR_GSEM, "MDD_shared.premunge.txt"))
pre_ppd <- prep_for_munge(ppd_resid,  file.path(DIR_GSEM, "PPD_resid.premunge.txt"))

# 5b. munge (writes <trait>.sumstats.gz to the working directory)
setwd(DIR_GSEM)
munge(c(pre_mdd, pre_ppd), HM3_SNPS,
      trait.names = c("MDD_shared", "PPD_resid"),
      N           = c(N_mdd_shared, N_ppd_resid),
      maf.filter  = 0.01)

# 5c. univariate LDSC per component
ldsc_h2 <- function(munged_file, name) {
  out <- ldsc(traits = munged_file,
              sample.prev = NA, population.prev = NA,
              ld = REF_LD_CHR, wld = REF_LD_CHR,
              trait.names = name)
  se <- sqrt(out$V[1, 1])
  tibble(trait = name, scale = "observed",
         h2 = out$S[1, 1], se = se, z = out$S[1, 1] / se)
}

h2_components <- bind_rows(
  ldsc_h2(file.path(DIR_GSEM, "MDD_shared.sumstats.gz"), "MDD_shared"),
  ldsc_h2(file.path(DIR_GSEM, "PPD_resid.sumstats.gz"),  "PPD_resid")
)
print(h2_components)

write.table(h2_components, file.path(DIR_H2, "gsem_components_h2.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Done. Results ->", DIR_GSEM, "| component h2 ->",
    file.path(DIR_H2, "gsem_components_h2.txt"), "\n")