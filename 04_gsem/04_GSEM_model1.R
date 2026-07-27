# =============================================================================
# Genomic SEM — model 1: GWAS-by-subtraction (MDD_shared / PPD_resid)
# Author : M. Shahisavandi
# Updated: 2026-02-19   |   R 4.5.1
#
# Assumes the GWAS were already munged in stage 01 (01_trait_h2.R), so there is
# NO munge step here — this script reads the munged files from DIR_MUNGE.
#
# Steps:
#   1. LDSC covariance structure (liability scale) for the SEM
#      + save per-trait liability-scale heritability to the h2 results
#   2. SEM without SNPs (common-factor / Cholesky)
#   3. Prepare SNP-level sumstats
#   4. GWAS-by-subtraction (userGWAS)
#
# All paths/constants come from config/config.R — edit them there, not here.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
options(stringsAsFactors = FALSE)
repos <- "https://cloud.r-project.org"

cran_pkgs <- c("tidyverse", "data.table", "here")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}
if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = repos)
  devtools::install_github("GenomicSEM/GenomicSEM")
}
library(tidyverse)
library(data.table)
library(GenomicSEM)

# ---- Config -----------------------------------------------------------------
# ---- Config (single source of truth for paths & constants) ------------------
source(here::here("GSEM","scripts", "config.R"))

# Make sure stage output dirs exist

for (d in c(DIR_GSEM, DIR_LDSC_R, DIR_H2)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Munged inputs produced by stage 01 (order: MDD first, PPD second)
munged_mdd_f <- file.path(DIR_MUNGE, "mdd_f.sumstats.gz")
munged_ppd   <- file.path(DIR_MUNGE, "ppd.sumstats.gz")
traits       <- c(munged_mdd_f, munged_ppd)
trait.names  <- c("MDD", "PPD")

# Prevalences in the SAME order as `traits` (MDD, PPD)
samp_prev <- c(SAMP_PREV["MDD"], SAMP_PREV["PPD"])   # c(MDD=0.44,PPD=0.24)
pop_prev  <- c(POP_PREV["MDD"],  POP_PREV["PPD"])    # MDD 0.21 (female), PPD 0.18

stopifnot(all(file.exists(traits)))   # fail early if stage 01 didn't run

# ---- helper: pull per-trait h2 (+ SE, Z) out of an ldsc() output ------------
extract_h2 <- function(ldsc_out, scale_label) {
  S  <- ldsc_out$S
  V  <- ldsc_out$V
  nm <- colnames(S); if (is.null(nm)) nm <- paste0("trait", seq_len(nrow(S)))
  # V is ordered as the column-major lower triangle of S (incl. diagonal)
  rc      <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
  se_all  <- sqrt(diag(V))
  is_diag <- rc[, 1] == rc[, 2]
  tibble(
    trait = nm,
    scale = scale_label,
    h2    = diag(S),
    se    = se_all[is_diag],
    z     = diag(S) / se_all[is_diag]
  )
}

# =============================================================================
# 1. LDSC covariance structure
# ============================================================================= 
# Liability scale (prevalences supplied) — used by the SEM
LDSCoutput <- ldsc(
  traits          = traits,
  sample.prev     = samp_prev,
  population.prev = pop_prev,
  ld              = REF_LD_CHR,
  wld             = REF_LD_CHR,
  trait.names     = trait.names
)
save(LDSCoutput, file = file.path(DIR_LDSC_R, "LDSCoutputMDDPPD.RData"))

# ---- Save liability-scale heritability into the h2 results folder -----------
h2_gsem <- extract_h2(LDSCoutput, "liability")
print(h2_gsem)
write_delim(h2_gsem, file.path(DIR_H2, "gsem_ldsc_h2.txt"), delim = "\t")

# =============================================================================
# 2. SEM without SNPs
# =============================================================================

load(file="LDSCoutputMDDPPD.RData")


model<-'MDD_shared=~NA*PPD + start(0.4)*MDD
        PPD_resid=~NA*PPD
         
         PPD_resid~~1*PPD_resid
         MDD_shared~~1*MDD_shared
         MDD_shared~~0*PPD_resid

         MDD ~~ 0*PPD
         MDD~~0*MDD
         PPD~~0*PPD'


output<-usermodel(LDSCoutput,estimation="DWLS",model=model)

output
save(output, file = file.path(DIR_GSEM, "Modeloutput.Rdata"))


# =============================================================================
# 3. Prepare SNP-level summary statistics
# =============================================================================
p_sumstats <- sumstats(
  files       = c(SUMSTATS_MDD_F, SUMSTATS_PPD),
  ref         = REF_1G,
  trait.names = c("MDD", "PPD"),
  se.logit    = c(TRUE, TRUE),
  info.filter = 0.6,
  maf.filter  = 0.01,
  OLS         = c(FALSE, FALSE),
  linprob     = c(FALSE, FALSE),
  N           = c(N_MDD_F, N_PPD),
  betas       = c("BETA", "BETA")
)
save(p_sumstats, file = file.path(DIR_GSEM, "Sumstats.RData"))

# =============================================================================
# 4. GWAS-by-subtraction
# =============================================================================
model <- '
  MDD_shared =~ NA*PPD + start(0.2)*PPD + start(0.5)*MDD
  PPD_resid  =~ NA*PPD + start(0.2)*PPD
  MDD_shared ~ SNP
  PPD_resid  ~ SNP

  PPD_resid  ~~ 1*PPD_resid
  MDD_shared ~~ 1*MDD_shared
  MDD_shared ~~ 0*PPD_resid

  MDD ~~ 0*PPD
  MDD ~~ 0*MDD
  PPD ~~ 0*PPD
  SNP ~~ SNP
'

outputGWAS <- userGWAS(
  covstruc   = LDSCoutput,
  SNPs       = p_sumstats,
  estimation = "DWLS",
  model      = model,
  sub        = c("MDD_shared~SNP", "PPD_resid~SNP"),
  cor        = 4
)
save(outputGWAS, file = file.path(DIR_GSEM, "outputGWAS.RData"))

head(outputGWAS[[2]][, 1:16])

cat("GSEM model 1 done. h2 -> gsem_ldsc_h2.txt; GWAS-by-subtraction -> outputGWAS.RData\n")