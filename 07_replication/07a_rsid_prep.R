# =============================================================================
# 07 prep — genotype file for the replication SNP (12:68017195)
# Author : M. Shahisavandi
#
# Reads the raw REF/ALT, dosage (DS) and hardcall (GT) extracts, flips both so
# the effect allele is C (= REF here, A1), and writes the merged file consumed
# by 07_rsid.R.  Paths come from config/config.R.
# =============================================================================

repos <- "https://cloud.r-project.org"
if (!requireNamespace("here", quietly = TRUE)) install.packages("here", repos = repos)

source(here::here("GSEM","scripts", "config.R"))
if (!dir.exists(DIR_SNP_BASE)) dir.create(DIR_SNP_BASE, recursive = TRUE)

# ---- 1. Confirm REF/ALT (expected REF = C, ALT = T -> A1 = C = REF) ----------
ref_alt <- read.table(file.path(DIR_SNP_BASE, "SNP_ref_alt.txt"),
                      col.names = c("CHROM", "POS", "ID", "REF", "ALT"))
print(ref_alt)

# ---- 2. Dosage (DS): flip so dosage counts the C (A1) effect allele ----------
df_ds <- read.table(file.path(DIR_SNP_BASE, "SNP_dosage.txt"),
                    col.names = c("SAMPLE", "DS"))
df_ds$DS_C <- 2 - df_ds$DS

# ---- 3. Hard calls (GT): sum alleles, then flip to A1 = C --------------------
df_gt <- read.table(file.path(DIR_SNP_BASE, "SNP_hardcalls.txt"),
                    col.names = c("SAMPLE", "GT"))
df_gt$hardcall_T <- sapply(df_gt$GT, function(g) {
  sum(as.integer(strsplit(g, "[/|]")[[1]]))   # handles 0/1, 1|0, 1|1, etc.
})
df_gt$hardcall_C <- 2 - df_gt$hardcall_T

# ---- 4. Merge + write --------------------------------------------------------
df <- merge(df_ds, df_gt, by = "SAMPLE")

write.table(df, SNP_DOSAGE_FILE, row.names = FALSE, quote = FALSE, sep = "\t")
head(df)

cat("Wrote effect-allele genotype file ->", SNP_DOSAGE_FILE, "\n")
