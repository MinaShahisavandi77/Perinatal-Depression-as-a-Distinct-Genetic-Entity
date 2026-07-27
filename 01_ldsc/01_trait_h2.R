# =============================================================================
# LDSC — univariate heritability (h2) + bivariate genetic correlation (rg)
# Author : M. Shahisavandi
# Updated: 2026-02-18
#
# 1. munge the QC'd GWAS for LDSC
# 2. univariate LDSC  -> per-trait h2            -> h2.txt
# 3. bivariate  LDSC  -> rg + gcov_int per pair  -> rg.txt
# 4. assemble cross-trait intercepts             -> sample_overlap.txt
#
# All paths/constants come from config/config.R — edit them there, not here.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
options(stringsAsFactors = FALSE)
repos <- "https://cloud.r-project.org"

cran_pkgs <- c("tidyverse", "data.table", "optparse", "janitor", "here")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
}
if (!requireNamespace("GenomicSEM", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos = repos)
  devtools::install_github("GenomicSEM/GenomicSEM")
}
library(tidyverse)   # dplyr, stringr, tibble, readr, tidyr, purrr
library(data.table)
library(janitor)
library(GenomicSEM)

# ---- Config (single source of truth for paths & constants) ------------------
source(here::here("GSEM","scripts", "config.R"))

# Make sure stage output dirs exist
for (d in c(DIR_MUNGE, DIR_H2, DIR_RG, DIR_OVERLAP)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# =============================================================================
# 1. Munge summary statistics
#    munge() writes <trait>.sumstats.gz to the working directory.
# =============================================================================
setwd(DIR_MUNGE)

# PPD  — 18 EUR cohorts (17,339 cases / 53,426 controls)
munge(SUMSTATS_PPD, HM3_SNPS,
      trait.names = "ppd",
      N           = N_PPD,
      info.filter = 0.9,
      maf.filter  = 0.01)

# MDD female — gender-stratified (130,471 cases / 159,521 controls)
munge(SUMSTATS_MDD_F, HM3_SNPS,
      trait.names = "mdd_f",
      N           = N_MDD_F,
      info.filter = 0.9,
      maf.filter  = 0.01)

# =============================================================================
# 2. Locate munged files for the phenotypes of interest
# =============================================================================
gwas_details <- tibble(
  full_path = list.files(DIR_MUNGE, pattern = "\\.sumstats\\.gz$",
                         recursive = TRUE, full.names = TRUE)
) %>%
  mutate(file_name = basename(full_path) %>% str_remove("\\.sumstats\\.gz$")) %>%
  filter(str_detect(file_name, str_c(PHENOTYPES, collapse = "|")))

if (nrow(gwas_details) == 0) {
  stop("No munged GWAS matched the phenotype list. Check file names / PHENOTYPES.")
}
print(gwas_details)

# =============================================================================
# 3. UNIVARIATE LDSC: per-trait h2
# =============================================================================
for (i in seq_len(nrow(gwas_details))) {
  h2_arg <- str_c(
    " ", LDSC_PY,
    " --h2 ",         gwas_details$full_path[i],
    " --ref-ld-chr ", REF_LD_CHR,
    " --w-ld-chr ",   REF_LD_CHR,
    " --out ",        file.path(DIR_H2, str_c(gwas_details$file_name[i], "_h2"))
  )
  system2("python", args = h2_arg)
}

# ---- Parse h2 logs ----------------------------------------------------------
parse_h2_log <- function(path) {
  read_lines(path, skip = 25, n_max = 4) %>%
    str_split(": ") %>%
    lapply(function(v) {
      if (str_detect(v[2], "\\(")) {                       # "h2 (se)" style line
        tibble(
          name  = c(v[1], str_c(v[1], " se")),
          value = v[2] %>% str_remove(".*: ") %>% str_split(" ") %>%
                  unlist() %>% parse_number()
        )
      } else {
        tibble(name = v[1], value = parse_number(v[2]))
      }
    }) %>%
    bind_rows() %>%
    mutate(phen = basename(path) %>% str_remove("_h2.log"))
}

h2 <- list.files(DIR_H2, pattern = "_h2.log", full.names = TRUE) %>%
  map(parse_h2_log) %>%
  bind_rows() %>%
  pivot_wider(names_from = name, values_from = value) %>%
  clean_names() %>%
  mutate(z = total_observed_scale_h2 / total_observed_scale_h2_se)

write_delim(h2, file.path(DIR_H2, "h2.txt"), delim = "\t")

# =============================================================================
# 4. BIVARIATE LDSC: genetic correlation (rg) + cross-trait intercept
#    One run per unique unordered pair.
# =============================================================================
pairs <- as.data.frame(t(combn(gwas_details$file_name, 2)),
                       stringsAsFactors = FALSE) %>%
  setNames(c("phen_1", "phen_2")) %>%
  left_join(gwas_details, by = c("phen_1" = "file_name")) %>%
  rename(path_1 = full_path) %>%
  left_join(gwas_details, by = c("phen_2" = "file_name")) %>%
  rename(path_2 = full_path)

for (i in seq_len(nrow(pairs))) {
  rg_arg <- str_c(
    " ", LDSC_PY,
    " --rg ",         pairs$path_1[i], ",", pairs$path_2[i],
    " --ref-ld-chr ", REF_LD_CHR,
    " --w-ld-chr ",   REF_LD_CHR,
    " --out ",        file.path(DIR_RG, str_c(pairs$phen_1[i], "_", pairs$phen_2[i], "_rg"))
  )
  system2("python", args = rg_arg)
}

# ---- Parse rg logs (robust: find the summary table, don't hard-code a skip) -
parse_rg_log <- function(path) {
  lines <- readr::read_lines(path)
  hdr   <- which(stringr::str_detect(lines, "^\\s*p1\\s+p2\\s+rg\\b"))
  if (length(hdr) == 0) { warning("No rg summary table in: ", basename(path)); return(NULL) }
  hdr   <- hdr[1]
  after <- lines[hdr:length(lines)]
  blank <- which(!nzchar(trimws(after)))
  last  <- if (length(blank)) hdr + blank[1] - 2 else length(lines)
  read.table(text = paste(lines[hdr:last], collapse = "\n"),
             header = TRUE, stringsAsFactors = FALSE)
}

rg <- list.files(DIR_RG, pattern = "_rg\\.log$", full.names = TRUE) |>
  purrr::map(parse_rg_log) |>
  dplyr::bind_rows() |>
  dplyr::mutate(
    p1 = basename(p1) |> stringr::str_remove("\\.sumstats\\.gz$"),
    p2 = basename(p2) |> stringr::str_remove("\\.sumstats\\.gz$")
  )

rg   # should now show one row: mdd_f, ppd, rg = 0.9978, ...
write_delim(rg, file.path(DIR_RG, "rg.txt"), delim = "\t")

# =============================================================================
# 5. Sample-overlap matrix from cross-trait intercepts (gcov_int)
#    Off-diagonal = gcov_int; diagonal = 1; then standardise with cov2cor().
# =============================================================================
n <- length(PHENOTYPES)
covar_matrix <- matrix(NA_real_, n, n,
                       dimnames = list(PHENOTYPES, PHENOTYPES))

for (k in seq_len(nrow(rg))) {
  i <- rg$p1[k]; j <- rg$p2[k]
  if (i %in% PHENOTYPES && j %in% PHENOTYPES) {
    covar_matrix[i, j] <- rg$gcov_int[k]
    covar_matrix[j, i] <- rg$gcov_int[k]   # symmetric
  }
}
diag(covar_matrix) <- 1

if (any(is.na(covar_matrix))) {
  stop("sample-overlap matrix has NA entries — a trait pair is missing from rg.txt.")
}

covar_matrix <- round(cov2cor(covar_matrix), 5)

write.table(covar_matrix,
            file = file.path(DIR_OVERLAP, "sample_overlap.txt"),
            quote = FALSE)

cat("LDSC done: h2.txt, rg.txt, sample_overlap.txt written.\n")