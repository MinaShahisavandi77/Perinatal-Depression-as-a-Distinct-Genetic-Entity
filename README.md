# Perinatal Depression as a Distinct Genetic Entity

Exploring genetic architecture, biological pathways, and polygenic prediction independent of major depression in females.

## Overview

This repository contains code and scripts used to design and run the analysis pipeline for disentangling the genetic architecture of **perinatal depression (PND)** from **major depression** in females.

## Pipeline summary

The pipeline consists of:

- **LDSC**
  - Univariate analyses
  - Bivariate analyses

- **MiXeR**
  - Univariate analyses
  - Bivariate analyses

- **Genomic Structural Equation Modeling (Genomic SEM)**
  - GWAS-by-subtraction to isolate PND-unique genetic signal and MDD_shared signal

- **Polygenic score (PGS) construction**
  - PGS generation using **SBayesR** (Bayesian regression)

- **PGS evaluation**
  - Testing PGS performance using **linear mixed models (LMMs)**

## Repository organization

Scripts are organized by analysis step (LDSC, MiXeR, Genomic SEM, PGS construction, and PGS testing).

## QC and summary statistics harmonization

Before running LDSC, MiXeR, Genomic SEM (GWAS-by-subtraction), and polygenic score analyses, GWAS summary statistics are QC’ed and harmonized into a consistent format.

For example, `scripts/qc/qc_mdd_female_sumstats.R` performs:
- Loading raw female MDD GWAS summary statistics
- Computing MAF from the effect allele frequency (EAF)
- Filtering variants by:
  - valid BETA/SE
  - MAF threshold (default ≥ 0.01)
  - heterogeneity filter (HetISq ≤ 80 )
  - minimum effective sample size (default ≥ 60% of maximum N)
- Sanity check: recompute Z and compare P-values derived from Z against reported P-values
- Output a harmonized table with columns suitable for downstream tools:
  `SNP, CHR, BP, A1, A2, BETA, SE, P, N, MAF,Z`

> Note: The QC script writes a gzipped TSV when the output path ends with `.gz`.


## LDSC (munge, univariate h2, and bivariate rg)

`scripts/ldsc/01_munge_and_run_ldsc_h2_rg.R`:
- munges QC’d summary statistics into LDSC-ready format (via `GenomicSEM::munge`)
- runs LDSC univariate SNP-heritability (`--h2`) for each trait
- optionally runs LDSC bivariate genetic correlation (`--rg`) for all trait pairs
- extracts results into `h2.tsv` (and `rg.tsv` if enabled)

Paths are provided via environment variables :

```bash
export GWAS_QC_DIR="/path/to/GWAS_QC/"
export MUNG_DIR="/path/to/mung/"
export LDSC_OUT_DIR="/path/to/ldsc_out/"
export LDSC_PY="/path/to/ldsc/ldsc.py"
export REF_LD_DIR="/path/to/eur_w_ld_chr/"
export HM3_SNPLIST="/path/to/eur_w_ld_chr/w_hm3.snplist"

# enable bivariate rg
export RUN_RG=1

Rscript scripts/ldsc/01_munge_and_run_ldsc_h2_rg.R
```
## MiXeR (univariate + bivariate)

This step runs **MiXeR** (via the GSA-MiXeR container) to estimate polygenicity/overlap using:
- **Univariate** models: `fit1` + `test1` for each trait
- **Bivariate** model: `fit2` + `test2` for the trait pair (PPD vs MDD_F)

### Inputs
- QC’d summary statistics for:
  - PPD
  - MDD (female)
- MiXeR LD reference files:
  - BIM pattern file (often uses `@` as the chromosome placeholder)
  - LD pattern file (often uses `@` as the chromosome placeholder)
- Apptainer/Singularity `.sif` image containing MiXeR

### Run (SLURM)

Set environment variables (paths are not hard-coded in the repository):

```bash
export MIXER_OUTDIR="/path/to/results/mix/ppd_mdd_mix"
export MIXER_SIF="/path/to/gsa-mixer_2.2.1.sif"

export SUMSTATS_PPD="/path/to/ppd.qc.sumstats"
export SUMSTATS_MDD_F="/path/to/mdd_f.qc.sumstats"

export MIXER_BIM_FILE="/path/to/ld/1000G.EUR.QC.@.bim"
export MIXER_LD_FILE="/path/to/ld/1000G.EUR.QC.@.run4.ld"

# optional
export CHR2USE="1-22"

sbatch scripts/mixer/run_mixer_ppd_mdd.slurm.sh
```

### Outputs
All MiXeR outputs are written under:

- `${MIXER_OUTDIR}/ppd_fit1.*`
- `${MIXER_OUTDIR}/ppd_test1.*`
- `${MIXER_OUTDIR}/mdd_f_fit1.*`
- `${MIXER_OUTDIR}/mdd_f_test1.*`
- `${MIXER_OUTDIR}/ppd_mdd_f_fit2.*`
- `${MIXER_OUTDIR}/ppd_mdd_f_test2.*`



sbatch scripts/mixer/run_mixer_ppd_mdd.slurm.sh


## Genomic SEM (Model 1: GWAS-by-subtraction)

This step runs **Genomic SEM** to:
1. **Munge** QC’d summary statistics into LDSC-ready format
2. Run **LDSC** to estimate the genetic covariance structure (using sample/population prevalences for liability-scale transforms where applicable)
3. Fit the **SEM model** without SNPs (`usermodel`)
4. Prepare cleaned SNP-level inputs (`sumstats`)
5. Run **GWAS-by-subtraction** (`userGWAS`) to estimate SNP effects on:
   - `MDD_shared`
   - `PPD_resid`

### Run

Set paths via environment variables (no hard-coded local paths are stored in the repository):

```bash
export GSEM_WORKDIR="/path/to/results/model1/"
export HM3_SNPLIST="/path/to/w_hm3.snplist"
export LD_DIR="/path/to/eur_w_ld_chr/"
export WLD_DIR="/path/to/eur_w_ld_chr/"
export REF_MAF="/path/to/reference.1000G.maf.0.005.txt"

export SUMSTATS_PPD_RAW="/path/to/ppd.qc.sumstats.gz"
export SUMSTATS_MDD_RAW="/path/to/mdd.qc.sumstats.gz"

# For liability-scale LDSC transforms (binary traits):
# sample prevalence = cases / (cases + controls) in the GWAS sample
# population prevalence = assumed prevalence in the population
export PPD_SAMPLE_PREV="0.5"
export MDD_SAMPLE_PREV="0.5"
export PPD_POP_PREV="0.18"
export MDD_POP_PREV="0.21"

Rscript scripts/genomicsem/model1_gwas_by_subtraction.R
```

### Outputs
Files are written to `${GSEM_WORKDIR}`, including:
- `LDSCoutput_MDD_PPD.RData`
- `Modeloutput.RData`
- `Sumstats.RData`
- `outputGWAS.RData`

The script also prints a preview of the GWAS results (`head(outputGWAS[[2]][, 1:16])`) at the end.

## Post-processing Genomic SEM GWAS results (Neff + significant SNPs)

`scripts/genomicsem/postprocess_outputGWAS_effective_N.R`:
- splits `outputGWAS.RData` into `MDD_shared` and `PPD_resid`
- exports full result tables
- exports genome-wide significant SNPs (p < 5e-8)
- computes per-SNP effective sample size (Neff) and summarizes Neff after MAF trimming

Run:

```bash
export GSEM_POSTPROC_DIR="/path/to/gsem_model/"
# optional, if the file name differs
# export OUTPUTGWAS_FILE="/path/to/gsem_model/outputGWAS.RData"

# constants used for Neff calculation 
export PPD_H2_COMPONENT="0.001249188"  from SEM model fit 
export MDD_H2_COMPONENT="0.09347482"

Rscript scripts/genomicsem/postprocess_outputGWAS_effective_N.R
```

## SBayesRC: prepare COJO-format summary statistics

SBayesRC (via GCTB) typically expects a COJO-style summary statistics file with columns:

`SNP  A1  A2  freq  b  se  p  N`

`scripts/sbayesrc/make_cojo_sumstats.R` converts:
- Genomic SEM GWAS-by-subtraction outputs (using `est`, `SE`, `Pval_Estimate`, and `Neff` to populate `N`)
- Original GWAS summary statistics (using `BETA`, `SE`, `P`, `N`)

Run:

```bash
export COJO_OUTDIR="/path/to/results/SBayesRC/COJO"

# Genomic SEM result table locations (optional but recommended)
export GSEM_MDD_F_DIR="/path/to/results_model1/MDD_F"
export GSEM_MDD_DIR="/path/to/results_model1/MDD"

# Original GWAS location
export ORIG_GWAS_DIR="/path/to/GWAS"

Rscript scripts/sbayesrc/make_cojo_sumstats.R
```

## PRS weights with SBayesRC (example SLURM script)

`scripts/sbayesrc/run_sbayesrc_prs.slurm.sh` runs the standard SBayesRC workflow:
1. `SBayesRC::tidy()` (harmonize/QC summary statistics to the LD panel)
2. `SBayesRC::impute()` (impute missing SNPs to the LD SNP panel)
3. `SBayesRC::sbayesrc()` (fit SBayesRC and output SNP weights)

### Running the script

Set required inputs as environment variables:

```bash
export SBRC_MA_FILE="/path/to/trait_sumstats_cojo.ma"
export SBRC_LD_DIR="/path/to/ukbEUR_HM3/"
export SBRC_ANNOT_FILE="/path/to/annot_baseline2.2.txt"
export SBRC_OUT_PREFIX="/path/to/results/sbayesrc/trait/trait_eur"

sbatch scripts/sbayesrc/run_sbayesrc_prs.slurm.sh
```

### Tuning vs predefined threshold (study-specific)

- **PPD_unique (PPD_resid)** had **very low SNP heritability**, so we used a **predefined threshold**:
  - `SBRC_MODE=low_h2`
  - `thresh=0.95`
  - `bTune=FALSE`

Example:

```bash
export SBRC_MODE="low_h2"
export SBRC_THRESH="0.95"
sbatch scripts/sbayesrc/run_sbayesrc_prs.slurm.sh
```

- For **other summary statistics / traits**, we enabled **tuning**:
  - `SBRC_MODE=tune`
  - `bTune=TRUE`

Example:

```bash
export SBRC_MODE="tune"
sbatch scripts/sbayesrc/run_sbayesrc_prs.slurm.sh
```

### Output
The main PRS weights file is:

- `${SBRC_OUT_PREFIX}_sbrc.weight`

## PRS evaluation (LMM)

`scripts/pgs/02_gsem_prs_lmm.R`:
- creates a one-pregnancy-per-IID dataset
- exports correlation matrices (Excel)
- reshapes outcomes to long format
- fits LMMs with random intercept `(1|IID)` across multiple PRS/time/history model specifications
- runs a sensitivity analysis stratified by `history_dep`

Run:

```bash
export PGS_INPUT_CSV="/path/to/dataseteuropaen_genr_mothers.csv"
export PGS_OUT_DIR="/path/to/PGS/result"

Rscript scripts/pgs/02_gsem_prs_lmm.R
```

Optional Word export (requires `officer` and `flextable` installed):

```bash
export EXPORT_WORD=1
Rscript scripts/pgs/02_gsem_prs_lmm.R
```



