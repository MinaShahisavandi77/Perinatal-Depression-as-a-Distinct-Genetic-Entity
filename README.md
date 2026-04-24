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
  - optional heterogeneity filter (HetISq ≤ 80 if available)
  - minimum effective sample size (default ≥ 60% of maximum N)
- Sanity check: recompute Z and compare P-values derived from Z against reported P-values
- Output a harmonized table with columns suitable for downstream tools:
  `SNP, CHR, BP, A1, A2, BETA, SE, P, N, MAF`

> Note: The QC script writes a gzipped TSV when the output path ends with `.gz`.
