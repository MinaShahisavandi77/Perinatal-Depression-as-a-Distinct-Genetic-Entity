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

