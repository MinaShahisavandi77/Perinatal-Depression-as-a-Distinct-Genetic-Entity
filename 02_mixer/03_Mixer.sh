#!/usr/bin/env bash
#SBATCH --mem=25GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 

set -euo pipefail

# FOR MORE INFORMATION LOOK AT https://github.com/precimed/mixer
###############################################################################
# User inputs (EDIT THESE)
###############################################################################
# Output directory
OUTDIR="/home/m.shahisavandi/GSEM/project/02_MiXeR/"
mkdir -p "$OUTDIR"
SIF="/mixer/gsa-mixer_2.2.1.sif"

SUM1="/GWAS_QC/ppd.qc.sumstats.gz"
SUM2="/GWAS_QC/mdd_f.qc.sumstats.gz"
BIM_FILE="/ld/1000G.EUR.QC.@.bim"
LD_FILE="/ld/1000G.EUR.QC.@.run4.ld"
###############################################################################
# Wrapper
###############################################################################
MIXER="apptainer exec --bind $PWD:/work --pwd /work $SIF python3 /tools/mixer/precimed/mixer.py"

echo "[INFO] MiXeR version:"
$MIXER --version

###############################################################################
# 0) Prepare LD files (run once per LD reference)
# NOTE: The exact flags for 'ld' depend on MiXeR help on your image.
# Check: $MIXER ld --help
###############################################################################
#echo "[INFO] Preparing LD files..."
#$MIXER ld \
 # --bfile "$LD_BFILE" \
  #--out "$OUTDIR/ld"

###############################################################################
# 1) UNIVARIATE: fit1 + test1 for each trait
# NOTE: Exact flags depend on: $MIXER fit1 --help and $MIXER test1 --help
###############################################################################
echo "[INFO] Univariate fit1: PPD"
$MIXER fit1 \
  --trait1-file "$SUM1" \
  --bim-file "$BIM_FILE" \
  --ld-file "$LD_FILE" \
  --chr2use 1-22 \
  --out "$OUTDIR/ppd_fit1" \
  --threads "$SLURM_CPUS_PER_TASK"

echo "[INFO] Univariate test1: PPD"
$MIXER test1 \
  --trait1-file "$SUM1" \
  --bim-file "$BIM_FILE" \
  --ld-file "$LD_FILE" \
  --chr2use 1-22 \
  --load-params-file "$OUTDIR/ppd_fit1.json" \
  --out "$OUTDIR/ppd_test1" \
  --threads "$SLURM_CPUS_PER_TASK"

echo "[INFO] Univariate fit1: MDD"
$MIXER fit1 \
  --trait1-file "$SUM2" \
  --bim-file "$BIM_FILE" \
  --ld-file "$LD_FILE" \
  --chr2use 1-22 \
  --out "$OUTDIR/mdd_fit1" \
  --threads "$SLURM_CPUS_PER_TASK"

#echo "[INFO] Univariate test1: MDD"
$MIXER test1 \
  --trait1-file "$SUM2" \
  --bim-file "$BIM_FILE" \
  --ld-file "$LD_FILE" \
  --chr2use 1-22 \
  --load-params-file "$OUTDIR/mdd_fit1.json" \
  --out "$OUTDIR/mdd_test1" \
  --threads "$SLURM_CPUS_PER_TASK"

###############################################################################
# 2) BIVARIATE: fit2 + test2 (cross-trait)
# NOTE: Exact flags depend on: $MIXER fit2 --help and $MIXER test2 --help
###############################################################################

echo "[INFO] Bivariate fit2: PPD vs MDD"
$MIXER fit2 \
  --trait1-file "$SUM1" \
  --trait2-file "$SUM2" \
  --bim-file "$BIM_FILE" \
  --ld-file "$LD_FILE" \
  --chr2use 1-22 \
  --trait1-params-file "$OUTDIR/mdd_fit1.json"\
  --trait2-params-file "$OUTDIR/ppd_fit1.json"\
  --out "$OUTDIR/ppd_mdd_fit2" \
  --threads "$SLURM_CPUS_PER_TASK"

echo "[INFO] Bivariate test2: PPD vs MDD"
$MIXER test2 \
  --trait1-file "$SUM1" \
  --trait2-file "$SUM2" \
  --bim-file "$BIM_FILE" \
  --ld-file "$LD_FILE" \
  --chr2use 1-22 \
  --load-params-file "$OUTDIR/ppd_mdd_fit2.json" \
  --out "$OUTDIR/ppd_mdd_test2" \
  --threads "$SLURM_CPUS_PER_TASK"
