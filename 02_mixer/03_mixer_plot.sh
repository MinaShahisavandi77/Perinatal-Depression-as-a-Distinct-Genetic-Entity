#!/usr/bin/env bash
#SBATCH --mem=25GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

set -euo pipefail

###############################################################################
# User inputs
###############################################################################
OUTDIR="results_final/mix/ppd_mdd_mix"
mkdir -p "$OUTDIR"

SIF="/mixer/gsa-mixer_2.2.1.sif"

# Sumstats
SUM1="/GWAS/ppd.qc.lava.sumstats"
SUM2="/GWAS/mdd_f.qc.lava.sumstats"

# Reference + LD
BIM_FILE="mix/ppd_mdd_mix/ld/1000G.EUR.QC.@.bim"
LD_FILE="$OUTDIR/ld/1000G.EUR.QC.@.run4.ld"

# Where to write figures
FIGDIR="$OUTDIR/figures"
mkdir -p "$FIGDIR"

###############################################################################
# Wrappers
###############################################################################
MIXER="apptainer exec --bind $PWD:/work --pwd /work $SIF python3 /tools/mixer/precimed/mixer.py"
FIGS="apptainer exec --bind $PWD:/work --pwd /work $SIF python3 /tools/mixer/precimed/mixer_figures.py"

echo "[INFO] MiXeR version:"
$MIXER --version

###############################################################################
# 2) VISUALIZATION: univariate figures from fit1 + test1 outputs
# NOTE: mixer_figures.py 'two' needs:
#   --json-fit  (fit1 params json)
#   --json-test (test1 json)
###############################################################################
# Update these if your files are named differently.
PPD_FIT_JSON="$OUTDIR/ppd_fit1.json"
PPD_TEST_JSON="$OUTDIR/ppd_test1.json"

MDD_FIT_JSON="$OUTDIR/mdd_fit1.json"
MDD_TEST_JSON="$OUTDIR/mdd_test1.json"

PPD_MDD_FIT2_JSON="$OUTDIR/ppd_mdd_fit2.json"
PPD_MDD_TEST2_JSON="$OUTDIR/ppd_mdd_test2.json"

# Sanity checks (fail fast with a clear error)
for f in "$PPD_FIT_JSON" "$PPD_TEST_JSON" "$MDD_FIT_JSON" "$MDD_TEST_JSON"; do
  if [[ ! -s "$f" ]]; then
    echo "[ERROR] Missing or empty file: $f"
    echo "        Check your --out prefixes for fit1/test1 and the produced filenames."
    exit 1
  fi
done

echo "[INFO] Plot univariate results: PPD"
$FIGS one \
  --statistic mean std \
  --ext png svg \
  --trait1 PPD \
  --json "$PPD_TEST_JSON" \
  --out "$FIGDIR/ppd"

echo "[INFO] Plot univariate results: MDD"
$FIGS one \
  --statistic mean std \
  --ext png svg \
  --trait1 MDD \
  --json "$MDD_TEST_JSON" \
  --out "$FIGDIR/mdd"

echo "[INFO] Done. Figures written to: $FIGDIR"
echo "[INFO] Generated files:"
ls -lh "$FIGDIR"



echo "[INFO] Plot bivariate results: PPD vs MDD (fit2/test2)"
$FIGS two \
  --statistic point_estimate \
  --ext png svg \
  --trait1 PPD \
  --trait2 MDD \
  --json-fit "$PPD_MDD_FIT2_JSON" \
  --json-test "$PPD_MDD_TEST2_JSON" \
  --out "$FIGDIR/ppd_vs_mdd"