#!/usr/bin/env bash
#SBATCH --mem=50GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 


##############################################
# User variables: EDIT THESE
##############################################

# GWAS summary in COJO format (European ancestry)
ma_file="SBayesRC/COJO/mdd_f_orginal_sumstats_cojo.ma"    

# LD reference panel directory (European LD, from SBayesRC Resources)
ld_folder="SBayesRC/ukbEUR_HM3/"   

# Functional annotation file (matching the LD panel, European-based)
annot="SBayesRC/annot_baseline2.2.txt"   

# Output prefix (directory + base name, no extension)
out_prefix="SBayesRC/mdd_f_org/mdd_f_org_eur"          

# Number of CPU cores
threads=4


##############################################
# Usually no change needed below
##############################################

# Export OpenMP threads for underlying C / Fortran code
export OMP_NUM_THREADS=$threads

echo "=================================================="
echo "SBayesRC pipeline for EU GWAS – weights only"
echo "Input MA file : $ma_file"
echo "LD folder     : $ld_folder"
echo "Annotation    : $annot"
echo "Out prefix    : $out_prefix"
echo "Threads       : $threads"
echo "=================================================="

##############################################
# Step 1 – Tidy summary statistics (optional but recommended)
##############################################
echo "[1/3] Tidy GWAS summary statistics..."

Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
                  output='${out_prefix}_tidy.ma', log2file=TRUE)"

echo "Tidy step finished."
echo "Check the log file: ${out_prefix}_tidy.ma.log"


##############################################
# Step 2 – Impute missing SNPs (optional but recommended)
##############################################
echo "[2/3] Imputing summary statistics to SNP panel..."

Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
                  output='${out_prefix}_imp.ma', log2file=TRUE)"

echo "Impute step finished."
echo "Check the log file: ${out_prefix}_imp.ma.log"


##############################################
# Step 3 – Run SBayesRC (main model, produces weights)
##############################################
echo "[3/3] Running SBayesRC with functional annotation..."

Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
                  outPrefix='${out_prefix}_sbrc', annot='$annot', log2file=TRUE)"

echo "SBayesRC finished."
echo "Check the log file: ${out_prefix}_sbrc.log"


##############################################
# Output: WEIGHTS
##############################################
echo "=================================================="
echo "Done."
echo "Your SNP weight file from SBayesRC is:"
echo "    ${out_prefix}_sbrc.weight"
echo "Use this file as the polygenic score weights."
echo "=================================================="


