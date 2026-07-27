# =============================================================================
# 07_rsid_plot — PLOT (single-SNP genotype trajectories)
# Author : M. Shahisavandi
#
# Loads the saved model-4 (Dosage×time) objects from 07_rsid.R and draws
# predicted depressive-symptom trajectories by genotype (TT / TC / CC), for
# the full sample and both history strata, side by side.
# Run 07_rsid.R first.  Paths come from config/config.R.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
repos <- "https://cloud.r-project.org"
pkgs <- c("dplyr", "tidyr", "ggplot2", "ggeffects", "lme4", "lmerTest", "here")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
library(dplyr); library(tidyr); library(ggplot2); library(ggeffects)
library(lme4); library(lmerTest)

# ---- Config -----------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R"))
if (!dir.exists(DIR_SNP_FIGS)) dir.create(DIR_SNP_FIGS, recursive = TRUE)

outcomes <- OUTCOMES
time_labels <- c(dep_pre_m = "Pregnancy", dep_2m_m = "2 months",
                 dep_6m_m = "6 months", dep_36m_m = "36 months")

# ---- Load model-4 objects + data --------------------------------------------
m_full <- readRDS(file.path(DIR_SNP_RES, "snp_models.rds"))[["model4"]]
m_yes  <- readRDS(file.path(DIR_SNP_RES, "historyYES_lmm_models.rds"))[["model4"]]
m_no   <- readRDS(file.path(DIR_SNP_RES, "historyNO_lmm_models.rds"))[["model4"]]

mother_long   <- read.csv(file.path(DIR_SNP_DSET, "mother_long.csv"))
dosage_levels <- sort(unique(mother_long$Dosage))

# ---- Genotype labelling (effect allele C: 0=TT, 1=TC, 2=CC) ------------------
geno_map    <- c("0" = "TT", "1" = "TC", "2" = "CC")
geno_levels <- unname(geno_map[as.character(dosage_levels)])
geno_colors <- c(TT = "#2c7fb8", TC = "#7f7f7f", CC = "#d62728")

stratum_levels <- c("All participants",
                    "Mothers with depression history",
                    "Mothers without depression history")

# ---- Predicted trajectories per stratum -------------------------------------
term_with_levels <- function(var, levels) paste0(var, " [", paste(levels, collapse = ","), "]")

get_pred_dosage_by_time <- function(model, stratum_label, dosage_levels) {
  ggpredict(model, terms = c("time", term_with_levels("Dosage", dosage_levels))) %>%
    as.data.frame() %>%
    mutate(stratum   = stratum_label,
           time      = factor(x, levels = outcomes),
           dosage_lv = as.numeric(as.character(group)))
}

pred_df <- bind_rows(
  get_pred_dosage_by_time(m_full, "All participants",                    dosage_levels),
  get_pred_dosage_by_time(m_yes,  "Mothers with depression history",     dosage_levels),
  get_pred_dosage_by_time(m_no,   "Mothers without depression history",  dosage_levels)
) %>%
  mutate(stratum  = factor(stratum, levels = stratum_levels),
         genotype = factor(geno_map[as.character(dosage_lv)], levels = geno_levels))

# ---- Plot (delete geom_ribbon() to drop the confidence band) ----------------
p_compare <- ggplot(pred_df,
                    aes(x = time, y = predicted, color = genotype, fill = genotype, group = genotype)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ stratum, ncol = 3, scales = "free_y") +
  scale_color_manual(values = geno_colors) +
  scale_fill_manual(values  = geno_colors) +
  scale_x_discrete(labels = time_labels) +
  labs(
    title    = "Depressive symptoms during and after pregnancy by 12:68017195 genotype",
    subtitle = "Facets = sample; ribbon = 95% CI. Effect allele C (TT / TC / CC).",
    x = "Time", y = "Predicted depressive symptoms (BSI score)",
    color = "Genotype", fill = "Genotype"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(DIR_SNP_FIGS, "interaction_model4_completecase_vs_stratified.png"),
       p_compare, width = 13, height = 5.2, dpi = 300)

cat("SNP genotype trajectory plot written to:", DIR_SNP_FIGS, "\n")
