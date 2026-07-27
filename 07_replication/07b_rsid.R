# =============================================================================
# 07_rsid — ANALYSIS (single-SNP replication, no plots)
# Author : M. Shahisavandi
#
# Merges the prepared genotype file (07_rsid_prep.R) with the one-pregnancy
# cohort, fits Dosage LMMs (full + history-stratified), and writes the result
# tables + the Dosage×time interaction p-value comparison.
# Plot lives in 07_rsid_plot.R.  Paths come from config/config.R.
# =============================================================================

# ---- Packages ---------------------------------------------------------------
repos <- "https://cloud.r-project.org"
pkgs <- c("data.table", "dplyr", "tidyr", "purrr",
          "lme4", "lmerTest", "broom", "broom.mixed", "here")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
library(data.table); library(dplyr); library(tidyr); library(purrr)
library(lme4); library(lmerTest); library(broom); library(broom.mixed)

# ---- Config -----------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R"))
for (d in c(DIR_SNP_DSET, DIR_SNP_RES)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

outcomes <- OUTCOMES
pc_vars  <- PC_VARS
# ---- 1. Load + merge --------------------------------------------------------
dosage <- fread(SNP_DOSAGE_FILE) %>%
  rename(IID = SAMPLE, Dosage = hardcall_C)

mother_single <- read.csv(PGS_ONE, header = TRUE)
mother_dosage <- merge(mother_single, dosage, by = "IID")

# ---- 2. Long-format reshape -------------------------------------------------
mother_long <- mother_dosage %>%
  pivot_longer(cols = all_of(outcomes), names_to = "time", values_to = "dep_value") %>%
  mutate(time = factor(time, levels = outcomes))

write.csv(mother_long, file.path(DIR_SNP_DSET, "mother_long.csv"), row.names = FALSE)

# ---- 3. Model specifications ------------------------------------------------


covariates <- c("h_m_age", "h_m_education", pc_vars)
cov_string <- paste(covariates, collapse = " + ")

plus      <- function(...) paste(c(...), collapse = " + ")
time_term <- "time"
hx_term   <- "history_dep"
int_term  <- "Dosage:time"

model_specs <- list(
  model1 = plus("Dosage"),
  model2 = plus("Dosage", time_term),
  model3 = plus("Dosage", time_term, hx_term),
  model4 = plus("Dosage", time_term, int_term),
  model5 = plus("Dosage", time_term, int_term, hx_term)
)

# ---- 4. Fitting / extraction helpers ----------------------------------------
build_lmm_formula <- function(fx, cov_string) {
  as.formula(paste("dep_value ~", fx, "+", cov_string, "+ (1 | IID)"))
}
fit_lmm_set <- function(data, model_specs, cov_string) {
  formulas <- imap(model_specs, ~ build_lmm_formula(.x, cov_string))
  models   <- map(formulas, ~ lmer(.x, data = data))
  list(formulas = formulas, models = models)
}
extract_lmm_results <- function(models) {
  fixed <- imap_dfr(models, ~ broom.mixed::tidy(.x, effects = "fixed") %>% mutate(model = .y), .id = "model_id")
  fit   <- imap_dfr(models, ~ broom.mixed::glance(.x) %>% mutate(model = .y), .id = "model_id")
  list(fixed = fixed, fit = fit)
}
fit_stratum <- function(dat, model_specs, cov_string) {
  lmm <- fit_lmm_set(dat, model_specs, cov_string)
  res <- extract_lmm_results(lmm$models)
  list(formulas = lmm$formulas, models = lmm$models, fixed = res$fixed, fit = res$fit)
}

# ---- 5. Full-sample (complete-case) models ----------------------------------
lmm_full <- fit_lmm_set(mother_long, model_specs, cov_string)
res_full <- extract_lmm_results(lmm_full$models)

write.csv(res_full$fixed, file.path(DIR_SNP_RES, "lmm_snp.csv"),     row.names = FALSE)
write.csv(res_full$fit,   file.path(DIR_SNP_RES, "lmm_snp_fit.csv"), row.names = FALSE)
saveRDS(lmm_full$models,  file.path(DIR_SNP_RES, "snp_models.rds"))

# ---- 6. Stratified models (history: Yes / No) -------------------------------
mother_long <- mother_long %>% mutate(history_dep = as.factor(history_dep))
print(table(mother_long$history_dep, useNA = "ifany"))

dat_yes <- mother_long %>% filter(history_dep == "Yes")
dat_no  <- mother_long %>% filter(history_dep == "No")

model_specs_stratum <- list(
  model1 = plus("Dosage"),
  model2 = plus("Dosage", time_term),
  model4 = plus("Dosage", time_term, int_term)
)

res_hist_yes <- fit_stratum(dat_yes, model_specs_stratum, cov_string)
res_hist_no  <- fit_stratum(dat_no,  model_specs_stratum, cov_string)

write.csv(res_hist_yes$fit,   file.path(DIR_SNP_RES, "historyYES_lmm_model_fit.csv"),     row.names = FALSE)
write.csv(res_hist_yes$fixed, file.path(DIR_SNP_RES, "historyYES_lmm_fixed_effects.csv"), row.names = FALSE)
write.csv(res_hist_no$fit,    file.path(DIR_SNP_RES, "historyNO_lmm_model_fit.csv"),      row.names = FALSE)
write.csv(res_hist_no$fixed,  file.path(DIR_SNP_RES, "historyNO_lmm_fixed_effects.csv"),  row.names = FALSE)
saveRDS(res_hist_yes$models,  file.path(DIR_SNP_RES, "historyYES_lmm_models.rds"))
saveRDS(res_hist_no$models,   file.path(DIR_SNP_RES, "historyNO_lmm_models.rds"))

# ---- 7. Dosage×time interaction p-values (model 4, three samples) -----------
dosage_interaction_pvals <- function(model, dosage_var = "Dosage", time_levels = outcomes) {
  b <- lme4::fixef(model); V <- as.matrix(vcov(model)); base <- time_levels[1]
  bind_rows(lapply(time_levels, function(t) {
    if (t == base) return(data.frame(time = t, p.value = NA_real_))
    term_int <- paste0(dosage_var, ":time", t)
    if (!term_int %in% names(b)) { warning("Interaction term not found: ", term_int); return(data.frame(time = t, p.value = NA_real_)) }
    est <- unname(b[term_int]); se <- sqrt(V[term_int, term_int]); z <- est / se
    data.frame(time = t, p.value = 2 * pnorm(abs(z), lower.tail = FALSE))
  })) %>%
    mutate(time = factor(time, levels = time_levels),
           sig = case_when(is.na(p.value) ~ "", p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ ""),
           is_sig = !is.na(p.value) & p.value < 0.05)
}

interaction_compare <- bind_rows(
  dosage_interaction_pvals(lmm_full$models[["model4"]])     %>% mutate(stratum = "all mothers"),
  dosage_interaction_pvals(res_hist_no$models[["model4"]])  %>% mutate(stratum = "mothers without history of depression"),
  dosage_interaction_pvals(res_hist_yes$models[["model4"]]) %>% mutate(stratum = "mothers with history of depression")
) %>%
  select(stratum, time, p.value, sig, is_sig) %>%
  arrange(stratum, time)

write.csv(interaction_compare,
          file.path(DIR_SNP_RES, "interaction_model4_pvals_comparison.csv"),
          row.names = FALSE)
print(interaction_compare)

cat("SNP analysis done. Results ->", DIR_SNP_RES, "\n")