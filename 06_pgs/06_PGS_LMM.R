############################################################
# 06_PGS_LMM — ANALYSIS (no plots)
# PRS longitudinal mixed models in Generation R mothers.
# Data prep -> LMM fits (full + history-stratified) -> numeric results.
# Plots live in 06_PGS_plot.R, which reads the outputs of this script.
#
# Paths come from config/config.R — edit them there, not here.
# Date: 2026-04-02
############################################################

# ---- Packages ---------------------------------------------------------------
repos <- "https://cloud.r-project.org"
pkgs <- c("dplyr", "tidyr", "purrr", "openxlsx",
          "lme4", "lmerTest", "broom", "broom.mixed",
          "performance", "stringr", "officer", "flextable", "here")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
library(dplyr); library(tidyr); library(purrr); library(openxlsx)
library(lme4);  library(lmerTest); library(broom); library(broom.mixed)
library(performance); library(stringr); library(officer); library(flextable)

# ---- Config -----------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R"))
for (d in c(DIR_PGS_MODELS, DIR_PGS_LINREG, DIR_PGS_DATASET)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Output files
output_cor_xlsx        <- file.path(DIR_PGS_LINREG,  "correlations.xlsx")
output_mother_long_csv <- file.path(DIR_PGS_DATASET, "dataseteuropaen_genr_mothers_long.csv")
output_mother_csv      <- file.path(DIR_PGS_DATASET, "dataseteuropaen_genr_mothers_one.csv")
out_fixed_effects_csv  <- file.path(DIR_PGS_MODELS,  "lmm_fixed_effects_all_models_with_pvalues.csv")
out_model_fit_csv      <- file.path(DIR_PGS_MODELS,  "lmm_model_fit_all_models.csv")
out_models_rds         <- file.path(DIR_PGS_MODELS,  "lmm_models_list.rds")
out_r2_full_csv        <- file.path(DIR_PGS_MODELS,  "r2_full_models.csv")
out_r2_strat_csv       <- file.path(DIR_PGS_MODELS,  "r2_stratified_models.csv")

# ---- Variables --------------------------------------------------------------
outcomes <- OUTCOMES
pc_vars  <- PC_VARS

pgs_list <- c("MDD_female_org", "MDD_female",
              "PPD_female", "PPD_org")
prs_vars <- c("PPD_org", "MDD_female", "PPD_female", "MDD_female_org")

############################
# 1) Load data
############################
mother <- read.csv(PGS_DATA, header = TRUE)
print(table(mother$MULTIPLE))
print(table(table(mother$IID)))

############################
# 2) Keep one pregnancy per IID
############################
set.seed(123)  # reproducible sampling
mother_onePreg <- mother %>%
  group_by(IID) %>%
  filter(
    if (any(MULTIPLE == "Yes", na.rm = TRUE)) row_number() == sample(row_number(), 1)
    else TRUE
  ) %>%
  ungroup()

print(table(table(mother_onePreg$IID)))   # should all be 1
write.csv(mother_onePreg, output_mother_csv, row.names = FALSE)

############################
# 3) Correlations -> Excel
############################
outcome_cor <- round(cor(mother_onePreg[, outcomes], use = "pairwise.complete.obs"), 3)
pgs_cor     <- round(cor(mother_onePreg[, pgs_list], use = "pairwise.complete.obs"), 3)
cross_cor   <- round(cor(mother_onePreg[, outcomes], mother_onePreg[, pgs_list], use = "pairwise.complete.obs"), 3)

wb <- createWorkbook()
addWorksheet(wb, "Outcome_Correlations"); writeData(wb, "Outcome_Correlations", outcome_cor, rowNames = TRUE)
addWorksheet(wb, "PGS_Correlations");     writeData(wb, "PGS_Correlations",     pgs_cor,     rowNames = TRUE)
addWorksheet(wb, "Outcome_vs_PGS");       writeData(wb, "Outcome_vs_PGS",       cross_cor,   rowNames = TRUE)
saveWorkbook(wb, file = output_cor_xlsx, overwrite = TRUE)

############################
# 4) Long dataset + z-scored PRS
############################
mother_long <- mother_onePreg %>%
  pivot_longer(cols = all_of(outcomes), names_to = "time", values_to = "dep_value") %>%
  mutate(time = factor(time, levels = outcomes))

std_vars <- c(prs_vars)
mother_long <- mother_long %>%
  mutate(across(all_of(std_vars), ~ as.numeric(scale(.x)), .names = "{.col}_z"))

write.csv(mother_long, output_mother_long_csv, row.names = FALSE)

covariates <- c("h_m_age", "h_m_education", pc_vars)
cov_string <- paste(covariates, collapse = " + ")

############################
# 5) Model specs + fitting helpers
############################
plus <- function(...) paste(c(...), collapse = " + ")

time_term <- "time"
hx_term   <- "history_dep"
main_gsem    <- c("MDD_female_z", "PPD_female_z")
main_fem_org <- c("MDD_female_org_z", "PPD_org_z")
int_gsem     <- paste0(main_gsem,    ":time")
int_fem_org  <- paste0(main_fem_org, ":time")

model_specs <- list(
  model1  = plus("MDD_female_z"),
  model2  = plus("PPD_female_z"),
  model3  = plus("MDD_female_org_z"),
  model4  = plus("PPD_org_z"),
  model5  = plus("MDD_female_z", time_term),
  model6  = plus("PPD_female_z", time_term),
  model7  = plus(main_gsem,    time_term),
  model8  = plus(main_fem_org, time_term),
  model9  = plus(main_gsem,    time_term, hx_term),
  model10 = plus(main_fem_org, time_term, hx_term),
  model11 = plus(main_gsem,    time_term, int_gsem),
  model12 = plus(main_fem_org, time_term, int_fem_org),
  model13 = plus(main_gsem,    time_term, int_gsem, hx_term)
)

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
get_r2_df <- function(models) {
  imap_dfr(models, ~{
    r2 <- performance::r2(.x)
    data.frame(model_id = .y,
               R2_marginal    = round(r2$R2_marginal, 3),
               R2_conditional = round(r2$R2_conditional, 3))
  })
}

############################
# 6) Fit LMMs (full sample)
############################
lmm_full <- fit_lmm_set(mother_long, model_specs, cov_string)
res_full <- extract_lmm_results(lmm_full$models)

write.csv(res_full$fixed, out_fixed_effects_csv, row.names = FALSE)
write.csv(res_full$fit,   out_model_fit_csv,     row.names = FALSE)
saveRDS(lmm_full$models,  out_models_rds)

############################
# 7) Sensitivity: history_dep = Yes vs No
############################
mother_long <- mother_long %>% mutate(history_dep = as.factor(history_dep))
print(table(mother_long$history_dep, useNA = "ifany"))

dat_yes <- mother_long %>% filter(history_dep == "Yes")
dat_no  <- mother_long %>% filter(history_dep == "No")

# within-stratum specs (history_dep is constant within stratum -> excluded)
model_specs_stratum <- list(
  model_main_gsem = plus(main_gsem,    time_term),
  model_main_org  = plus(main_fem_org, time_term),
  model_int_gsem  = plus(main_gsem,    time_term, int_gsem),
  model_int_org   = plus(main_fem_org, time_term, int_fem_org)
)

fit_stratum <- function(dat, model_specs, cov_string) {
  lmm <- fit_lmm_set(dat, model_specs, cov_string)
  res <- extract_lmm_results(lmm$models)
  list(formulas = lmm$formulas, models = lmm$models, fixed = res$fixed, fit = res$fit)
}

res_hist_yes <- fit_stratum(dat_yes, model_specs_stratum, cov_string)
res_hist_no  <- fit_stratum(dat_no,  model_specs_stratum, cov_string)

write.csv(res_hist_yes$fit,   file.path(DIR_PGS_MODELS, "historyYES_lmm_model_fit.csv"),     row.names = FALSE)
write.csv(res_hist_yes$fixed, file.path(DIR_PGS_MODELS, "historyYES_lmm_fixed_effects.csv"), row.names = FALSE)
write.csv(res_hist_no$fit,    file.path(DIR_PGS_MODELS, "historyNO_lmm_model_fit.csv"),      row.names = FALSE)
write.csv(res_hist_no$fixed,  file.path(DIR_PGS_MODELS, "historyNO_lmm_fixed_effects.csv"),  row.names = FALSE)
saveRDS(res_hist_yes$models,  file.path(DIR_PGS_MODELS, "historyYES_lmm_models.rds"))
saveRDS(res_hist_no$models,   file.path(DIR_PGS_MODELS, "historyNO_lmm_models.rds"))

############################
# 8) Fit comparison (ΔAIC / ΔBIC within stratum × complexity)
############################
make_fit_comp <- function(fit_df, stratum_label) {
  fit_df %>%
    select(model_id, AIC, BIC) %>%
    mutate(
      stratum = stratum_label,
      type = ifelse(grepl("_gsem$", model_id), "GSEM", "Original"),
      complexity = case_when(
        grepl("^model_main", model_id) ~ "Main effects (+ time)",
        grepl("^model_int",  model_id) ~ "Main effects (+ time) + PRS×time",
        TRUE ~ "Other"
      )
    ) %>%
    group_by(stratum, complexity) %>%
    mutate(dAIC = AIC - min(AIC, na.rm = TRUE),
           dBIC = BIC - min(BIC, na.rm = TRUE)) %>%
    ungroup()
}

fit_comp <- bind_rows(
  make_fit_comp(res_hist_yes$fit, "History: Yes"),
  make_fit_comp(res_hist_no$fit,  "History: No")
) %>% arrange(stratum, complexity, dAIC)

write.csv(fit_comp, file.path(DIR_PGS_MODELS, "sensitivity_history_stratified_fit_comparison.csv"),
          row.names = FALSE)
print(fit_comp)

############################
# 9) R² results  (saved for the record + reused by the plot script)
############################
r2_full <- get_r2_df(lmm_full$models)
write.csv(r2_full, out_r2_full_csv, row.names = FALSE)

r2_strat <- bind_rows(
  get_r2_df(res_hist_yes$models) %>% mutate(stratum = "History: Yes"),
  get_r2_df(res_hist_no$models)  %>% mutate(stratum = "History: No")
)
write.csv(r2_strat, out_r2_strat_csv, row.names = FALSE)

print(r2_full)
print(r2_strat)

############################
# 10) Word results table (Models 1-13, selected terms)
############################
model_levels  <- paste0("model", 1:13)
prs_terms     <- c("MDD_female_z", "PPD_female_z", "MDD_female_org_z", "PPD_org_z")
time_regex    <- "^time"
prs_regex     <- paste0("^(", paste(prs_terms, collapse = "|"), ")$")
int_regex     <- paste0("(^(", paste(prs_terms, collapse = "|"), "):time)|(^time:(", paste(prs_terms, collapse = "|"), "))")
history_regex <- "^history_dep"

keep_term <- function(term) {
  str_detect(term, prs_regex) | str_detect(term, time_regex) |
    str_detect(term, int_regex) | str_detect(term, history_regex)
}
fmt_p <- function(p) case_when(is.na(p) ~ NA_character_,
                               p < 0.001 ~ "<0.001",
                               TRUE ~ formatC(p, format = "f", digits = 3))
term_block <- function(term) case_when(
  str_detect(term, prs_regex) ~ 1L,
  str_detect(term, time_regex) & !str_detect(term, int_regex) ~ 2L,
  str_detect(term, int_regex) ~ 3L,
  str_detect(term, history_regex) ~ 4L,
  TRUE ~ 99L
)

word_tbl_full_df <- res_full$fixed %>%
  filter(effect == "fixed") %>%
  filter(keep_term(term)) %>%
  mutate(model_id = factor(model_id, levels = model_levels),
         Beta = estimate, SE = std.error, p = fmt_p(p.value),
         block = term_block(term), term_sort = term) %>%
  arrange(model_id, block, term_sort) %>%
  select(Model = model_id, Term = term, Beta, SE, p)

doc <- read_docx()
doc <- body_add_par(doc, "LMM selected results (Models 1–13)", style = "heading 1")
doc <- body_add_par(
  doc,
  "Includes: MDD/PPD PRS main effects, time contrasts, PRS×time interaction contrasts, and history_dep. Excludes covariates (age/education/PCs).",
  style = "Normal"
)
ft <- theme_booktabs(autofit(flextable(word_tbl_full_df)))
doc <- body_add_flextable(doc, ft)
print(doc, target = file.path(DIR_PGS_MODELS, "lmm_selected_terms_fullsample_models1to13.docx"))

cat("Analysis done. Models, results, R², and the Word table written to:", DIR_PGS_MODELS, "\n")