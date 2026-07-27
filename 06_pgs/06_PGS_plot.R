############################################################
# 06_PGS_plot — PLOTS ONLY
# Reads the saved outputs of 06_PGS_LMM.R and produces:
#   - full-sample forest plot (Models 1-13) with R²
#   - stratified forest plots (history Yes / No)
#   - interaction-trajectory plots (full + both strata)
#
# Run 06_PGS_LMM.R first. Paths come from config/config.R.
# Date: 2026-04-17
############################################################

# ---- Packages ---------------------------------------------------------------
repos <- "https://cloud.r-project.org"
pkgs <- c("dplyr", "tidyr", "purrr", "stringr",
          "ggplot2", "ggeffects", "performance", "forcats",
          "lme4", "lmerTest", "here")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
library(dplyr); library(tidyr); library(purrr); library(stringr)
library(ggplot2); library(ggeffects); library(performance); library(forcats)
library(lme4); library(lmerTest)

# ---- Config -----------------------------------------------------------------
source(here::here("GSEM","scripts", "config.R"))
if (!dir.exists(DIR_PGS_FIGS)) dir.create(DIR_PGS_FIGS, recursive = TRUE)

# ---- Load analysis outputs --------------------------------------------------
models_full <- readRDS(file.path(DIR_PGS_MODELS, "lmm_models_list.rds"))
fixed_full  <- read.csv(file.path(DIR_PGS_MODELS, "lmm_fixed_effects_all_models_with_pvalues.csv"))

models_yes  <- readRDS(file.path(DIR_PGS_MODELS, "historyYES_lmm_models.rds"))
models_no   <- readRDS(file.path(DIR_PGS_MODELS, "historyNO_lmm_models.rds"))
fixed_yes   <- read.csv(file.path(DIR_PGS_MODELS, "historyYES_lmm_fixed_effects.csv"))
fixed_no    <- read.csv(file.path(DIR_PGS_MODELS, "historyNO_lmm_fixed_effects.csv"))

# ---- Shared constants -------------------------------------------------------
outcomes <- OUTCOMES
time_labels <- c(dep_pre_m = "Pregnancy", dep_2m_m = "2 months",
                 dep_6m_m = "6 months", dep_36m_m = "36 months")
prs_levels <- c(-1, 0, 1)

prs_terms <- c("MDD_female_z", "PPD_female_z", "MDD_female_org_z", "PPD_org_z")
int_regex <- paste0("(^(", paste(prs_terms, collapse = "|"), "):time)|",
                    "(^time:(", paste(prs_terms, collapse = "|"), "))")

my_colors <- c("Shared PRS" = "#0072B2", "Unique PRS" = "#009E73", "Original PRS" = "#D55E00",
               "Shared PRS × Time" = "#56B4E9", "Unique PRS × Time" = "#00A087",
               "Original PRS × Time" = "#CC79A7")

desired_order <- c(
  "MDD shared PRS", "PPD unique PRS", "MDD original PRS", "PPD original PRS",
  "MDD shared × 2 months", "PPD unique × 2 months", "MDD original × 2 months", "PPD original × 2 months",
  "MDD shared × 6 months", "PPD unique × 6 months", "MDD original × 6 months", "PPD original × 6 months",
  "MDD shared × 36 months", "PPD unique × 36 months", "MDD original × 36 months", "PPD original × 36 months"
)

get_r2_df <- function(models) {
  imap_dfr(models, ~{
    r2 <- performance::r2(.x)
    data.frame(model_id = .y,
               R2_marginal = round(r2$R2_marginal, 3),
               R2_conditional = round(r2$R2_conditional, 3))
  })
}

# clean PRS / interaction term labels and colour groups
add_term_labels <- function(df) {
  df %>%
    mutate(
      term_clean = case_when(
        term == "MDD_female_z" ~ "MDD shared PRS",
        term == "PPD_female_z" ~ "PPD unique PRS",
        term == "MDD_female_org_z" ~ "MDD original PRS",
        term == "PPD_org_z" ~ "PPD original PRS",
        str_detect(term, "MDD_female_z:timedep_2m_m")      ~ "MDD shared × 2 months",
        str_detect(term, "PPD_female_z:timedep_2m_m")      ~ "PPD unique × 2 months",
        str_detect(term, "MDD_female_org_z:timedep_2m_m")  ~ "MDD original × 2 months",
        str_detect(term, "PPD_org_z:timedep_2m_m")         ~ "PPD original × 2 months",
        str_detect(term, "MDD_female_z:timedep_6m_m")      ~ "MDD shared × 6 months",
        str_detect(term, "PPD_female_z:timedep_6m_m")      ~ "PPD unique × 6 months",
        str_detect(term, "MDD_female_org_z:timedep_6m_m")  ~ "MDD original × 6 months",
        str_detect(term, "PPD_org_z:timedep_6m_m")         ~ "PPD original × 6 months",
        str_detect(term, "MDD_female_z:timedep_36m_m")     ~ "MDD shared × 36 months",
        str_detect(term, "PPD_female_z:timedep_36m_m")     ~ "PPD unique × 36 months",
        str_detect(term, "MDD_female_org_z:timedep_36m_m") ~ "MDD original × 36 months",
        str_detect(term, "PPD_org_z:timedep_36m_m")        ~ "PPD original × 36 months",
        TRUE ~ term
      ),
      effect_group = case_when(
        str_detect(term_clean, "shared PRS$")   ~ "Shared PRS",
        str_detect(term_clean, "unique PRS$")   ~ "Unique PRS",
        str_detect(term_clean, "original PRS$") ~ "Original PRS",
        str_detect(term_clean, "shared ×")      ~ "Shared PRS × Time",
        str_detect(term_clean, "unique ×")      ~ "Unique PRS × Time",
        str_detect(term_clean, "original ×")    ~ "Original PRS × Time",
        TRUE ~ "Other"
      )
    )
}

# =============================================================================
# 1) FULL-SAMPLE forest plot (Models 1-13)
# =============================================================================
r2_full <- get_r2_df(models_full)

forest_df <- fixed_full %>%
  filter(effect == "fixed") %>%
  filter(term %in% prs_terms | str_detect(term, int_regex)) %>%
  mutate(
    conf.low  = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error,
    stars = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ ""),
    model_id = factor(model_id, levels = paste0("model", 1:13), ordered = TRUE)
  ) %>%
  add_term_labels()

forest_df$term_clean <- factor(forest_df$term_clean, levels = rev(desired_order), ordered = TRUE)

facet_labels <- r2_full %>%
  mutate(model_id = factor(model_id, levels = paste0("model", 1:13), ordered = TRUE),
         facet_lab = paste0(model_id, "\nRm=", R2_marginal, " | Rc=", R2_conditional))

forest_df <- forest_df %>%
  left_join(facet_labels %>% select(model_id, facet_lab), by = "model_id")
forest_df$facet_lab <- factor(forest_df$facet_lab, levels = facet_labels$facet_lab, ordered = TRUE)

p_full <- ggplot(forest_df,
                 aes(x = estimate, y = term_clean, xmin = conf.low, xmax = conf.high, color = effect_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, color = "gray40") +
  geom_pointrange(linewidth = 0.9, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = stars, x = conf.high + 0.03), size = 4, fontface = "bold", color = "black") +
  facet_wrap(~ facet_lab, scales = "free_y", ncol = 3) +
  scale_color_manual(values = my_colors) +
  labs(title = "PRS effects across longitudinal mixed models",
       subtitle = "Interaction effects are relative to baseline (dep_pre_m)",
       x = "Beta estimate (95% CI)", y = "", color = "Effect type") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(face = "bold", size = 17),
        plot.subtitle = element_text(size = 11),
        legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(DIR_PGS_FIGS, "forest_plot_FINAL_ordered.png"), p_full, width = 18, height = 12, dpi = 600)

# =============================================================================
# 2) STRATIFIED forest plots (history Yes / No)
# =============================================================================
make_stratified_forest_plot <- function(model_list, fixed_df, stratum_label, out_dir) {
  r2_df <- get_r2_df(model_list)
  strat_levels <- c("model_main_gsem", "model_main_org", "model_int_gsem", "model_int_org")
  facet_name_map <- c(model_main_gsem = "Main GSEM", model_main_org = "Main Original",
                      model_int_gsem = "Interaction GSEM", model_int_org = "Interaction Original")

  fdf <- fixed_df %>%
    filter(effect == "fixed") %>%
    filter(term %in% prs_terms | str_detect(term, int_regex)) %>%
    mutate(
      conf.low  = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error,
      stars = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ ""),
      model_id = factor(model_id, levels = strat_levels, ordered = TRUE)
    ) %>%
    add_term_labels()
  fdf$term_clean <- factor(fdf$term_clean, levels = rev(desired_order), ordered = TRUE)

  facet_labels <- r2_df %>%
    mutate(model_id = factor(model_id, levels = strat_levels, ordered = TRUE),
           facet_lab = paste0(unname(facet_name_map[as.character(model_id)]),
                              "\nRm=", R2_marginal, " | Rc=", R2_conditional))

  fdf <- fdf %>% left_join(facet_labels %>% select(model_id, facet_lab), by = "model_id")
  fdf$facet_lab <- factor(fdf$facet_lab, levels = facet_labels$facet_lab, ordered = TRUE)

  p <- ggplot(fdf,
              aes(x = estimate, y = term_clean, xmin = conf.low, xmax = conf.high, color = effect_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, color = "gray40") +
    geom_pointrange(linewidth = 0.9, position = position_dodge(width = 0.5)) +
    geom_text(aes(label = stars, x = conf.high + 0.03), size = 4, fontface = "bold", color = "black") +
    facet_wrap(~ facet_lab, scales = "free_y", ncol = 2) +
    scale_color_manual(values = my_colors) +
    labs(title = "PRS effects across stratified longitudinal mixed models",
         subtitle = paste0(stratum_label, " | Interaction effects are relative to baseline (dep_pre_m)"),
         x = "Beta estimate (95% CI)", y = "", color = "Effect type") +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 9, face = "bold"),
          axis.text.x = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(size = 11),
          legend.position = "bottom", panel.grid.minor = element_blank())

  safe_name <- gsub("[^A-Za-z0-9]+", "_", stratum_label)
  ggsave(file.path(out_dir, paste0("forest_plot_stratified_", safe_name, ".png")),
         p, width = 14, height = 10, dpi = 600)
  list(plot = p, data = fdf, r2 = r2_df)
}

forest_hist_yes <- make_stratified_forest_plot(models_yes, fixed_yes, "History: Yes", DIR_PGS_FIGS)
forest_hist_no  <- make_stratified_forest_plot(models_no,  fixed_no,  "History: No",  DIR_PGS_FIGS)

# =============================================================================
# 3) Interaction-trajectory plots (full + strata)
# =============================================================================
term_with_levels <- function(var, levels) paste0(var, " [", paste(levels, collapse = ","), "]")

get_pred_prs_by_time <- function(model, prs_var, prs_label, type_label) {
  ggeffects::ggpredict(model, terms = c("time", term_with_levels(prs_var, prs_levels))) %>%
    as.data.frame() %>%
    mutate(type = type_label, prs = prs_label,
           time = factor(x, levels = outcomes),
           prs_level = as.numeric(as.character(group)))
}

prs_interaction_pvals <- function(model, prs_var, time_levels) {
  b <- lme4::fixef(model); V <- as.matrix(vcov(model)); base <- time_levels[1]
  bind_rows(lapply(time_levels, function(t) {
    if (t == base) return(data.frame(time = t, p.value = NA_real_))
    term_int <- paste0(prs_var, ":time", t)
    if (!term_int %in% names(b)) { warning("Interaction term not found: ", term_int); return(data.frame(time = t, p.value = NA_real_)) }
    est <- unname(b[term_int]); se <- sqrt(V[term_int, term_int]); z <- est / se
    data.frame(time = t, p.value = 2 * pnorm(abs(z), lower.tail = FALSE))
  })) %>%
    mutate(time = factor(time, levels = time_levels),
           sig = case_when(is.na(p.value) ~ "", p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ ""),
           is_sig = !is.na(p.value) & p.value < 0.05)
}

# builds the MDD|PPD trajectory plot from a pair of interaction models
make_interaction_plot <- function(model_gsem, model_org, title, subtitle, out_png, star_size = 7) {
  sig_df <- bind_rows(
    prs_interaction_pvals(model_org,  "MDD_female_org_z", outcomes) %>% mutate(type = "Original", prs = "MDD"),
    prs_interaction_pvals(model_gsem, "MDD_female_z",     outcomes) %>% mutate(type = "GSEM",     prs = "MDD"),
    prs_interaction_pvals(model_org,  "PPD_org_z",        outcomes) %>% mutate(type = "Original", prs = "PPD"),
    prs_interaction_pvals(model_gsem, "PPD_female_z",     outcomes) %>% mutate(type = "GSEM",     prs = "PPD")
  ) %>% select(type, prs, time, sig, is_sig)

  pred_df <- bind_rows(
    get_pred_prs_by_time(model_org,  "MDD_female_org_z", "MDD", "Original"),
    get_pred_prs_by_time(model_gsem, "MDD_female_z",     "MDD", "GSEM"),
    get_pred_prs_by_time(model_org,  "PPD_org_z",        "PPD", "Original"),
    get_pred_prs_by_time(model_gsem, "PPD_female_z",     "PPD", "GSEM")
  ) %>%
    mutate(type = factor(type, levels = c("Original", "GSEM")),
           prs_level_f = factor(prs_level, levels = prs_levels, labels = c("-1 SD", "0", "+1 SD"))) %>%
    left_join(sig_df, by = c("type", "prs", "time")) %>%
    mutate(ribbon_alpha = ifelse(is_sig, 0.22, 0.08)) %>%
    group_by(prs) %>%
    mutate(predicted_z = as.numeric(scale(predicted)),
           conf.low_z  = (conf.low  - mean(predicted, na.rm = TRUE)) / sd(predicted, na.rm = TRUE),
           conf.high_z = (conf.high - mean(predicted, na.rm = TRUE)) / sd(predicted, na.rm = TRUE)) %>%
    ungroup()

  p <- ggplot(pred_df, aes(x = time, y = predicted_z, color = type, linetype = type,
                           group = interaction(type, prs_level_f))) +
    geom_ribbon(aes(ymin = conf.low_z, ymax = conf.high_z, fill = type, alpha = ribbon_alpha), color = NA) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    geom_text(data = pred_df %>% filter(prs_level == 0), aes(label = sig),
              vjust = -0.8, size = star_size, fontface = "bold", show.legend = FALSE) +
    facet_wrap(~ prs, ncol = 2, scales = "free_y") +
    scale_color_manual(values = c(Original = "#1f77b4", GSEM = "#d62728")) +
    scale_fill_manual(values  = c(Original = "#1f77b4", GSEM = "#d62728")) +
    scale_alpha_identity() +
    scale_x_discrete(labels = time_labels) +
    scale_y_continuous(breaks = c(-1, 0, 1)) +
    labs(title = title, subtitle = subtitle, x = "Time",
         y = "standardized predicted values within panel",
         color = NULL, linetype = NULL, fill = NULL) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(out_png, p, width = 11, height = 5.8, dpi = 300)
  p
}

# full sample (interaction models 11 = GSEM, 12 = Original)
make_interaction_plot(
  model_gsem = models_full[["model11"]],
  model_org  = models_full[["model12"]],
  title = "PRS×time interaction models: predicted trajectories with interaction-term significance",
  subtitle = paste0("Ribbon = 95% CI; y-axis standardized (z-score). Stars: PRS×time interaction significant (p<0.05). ",
                    "CIs rescaled consistently with predicted values for visualization."),
  out_png = file.path(DIR_PGS_FIGS, "interaction_trajectories_MDD_PPD_fullsample.png"),
  star_size = 7
)

# stratified
make_interaction_plot(
  model_gsem = models_yes[["model_int_gsem"]],
  model_org  = models_yes[["model_int_org"]],
  title = "PRS×time interaction: predicted trajectories (Psychiatric History: YES)",
  subtitle = "Ribbon = 95% CI; y-axis standardized (z-score). Stars: PRS×time interaction significant (p<0.05).",
  out_png = file.path(DIR_PGS_FIGS, "interaction_trajectories_historyYES.png"),
  star_size = 10
)
make_interaction_plot(
  model_gsem = models_no[["model_int_gsem"]],
  model_org  = models_no[["model_int_org"]],
  title = "PRS×time interaction: predicted trajectories (Psychiatric History: NO)",
  subtitle = "Ribbon = 95% CI; y-axis standardized (z-score). Stars: PRS×time interaction significant (p<0.05).",
  out_png = file.path(DIR_PGS_FIGS, "interaction_trajectories_historyNO.png"),
  star_size = 10
)

cat("All PGS figures written to:", DIR_PGS_FIGS, "\n")