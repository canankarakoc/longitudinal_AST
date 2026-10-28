# =============================================================================
# AMR_common.R  -  STABLE LAYER for the P. aeruginosa / CF resistance pipeline
#
# Author: Canan Karakoc
# Assembled: August 25, 2026, by merging the shared machinery
#
# This file holds the things that should NOT need editing to run a new
# analysis: packages, palette and themes, priors, label maps, and every helper
# function. AMR_pipeline.R sources it and does the actual work.
#
#   source("AMR_common.R")
#
# Nothing here reads data or fits a model on load. load_amr_data() is defined
# here but must be called explicitly by the pipeline.
# =============================================================================


# =============================================================================
# PACKAGES AND GLOBAL OPTIONS
# 
# =============================================================================

library(tidyverse)
library(brms)
library(bayesplot)
library(patchwork)
library(tidybayes)
library(ggdist)
library(markovchain)
library(rEDM)
library(readxl)
library(multispatialCCM)   # Clark et al. 2015 multispatial CCM (Section 13)

# NOTE: deliberately NOT using the conflicted package here. Synthesis.R and
# Talk_figures*.R never did, and activating it breaks brms mid-fit: rstan calls
# utils::apropos() while compiling a model, apropos() calls exists() over every
# attached package, and conflicted aborts on ddirichlet (extraDistr vs brms).
# The failure surfaces as "[conflicted] ddirichlet found in 2 packages" from
# inside brm(), nowhere near this line. Load order gives the dplyr verbs we
# want; qualify with :: anywhere it actually matters.

# Defensive: if conflicted is already attached - most likely because
# analytical_data_prep.R was sourced earlier in the SAME R session - detach it.
# Removing it from this file is not enough; once the package is on the search
# path it will abort brm() from inside rstan's compile step no matter which
# script attached it.
if ("package:conflicted" %in% search()) {
  detach("package:conflicted", unload = TRUE, force = TRUE)
  message("AMR_common.R: detached the conflicted package.\n",
          "  It aborts brms mid-fit (rstan -> apropos -> exists -> ddirichlet clash).\n",
          "  If you need it for interactive work, load it AFTER the models are fitted.")
}

options(mc.cores = parallel::detectCores())


# =============================================================================
# OUTPUT DIRECTORIES
# Union of every directory written to by the three source scripts.
# =============================================================================

for (d in c("figures", "results", "models",
            "figures/heterogeneity", "results/heterogeneity", "models/heterogeneity",
            "figures/noise",         "results/noise",         "models/noise",
            "figures/talk",          "figures/ccm",           "results/ccm",
            "figures/perclass",      "results/perclass",      "models/perclass",
            "results/deep")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}


# =============================================================================
# PALETTE AND THEMES
# Extracted from the EEID 2026 poster (Talk_figures_v2.R).
# =============================================================================

# ---- poster palette + theme (extracted from EEID 2026 poster) --------------
pal <- c(
  ink   = "#3D4F52",  # dark slate-teal : text, structural lines, dark category
  black = "black",    # text might look better
  paper = "#F8EAD0",  # cream           : poster panel background
  red   = "#B85C44",  # terracotta      : resistance / positive / "up"
  blue  = "#207E9D",  # teal-blue       : protective / negative / route / converges
  gold  = "#E8C547",  # gold            : drug / exposure / highlight
  plum  = "#90606C",  # muted mauve     : phenotype / carbapenem
  rose  = "#9F334F",  # deep rose       : strong / alternate red
  sand  = "#AB9968"   # khaki           : neutral / baseline / inconclusive
)

# Text/axis colour. Default = ink, for the light/cream slide panels the poster
# uses. If you drop figures straight onto the dark slate background instead,
# set fg <- unname(pal["paper"]).
fg <- unname(pal["black"])

theme_talk <- function(base_size = 16, base_family = "Lato") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      axis.ticks.length = unit(.2, "cm"),
      plot.background   = element_rect(fill = "transparent", colour = NA),
      panel.background  = element_rect(fill = "transparent", colour = NA),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key        = element_rect(fill = "transparent", colour = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey80", linewidth = 0.3),
      text         = element_text(colour = fg, size = 18),
      axis.text    = element_text(colour = fg, size = 18),
      axis.title   = element_text(colour = fg, size = 18),
      legend.text  = element_text(size = 16),
      plot.title   = element_text(colour = fg, face = "bold", hjust = 0,
                                  margin = margin(b = 6), size = 18),
      plot.caption  = element_blank(),
      plot.subtitle = element_blank(),   # subtitles OFF (declutter); text stays in code, hidden
      legend.title  = element_blank(),
      legend.position = "top"
    )
}
# Non-interactive guard: some packages (rEDM) call plot() internally. Under
# Rscript that opens Rplots.pdf in the working directory, which breaks when two
# R processes share a directory and takes every later plot down with it.
# Route implicit plotting to a null device; ggsave() opens its own and is
# unaffected.
if (!interactive()) {
  options(device = function(...) grDevices::pdf(NULL))
}

theme_set(theme_talk())
for (d in c("figures/talk", "figures/heterogeneity", "figures/noise",
            "results/noise", "models", "models/noise", "models/heterogeneity")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Font
library(showtext)
font_add_google("Lato")     # fetches Lato (regular + bold + italic) from Google Fonts
showtext_auto()             # route all text through showtext
showtext_opts(dpi = 300)    # MUST match your ggsave dpi or text sizes come out wrong

# (theme_talk is defined once above, now with base_family = "Lato"; re-set it
#  here so the default theme renders with the registered Lato font.)
theme_set(theme_talk())

# =============================================================================
# PRIORS
# 
# =============================================================================


priors_beta_random <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(gamma(2, 0.1), class = "phi")
)

priors_bernoulli <- c(
  prior(normal(0, 2), class = "b"),
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(1), class = "sd")
)

priors_dist_betabin <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 2),   class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(normal(0, 0.5), class = "b",         dpar = "phi"),
  prior(normal(2, 1),   class = "Intercept", dpar = "phi"),
  prior(exponential(2), class = "sd",        dpar = "phi")
)


# =============================================================================
# LABEL AND COLOUR MAPS
# 
# =============================================================================


# Pretty names: kept short and theme-consistent for clean figure labels
predictor_recodes <- c(
  prev_ag_any_resistant         = "Prior AG resistance",
  prev_has_any_mucoid           = "Prior mucoid phenotype",
  prev_has_mrsa                 = "Prior MRSA",
  prev_has_mssa                 = "Prior MSSA",
  prev_aminoglycoside_prescribed = "Prior AG prescription"
)

# Theme grouping for color/legend (used in some figures)
predictor_theme <- c(
  "Prior AG resistance"     = "Phenotype",
  "Prior mucoid phenotype"  = "Phenotype",
  "Prior MRSA"              = "Species",
  "Prior MSSA"              = "Species",
  "Prior AG prescription"   = "Drug"
)
theme_colors <- c("Phenotype" = "#90606C", "Species" = "#207E9D", "Drug" = "#E8C547")

# =============================================================================
# HELPERS: labels and formatting
# 
# =============================================================================

lab_n <- function(df, patient_col = "individual_id") {
  NULL
  # sprintf("N = %s obs  |  %s patients",
  #         scales::comma(nrow(df)),
  #         dplyr::n_distinct(df[[patient_col]]))
}

cohort_facts <- function(visits = pa_visit_full, isolates = pa_isolates,
                         analysis = pa_analysis, rx = pa_prescription,
                         enrolled = 114L) {   # enrolled subjects: study design / raw file
  has <- function(df, c) c %in% names(df)
  pct <- function(x) sprintf("%.0f%%", 100 * mean(x, na.rm = TRUE))
  
  cat("\n=======================  COHORT FACTS  =======================\n")
  cat("SCALE\n")
  cat(sprintf("  Patients                 : %d\n", n_distinct(visits$individual_id)))
  cat(sprintf("  PA-positive visits       : %s\n", scales::comma(nrow(visits))))
  cat(sprintf("  PA isolates              : %s\n", scales::comma(nrow(isolates))))
  cr <- tryCatch(nrow(readr::read_csv("data/clinical_microbe_data_cleaned.csv",
                                      show_col_types = FALSE)),
                 error = function(e) NA_integer_)
  if (!is.na(cr)) cat(sprintf("  All culture reports      : %s\n", scales::comma(cr)))
  ipc <- isolates %>% count(individual_id, date_collected_from_subject)
  cat(sprintf("  Isolates per culture     : median %g, mean %.1f, max %d\n",
              median(ipc$n), mean(ipc$n), max(ipc$n)))
  
  cat("FOLLOW-UP\n")
  vpp <- visits %>% count(individual_id)
  cat(sprintf("  Visits per patient       : median %g, range %d-%d\n",
              median(vpp$n), min(vpp$n), max(vpp$n)))
  if (has(visits, "days_since_first_visit")) {
    fu <- visits %>% group_by(individual_id) %>%
      summarise(yr = max(days_since_first_visit, na.rm = TRUE) / 365.25, .groups = "drop")
    cat(sprintf("  Follow-up (years)        : median %.1f, max %.1f\n",
                median(fu$yr), max(fu$yr)))
  }
  
  if (has(isolates, "pseudomonas_aeruginosa_phenotype")) {
    cat("PA MORPHOTYPE (phenotype only - no strain typing)\n")
    ph <- isolates %>% count(pseudomonas_aeruginosa_phenotype, sort = TRUE)
    for (i in seq_len(nrow(ph)))
      cat(sprintf("  %-26s %s\n",
                  ifelse(is.na(ph[[1]][i]), "(missing)", ph[[1]][i]),
                  scales::comma(ph$n[i])))
  }
  
  cat("ANTIBIOTICS (PA-relevant), % of isolates with a result\n")
  drug_map <- tibble::tribble(
    ~cls,              ~drug,
    "Aminoglycoside",  "amikacin",
    "Aminoglycoside",  "gentamicin",
    "Aminoglycoside",  "tobramycin",
    "Carbapenem",      "imipenem",
    "Carbapenem",      "meropenem",
    "Cephalosporin",   "cefepime",
    "Cephalosporin",   "ceftazidime",
    "Fluoroquinolone", "ciprofloxacin",
    "Fluoroquinolone", "levofloxacin",
    "Penicillin",      "piperacillin_tazobactam",
    "Monobactam",      "aztreonam",
    "Polymyxin",       "colistin"
  ) %>% filter(drug %in% names(isolates))
  for (i in seq_len(nrow(drug_map)))
    cat(sprintf("  %-16s %-23s %s\n", drug_map$cls[i], drug_map$drug[i],
                pct(!is.na(isolates[[drug_map$drug[i]]]))))
  cat(sprintf("  => %d drugs across %d classes\n",
              nrow(drug_map), n_distinct(drug_map$cls)))
  
  cat("CO-COLONIZERS\n")
  org_col <- intersect(c("co_organisms", "microorganisms"), names(visits))[1]
  if (!is.na(org_col)) {
    sp <- visits %>%
      filter(!is.na(.data[[org_col]]), .data[[org_col]] != "") %>%
      mutate(org = str_split(.data[[org_col]], "[;,/]| and ")) %>%
      tidyr::unnest(org) %>%
      mutate(org = str_trim(org)) %>% filter(org != "") %>%
      count(org, sort = TRUE)
    cat(sprintf("  distinct labels          : %d (from '%s')\n", nrow(sp), org_col))
    for (i in seq_len(min(8, nrow(sp))))
      cat(sprintf("    %-42s %d visits\n", sp$org[i], sp$n[i]))
  }
  if (has(visits, "has_coinfection"))
    cat(sprintf("  visits polymicrobial     : %s\n", pct(visits$has_coinfection == 1)))
  
  cat("RESISTANCE (conservative R + I)\n")
  res_map <- c("Aminoglycoside"  = "ag_any_resistant_conservative",
               "Carbapenem"      = "carb_any_resistant_conservative",
               "Fluoroquinolone" = "fq_any_resistant_conservative",
               "Cephalosporin"   = "ceph_any_resistant_conservative")
  for (nm in names(res_map))
    if (has(visits, res_map[[nm]]))
      cat(sprintf("  %-16s %s\n", nm, pct(visits[[res_map[[nm]]]] == 1)))
  if (has(isolates, "multi_drug_resistant"))
    cat(sprintf("  MDR-PA (isolate-level)   : %s\n", pct(isolates$multi_drug_resistant == 1)))
  
  cat("DATA FUNNEL (raw -> what the models use)\n")
  cat(sprintf("  Enrolled subjects        : %d        [study design / raw file]\n", enrolled))
  if (!is.na(cr))
    cat(sprintf("  Cleaned isolate records  : %s     [prep drops ~20 empty/dateless rows]\n",
                scales::comma(cr)))
  cat(sprintf("  -> PA isolates           : %s     [rest are co-colonizers = ecology context]\n",
              scales::comma(nrow(isolates))))
  cat(sprintf("  PA-positive patients     : %d\n", n_distinct(visits$individual_id)))
  cat(sprintf("  PA visits                : %s\n", scales::comma(nrow(visits))))
  cat(sprintf("  -> analysis subset       : %s (%.0f%% of visits), %d pts  [drop first visits: no prior to lag]\n",
              scales::comma(nrow(analysis)), 100 * nrow(analysis) / nrow(visits),
              n_distinct(analysis$individual_id)))
  cat(sprintf("  -> drug/prescription set : %s (%.0f%% of visits), %d pts  [exacerbation visits only]\n",
              scales::comma(nrow(rx)), 100 * nrow(rx) / nrow(visits),
              n_distinct(rx$individual_id)))
  cat("  (raw record count is printed by the prep script: 'Raw data loaded: N rows')\n")
  cat("==============================================================\n\n")
  invisible(NULL)
}

pct_present <- function(df, cols) {
  cols <- cols[cols %in% names(df)]
  if (!length(cols)) return(tibble(item = character(), pct = numeric()))
  map_dfr(cols, ~tibble(item = .x, pct = mean(df[[.x]] == 1, na.rm = TRUE)))
}

pretty_term <- function(x) {
  x %>% str_remove("^b_") %>% str_remove("^prev_") %>%
    str_remove("_prescribed$") %>% str_remove("_exposure_z$") %>%
    str_remove("_exposure$") %>%
    recode(
      "any_resistance"  = "Prior resistance",
      "has_mrsa"        = "MRSA",
      "has_mssa"        = "MSSA",
      "has_any_mucoid"  = "Mucoid",
      "aminoglycoside"  = "Aminoglycoside",
      "carbapenem"      = "Carbapenem",
      "cephalosporin"   = "Cephalosporin",
      "fluoroquinolone" = "Fluoroquinolone",
      "fluoroquinolones"= "Fluoroquinolone",
      "colistin"        = "Colistin",
      "total_pressure_z"= "Total drug pressure",
      "has_exacerbation"= "Exacerbation"
    )
}

relabel_pred <- function(var, persist_col) {
  m <- setNames(c("Prior same-class R", common_labels), c(persist_col, common_vars))
  unname(m[var])
}

sev_recode <- function(x) factor(dplyr::case_when(
  x == "Clinically Stable"  ~ "Stable",
  x == "APE Outpatient"     ~ "APEout",
  x == "APE IV Antibiotics" ~ "APEiv"),
  levels = c("Stable", "APEout", "APEiv"))

fit_class_sev <- function(cfg_row) {
  nres <- cfg_row$nres; ntest <- cfg_row$ntest; persist <- cfg_row$persist
  d <- pa_visit_full %>%
    filter(.data[[ntest]] >= 2, !is.na(.data[[persist]]), !is.na(clinical_status),
           if_all(all_of(common_vars), ~ !is.na(.))) %>%
    mutate(clinical_status = sev_recode(clinical_status))
  cat(sprintf("\n[%-15s +status] N = %d visits (Stable/APEout/APEiv = %s)\n",
              cfg_row$label, nrow(d), paste(table(d$clinical_status), collapse = "/")))
  mu_str  <- sprintf("%s | trials(%s) ~ %s + %s + clinical_status + (1 | individual_id)",
                     nres, ntest, persist, common_rhs)
  phi_str <- sprintf("phi ~ %s + %s + (1 | individual_id)", persist, common_rhs)
  m <- brm(bf(as.formula(mu_str), as.formula(phi_str), family = beta_binomial()),
           data = d, prior = priors_dist_betabin,
           chains = 4, cores = 4, iter = 4000, warmup = 2000,
           control = list(adapt_delta = 0.95, max_treedepth = 12), seed = 42,
           file = paste0("models/perclass/m_dist_sev_", cfg_row$cls), file_refit = "on_change")
  list(cls = cfg_row$cls, label = cfg_row$label, persist = persist, m = m)
}


# =============================================================================
# HELPERS: data loading
# 
# =============================================================================

# Prefer the .rds files: the csv round-trip drops every factor (clinical_status
# and pa_phenotype come back as character with their level order lost).
# Falls back to csv with a warning if the rds are missing.
load_amr_data <- function(dir = "data", quiet = FALSE) {
  rds <- file.path(dir, c("pa_visit_full.rds", "pa_isolates.rds"))
  if (all(file.exists(rds))) {
    out <- list(pa_visit_full = readRDS(rds[1]), pa_isolates = readRDS(rds[2]))
  } else {
    warning("rds not found in ", dir, "; falling back to csv - factor levels will be lost.")
    out <- list(
      pa_visit_full = readr::read_csv(file.path(dir, "pa_visit_full.csv"), show_col_types = FALSE),
      pa_isolates   = readr::read_csv(file.path(dir, "pa_isolates.csv"),   show_col_types = FALSE)
    )
  }
  if (!quiet) {
    cat("Loaded pa_visit_full:", nrow(out$pa_visit_full), "visit-rows\n")
    cat("Loaded pa_isolates  :", nrow(out$pa_isolates),   "isolate-rows\n\n")
  }
  out
}


# =============================================================================
# HELPERS: model fitting
# 
# =============================================================================

fit_stratified_models <- function(data, class_name, 
                                  outcome_var, prev_resist_var, prev_rx_var,
                                  file_prefix) {
  
  df <- data %>%
    filter(!is.na(.data[[prev_resist_var]]), !is.na(.data[[outcome_var]]))
  
  # Predictors that actually enter each formula. brms drops rows with NA in
  # any of them, so the sample-size gate below must be applied to the
  # complete-case frame -- otherwise a stratum clears "N >= 50" and is then
  # fitted on far fewer rows (this is how M5A cleared the gate and fitted 33).
  preds_acq <- c("prev_any_resistance", prev_rx_var,
                 "prev_has_mrsa", "prev_has_mssa", "prev_has_any_mucoid")
  preds_per <- c(prev_rx_var,
                 "prev_has_mrsa", "prev_has_mssa", "prev_has_any_mucoid")
  
  df_susceptible_raw <- df %>% filter(.data[[prev_resist_var]] == 0)
  df_resistant_raw   <- df %>% filter(.data[[prev_resist_var]] == 1)
  
  df_susceptible <- df_susceptible_raw %>%
    filter(if_all(all_of(preds_acq), ~ !is.na(.)))
  df_resistant <- df_resistant_raw %>%
    filter(if_all(all_of(preds_per), ~ !is.na(.)))
  
  cat("\n", class_name, ":\n")
  cat("  Dropped for missing predictors: acq",
      nrow(df_susceptible_raw) - nrow(df_susceptible), "of", nrow(df_susceptible_raw),
      "| per", nrow(df_resistant_raw) - nrow(df_resistant), "of", nrow(df_resistant_raw), "\n")
  cat("  Susceptible → ? : N =", nrow(df_susceptible), 
      ", Events =", sum(df_susceptible[[outcome_var]]), "\n")
  cat("  Resistant → ? : N =", nrow(df_resistant), 
      ", Events =", sum(df_resistant[[outcome_var]]), "\n")
  
  # Acquisition model
  m_acquisition <- NULL
  if (nrow(df_susceptible) >= 50 & sum(df_susceptible[[outcome_var]]) >= 10) {
    cat("  Fitting acquisition model...\n")
    
    formula_acq <- as.formula(paste0(
      outcome_var, " ~ ",
      "prev_any_resistance + ",
      prev_rx_var, " + ",
      "prev_has_mrsa + prev_has_mssa + ",
      "prev_has_any_mucoid + ",
      # prev_has_exacerbation omitted: constant in pa_prescription (APE-only).
      "(1 | individual_id)"
    ))
    
    m_acquisition <- brm(
      formula_acq,
      data = df_susceptible,
      family = bernoulli(),
      prior = priors_bernoulli,
      chains = 4, cores = 4,
      iter = 2000, warmup = 1000,
      # brms defaults to adapt_delta = 0.8; every other fit in this pipeline
      # uses >= 0.95, and 0.8 produced divergent transitions in these strata.
      control = list(adapt_delta = 0.95),
      seed = 42,
      file = paste0("models/", file_prefix, "_acquisition"),
      file_refit = "on_change"
    )
  }
  
  # Persistence model
  m_persistence <- NULL
  if (nrow(df_resistant) >= 50 & sum(df_resistant[[outcome_var]]) >= 10) {
    cat("  Fitting persistence model...\n")
    
    formula_per <- as.formula(paste0(
      outcome_var, " ~ ",
      prev_rx_var, " + ",
      "prev_has_mrsa + prev_has_mssa + ",
      "prev_has_any_mucoid + ",
      # prev_has_exacerbation omitted: constant in pa_prescription (APE-only).
      "(1 | individual_id)"
    ))
    
    m_persistence <- brm(
      formula_per,
      data = df_resistant,
      family = bernoulli(),
      prior = priors_bernoulli,
      chains = 4, cores = 4,
      iter = 2000, warmup = 1000,
      control = list(adapt_delta = 0.95),
      seed = 42,
      file = paste0("models/", file_prefix, "_persistence"),
      file_refit = "on_change"
    )
  }
  
  return(list(
    acquisition = m_acquisition,
    persistence = m_persistence,
    data_susceptible = df_susceptible,
    data_resistant = df_resistant
  ))
}

fit_cross_model <- function(outcome_var, prev_outcome_var, 
                            cross_rx_var, same_rx_var,
                            model_name, file_name) {
  
  formula_str <- paste0(
    outcome_var, " ~ ",
    prev_outcome_var, " + ",
    cross_rx_var, " + ",
    same_rx_var, " + ",
    "prev_has_mrsa + prev_has_mssa + ",
    "prev_has_any_mucoid + ",
    # prev_has_exacerbation omitted: constant in pa_prescription (APE-only).
    "(1 | individual_id)"
  )
  
  cat("Fitting", model_name, "...\n")
  
  brm(
    as.formula(formula_str),
    data = pa_cross,
    family = bernoulli(),
    prior = priors_bernoulli,
    chains = 4, cores = 4,
    iter = 2000, warmup = 1000,
    control = list(adapt_delta = 0.95),
    seed = 42,
    file = paste0("models/", file_name),
    file_refit = "on_change"
  )
}

# Mean pairwise correlation between the AG drug calls made on the SAME colony.
# Used twice: to deflate the drug-test denominator in the measurement-error
# model (Section 4), and to size the design effect in the depth-reliability
# curve (Section 14). Defined once here so the two cannot drift apart.
ag_within_colony_rho <- function(iso_df,
                                 drugs = c("amikacin", "tobramycin", "gentamicin"),
                                 nonsusceptible = c("Resistant", "Intermediate")) {
  dd <- iso_df %>%
    filter(if_all(all_of(drugs), ~ .x %in% c("Susceptible", "Intermediate", "Resistant"))) %>%
    transmute(across(all_of(drugs), ~ .x %in% nonsusceptible))
  if (nrow(dd) < 10) return(NA_real_)
  cmb <- utils::combn(drugs, 2)
  mean(apply(cmb, 2, function(p) suppressWarnings(cor(dd[[p[1]]], dd[[p[2]]]))),
       na.rm = TRUE)
}

fit_variance_decomp <- function(iso_df, drug_cols, class_label, file_tag) {
  d <- iso_df %>%
    rowwise() %>%
    mutate(
      .vals_tested = sum(c_across(all_of(drug_cols)) %in%
                           c("Resistant", "Intermediate", "Susceptible")),
      R = if_else(
        .vals_tested == 0, NA_integer_,
        as.integer(any(c_across(all_of(drug_cols)) %in%
                         c("Resistant", "Intermediate")))   # CONSERVATIVE
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(R))
  
  m <- brm(
    R ~ 1 + (1 | individual_id) + (1 | individual_id:visit_number),
    data = d, family = bernoulli(),
    prior = c(prior(normal(0, 2), class = "Intercept"),
              prior(exponential(1), class = "sd")),
    chains = 4, cores = 4, iter = 4000, warmup = 2000,
    control = list(adapt_delta = 0.95), seed = 42,
    file = paste0("models/noise/vardecomp_", file_tag),
    file_refit = "on_change"
  )
  
  dr      <- as_draws_df(m)
  v_pat   <- dr[["sd_individual_id__Intercept"]]^2
  v_vis   <- dr[["sd_individual_id:visit_number__Intercept"]]^2
  v_resid <- pi^2 / 3
  tot     <- v_pat + v_vis + v_resid
  
  mk <- function(v, comp) tibble(
    class = class_label, component = comp,
    median = median(v),
    lo = as.numeric(quantile(v, .025)),
    hi = as.numeric(quantile(v, .975))
  )
  bind_rows(
    mk(v_pat   / tot, "between_patient"),
    mk(v_vis   / tot, "between_visit_within_patient"),
    mk(v_resid / tot, "between_isolate_within_visit")
  )
}

fit_staph_stratum <- function(data, staph_value, label) {
  d <- data %>% filter(has_any_staph == staph_value)
  if (nrow(d) < 80) {
    cat("Skipping", label, "- only", nrow(d), "visits\n")
    return(NULL)
  }
  cat("\nFitting Staph stratum:", label,
      "- N =", nrow(d), "visits,", n_distinct(d$individual_id), "patients\n")
  
  # Within stratum, has_any_staph is constant, so we drop it from formula.
  # MRSA and MSSA also drop out of the Staph- stratum (always 0). To keep
  # the same model shape across strata, we use only the predictors that
  # vary in BOTH strata: prev_ag_any_resistant, prev_has_any_mucoid,
  # prev_aminoglycoside_prescribed.
  brm(
    bf(
      ag_n_resistant | trials(ag_n_tested) ~
        prev_ag_any_resistant +
        prev_has_any_mucoid +
        prev_aminoglycoside_prescribed +
        (1 | individual_id),
      phi ~
        prev_ag_any_resistant +
        prev_has_any_mucoid +
        prev_aminoglycoside_prescribed +
        (1 | individual_id),
      family = beta_binomial()
    ),
    data = d,
    prior = priors_dist_betabin,
    chains = 4, cores = 4,
    iter = 4000, warmup = 2000,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    seed = 42,
    file = paste0("models/heterogeneity/m_staph_v5_", label),
    file_refit = "on_change"
  )
}

fit_one_class <- function(cfg_row) {
  nres <- cfg_row$nres; ntest <- cfg_row$ntest; persist <- cfg_row$persist
  d <- pa_visit_full %>%
    filter(.data[[ntest]] >= 2, !is.na(.data[[persist]]),
           if_all(all_of(common_vars), ~ !is.na(.)))
  # NB: ">= 2 drug-tests", not ">= 2 isolates". For a two-agent class one colony
  # tested on both agents satisfies it; see Section 10.2.
  cat(sprintf("\n[%-15s] N = %d visits, %d patients (>=2 %s DRUG-TESTS, not colonies)\n",
              cfg_row$label, nrow(d), n_distinct(d$individual_id), cfg_row$cls))
  
  mu_str  <- sprintf("%s | trials(%s) ~ %s + %s + (1 | individual_id)",
                     nres, ntest, persist, common_rhs)
  phi_str <- sprintf("phi ~ %s + %s + (1 | individual_id)", persist, common_rhs)
  
  m_naive <- brm(bf(as.formula(mu_str), family = beta_binomial()),
                 data = d, prior = priors_beta_random,
                 chains = 4, cores = 4, iter = 4000, warmup = 2000,
                 control = list(adapt_delta = 0.95), seed = 42,
                 file = paste0("models/perclass/m_naive_", cfg_row$cls), file_refit = "on_change")
  
  m_dist <- brm(bf(as.formula(mu_str), as.formula(phi_str), family = beta_binomial()),
                data = d, prior = priors_dist_betabin,
                chains = 4, cores = 4, iter = 4000, warmup = 2000,
                control = list(adapt_delta = 0.95, max_treedepth = 12), seed = 42,
                file = paste0("models/perclass/m_dist_", cfg_row$cls), file_refit = "on_change")
  
  list(cls = cfg_row$cls, label = cfg_row$label, persist = persist,
       n = nrow(d), np = n_distinct(d$individual_id),
       m_naive = m_naive, m_dist = m_dist)
}

fit_class_sev <- function(cfg_row) {
  nres <- cfg_row$nres; ntest <- cfg_row$ntest; persist <- cfg_row$persist
  d <- pa_visit_full %>%
    filter(.data[[ntest]] >= 2, !is.na(.data[[persist]]), !is.na(clinical_status),
           if_all(all_of(common_vars), ~ !is.na(.))) %>%
    mutate(clinical_status = sev_recode(clinical_status))
  cat(sprintf("\n[%-15s +status] N = %d visits (Stable/APEout/APEiv = %s)\n",
              cfg_row$label, nrow(d), paste(table(d$clinical_status), collapse = "/")))
  mu_str  <- sprintf("%s | trials(%s) ~ %s + %s + clinical_status + (1 | individual_id)",
                     nres, ntest, persist, common_rhs)
  phi_str <- sprintf("phi ~ %s + %s + (1 | individual_id)", persist, common_rhs)
  m <- brm(bf(as.formula(mu_str), as.formula(phi_str), family = beta_binomial()),
           data = d, prior = priors_dist_betabin,
           chains = 4, cores = 4, iter = 4000, warmup = 2000,
           control = list(adapt_delta = 0.95, max_treedepth = 12), seed = 42,
           file = paste0("models/perclass/m_dist_sev_", cfg_row$cls), file_refit = "on_change")
  list(cls = cfg_row$cls, label = cfg_row$label, persist = persist, m = m)
}


# =============================================================================
# HELPERS: posterior extraction
# 
# =============================================================================

extract_effects <- function(model, model_name, params_map) {
  if (is.null(model)) return(NULL)
  
  posterior_df <- as_draws_df(model)
  
  effects <- map_dfr(names(params_map), function(param) {
    if (param %in% names(posterior_df)) {
      draws <- posterior_df[[param]]
      tibble(
        Model = model_name,
        Variable = params_map[[param]],
        Parameter = param,
        Estimate = mean(draws, na.rm = TRUE),
        SD = sd(draws, na.rm = TRUE),
        CI_lower = quantile(draws, 0.025, na.rm = TRUE),
        CI_upper = quantile(draws, 0.975, na.rm = TRUE),
        Prob_positive = mean(draws > 0, na.rm = TRUE),
        OR = exp(mean(draws, na.rm = TRUE)),
        OR_lower = exp(quantile(draws, 0.025, na.rm = TRUE)),
        OR_upper = exp(quantile(draws, 0.975, na.rm = TRUE))
      )
    }
  })
  
  if (is.null(effects) || nrow(effects) == 0) return(NULL)
  
  effects %>%
    mutate(
      # Labelled as probability of direction, NOT as stars. A clinical reader
      # sees "*" and reads p < 0.05; the lowest tier here was pd > 0.90, which
      # is a one-sided p of about 0.10. Reporting pd explicitly removes the
      # false equivalence, and matches the pd > 0.975 cut used elsewhere in this
      # file. Significance_stars is kept for figures that still expect a glyph,
      # but its lowest tier now starts at 0.95.
      Significance = case_when(
        pmax(Prob_positive, 1 - Prob_positive) > 0.99  ~ "pd>0.99",
        pmax(Prob_positive, 1 - Prob_positive) > 0.975 ~ "pd>0.975",
        pmax(Prob_positive, 1 - Prob_positive) > 0.95  ~ "pd>0.95",
        pmax(Prob_positive, 1 - Prob_positive) > 0.90  ~ "pd>0.90",
        TRUE ~ ""
      ),
      Significance_stars = case_when(
        pmax(Prob_positive, 1 - Prob_positive) > 0.99  ~ "***",
        pmax(Prob_positive, 1 - Prob_positive) > 0.975 ~ "**",
        pmax(Prob_positive, 1 - Prob_positive) > 0.95  ~ "*",
        TRUE ~ ""
      ),
      Direction = case_when(
        Prob_positive > 0.95 ~ "↑ Resistance",
        Prob_positive < 0.05 ~ "↓ Resistance",
        TRUE ~ "Uncertain"
      )
    )
}

extract_diagnostics <- function(model, model_name) {
  if (is.null(model)) return(NULL)
  
  
  # Get summary
  model_summary <- summary(model)
  
  
  # Extract fixed effects diagnostics
  fixed_df <- as.data.frame(model_summary$fixed)
  
  # Get number of divergent transitions
  np <- nuts_params(model)
  n_divergent <- if (!is.null(np)) sum(np$Value[np$Parameter == "divergent__"]) else NA
  
  # Get sampling diagnostics
  sampler_params <- rstan::get_sampler_params(model$fit, inc_warmup = FALSE)
  
  # Calculate max treedepth hits
  max_treedepth_hits <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] >= 10)))
  
  tibble(
    Model = model_name,
    N_obs = nobs(model),
    N_patients = n_distinct(model$data$individual_id),
    N_params = nrow(fixed_df),
    
    # Convergence: R-hat (should be < 1.01)
    Rhat_max = max(fixed_df$Rhat, na.rm = TRUE),
    Rhat_mean = mean(fixed_df$Rhat, na.rm = TRUE),
    Rhat_ok = all(fixed_df$Rhat < 1.01, na.rm = TRUE),
    
    # Effective sample size (should be > 400 for reliable inference)
    ESS_bulk_min = min(fixed_df$Bulk_ESS, na.rm = TRUE),
    ESS_tail_min = min(fixed_df$Tail_ESS, na.rm = TRUE),
    ESS_ok = min(fixed_df$Bulk_ESS, na.rm = TRUE) > 400 & 
      min(fixed_df$Tail_ESS, na.rm = TRUE) > 400,
    
    # Sampling issues
    N_divergent = n_divergent,
    N_max_treedepth = max_treedepth_hits,
    Sampling_ok = n_divergent == 0 & max_treedepth_hits == 0,
    
    # Overall quality
    Quality = case_when(
      Rhat_ok & ESS_ok & Sampling_ok ~ "✓ Good",
      !Sampling_ok ~ "⚠ Sampling issues",
      !Rhat_ok ~ "⚠ Convergence issues",
      !ESS_ok ~ "⚠ Low ESS",
      TRUE ~ "⚠ Check"
    )
  )
}

extract_draws_long <- function(model, var_recodes) {
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value") %>%
    mutate(
      submodel = if_else(str_detect(param, "^b_phi_"), "phi", "mu"),
      variable = str_remove(param, "^b_phi_") %>% str_remove("^b_")
    ) %>%
    filter(variable != "Intercept") %>%
    mutate(pretty_var = recode(variable, !!!var_recodes))
}

extract_disp_effects <- function(model, var_recodes) {
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value") %>%
    mutate(
      submodel = if_else(str_detect(param, "^b_phi_"), "phi", "mu"),
      variable = str_remove(param, "^b_phi_") %>% str_remove("^b_")
    ) %>%
    filter(variable != "Intercept") %>%
    group_by(submodel, variable) %>%
    summarise(
      estimate = median(value),
      ci_low   = quantile(value, 0.025),
      ci_high  = quantile(value, 0.975),
      pd       = max(mean(value > 0), mean(value < 0)),
      .groups  = "drop"
    ) %>%
    mutate(
      significant = pd > 0.975,
      direction   = if_else(estimate > 0, "+", "-"),
      pretty_var  = recode(variable, !!!var_recodes)
    )
}

get_effects <- function(model, model_label) {
  ps <- posterior_summary(model)
  tibble(param = rownames(ps),
         est = ps[, "Estimate"], lo = ps[, "Q2.5"], hi = ps[, "Q97.5"]) %>%
    filter(str_starts(param, "b_"), param != "b_Intercept") %>%
    mutate(term = pretty_term(param), model = model_label)
}

grab_comp <- function(model, label) {
  as_draws_df(model) %>%
    select(starts_with("b_"), starts_with("bsp_")) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value") %>%
    filter(!str_detect(param, "Intercept")) %>%
    mutate(variable = if_else(str_detect(param, "prev_theta_pm"),
                              "prev_any_resistance",
                              str_remove(param, "^b_") %>% str_remove("^bsp_")),
           model = label) %>%
    filter(variable %in% names(comp_lab)) %>%
    mutate(pretty_var = recode(variable, !!!comp_lab))
}

eff_from <- function(model, label, persist_col, kind) {
  as_draws_df(model) %>%
    select(starts_with("b_"), starts_with("bsp_")) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value") %>%
    filter(!str_detect(param, "Intercept")) %>%
    mutate(submodel = if_else(str_detect(param, "^b_phi_"), "phi", "mu"),
           var = param %>% str_remove("^b_phi_") %>% str_remove("^b_") %>% str_remove("^bsp_")) %>%
    filter(var %in% c(persist_col, common_vars)) %>%
    group_by(submodel, var) %>%
    summarise(est = median(value), lo = quantile(value, .025),
              hi = quantile(value, .975),
              pd = max(mean(value > 0), mean(value < 0)), .groups = "drop") %>%
    mutate(class = label, kind = kind,
           pred = relabel_pred(var, persist_col),
           direction = if_else(est > 0, "+", "-"),
           survives = pd > 0.975)
}

eff_sev <- function(model, label, persist_col) {
  as_draws_df(model) %>%
    select(starts_with("b_")) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value") %>%
    filter(!str_detect(param, "Intercept"), !str_detect(param, "^b_phi_")) %>%
    mutate(var = str_remove(param, "^b_")) %>%
    group_by(var) %>%
    summarise(est = median(value), lo = quantile(value, .025),
              hi = quantile(value, .975),
              pd = max(mean(value > 0), mean(value < 0)), .groups = "drop") %>%
    mutate(class = label,
           pred = dplyr::case_when(
             var == persist_col ~ "Prior same-class R",
             var == "clinical_statusAPEout" ~ "Status: APE Outpatient",
             var == "clinical_statusAPEiv"  ~ "Status: APE IV",
             var %in% common_vars ~ unname(setNames(common_labels, common_vars)[var]),
             TRUE ~ var),
           survives = pd > 0.975)
}


# =============================================================================
# HELPERS: plotting
# plot_posteriors() and plot_stratum() were duplicated in Synthesis.R and
# Talk_figures_v2.R, identical except for the +/- fill colours. Merged here
# with those two colours as arguments, defaulting to the poster palette.
# =============================================================================

create_forest_plot <- function(effects_df, title, subtitle = NULL) {
  effects_df %>%
    mutate(Variable = fct_reorder(Variable, Estimate)) %>%
    ggplot(aes(x = Estimate, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                   height = 0.2, linewidth = 0.8, color = "grey40") +
    geom_point(aes(color = Direction, size = abs(Estimate)), alpha = 0.8) +
    geom_text(aes(label = Significance), hjust = -0.5, vjust = 0.5, size = 5) +
    scale_color_manual(values = c("↑ Resistance" = "#E41A1C", 
                                  "↓ Resistance" = "#377EB8",
                                  "Uncertain" = "grey50"),
                       guide = "none") +
    scale_size_continuous(range = c(3, 6), guide = "none") +
    labs(title = title, subtitle = subtitle, x = "Log-odds (effect on resistance)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.title.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30")
    )
}

barpanel <- function(d, title, subtitle, fill) {
  ggplot(d, aes(reorder(item, pct), pct)) +
    geom_col(fill = fill, width = 0.7) +
    geom_text(aes(label = scales::percent(pct, accuracy = 1)), hjust = -0.15, size = 5) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.20))) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme(plot.title = element_text(face = "bold", size = 18),
          plot.subtitle = element_blank())
}

rx_by_prior <- function(df, rx_col, prev_col, label) {
  if (!all(c(rx_col, prev_col) %in% names(df))) return(NULL)
  df %>%
    filter(!is.na(.data[[rx_col]]), !is.na(.data[[prev_col]])) %>%
    transmute(prior = if_else(.data[[prev_col]] == 1, "Prior R", "Prior S"),
              rx    = as.integer(.data[[rx_col]] > 0)) %>%
    group_by(prior) %>%
    summarise(rate = mean(rx), .groups = "drop") %>%
    mutate(class = label)
}

perclass_forest <- function(sm, xlab, title) {
  dist_eff %>%
    filter(submodel == sm) %>%
    mutate(pred  = factor(pred, levels = rev(pred_order)),
           class = factor(class, levels = class_cfg$label)) %>%
    ggplot(aes(est, pred)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi, color = direction, alpha = survives),
                    linewidth = 0.7, size = 0.4) +
    facet_wrap(~ class, nrow = 1) +
    scale_color_manual(values = c("+" = pal[["red"]], "-" = pal[["blue"]]), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), guide = "none") +
    labs(x = xlab, y = NULL, title = title,
         caption = "Distributional beta-binomial; solid = 95% CI excludes 0.") +
    theme(plot.title = element_text(face = "bold", size = 15),
          strip.text = element_text(face = "bold"))
}

# CCM coupling diagram. Nodes = resistance classes, one bold bidirectional edge
# per pair flagged `converges`, faint edges for pairs tested but flat.
#
# Factored out so the SAME renderer serves both the raw Section 7.5 run and the
# surrogate-tested Section 12.1 run. Whatever sets `converges` decides which
# edges are bold, so the figure follows the evidence instead of being pinned to
# whichever CCM happened to run first.
#
# edges: a, b, rho, converges (one row per unordered pair)
# xy   : lab, x, y (node layout)
ccm_coupling_diagram <- function(edges, xy, title, caption, r_clear = 0.85,
                                 edge_label = "cross-map skill %.2f") {
  # Callers may hand this an edge table that already carries node coordinates
  # (Section 7.5 builds them inline). Drop any such columns first so the joins
  # below cannot produce duplicate x1/y1/x2/y2 names.
  e <- edges %>%
    dplyr::select(-dplyr::any_of(c("x", "y", "x1", "y1", "x2", "y2"))) %>%
    dplyr::left_join(xy, by = c("a" = "lab")) %>% dplyr::rename(x1 = x, y1 = y) %>%
    dplyr::left_join(xy, by = c("b" = "lab")) %>% dplyr::rename(x2 = x, y2 = y)

  bold <- e %>% dplyr::filter(converges) %>%
    dplyr::mutate(d   = sqrt((x2 - x1)^2 + (y2 - y1)^2),
                  x1c = x1 + (x2 - x1) / d * r_clear, y1c = y1 + (y2 - y1) / d * r_clear,
                  x2c = x2 - (x2 - x1) / d * r_clear, y2c = y2 - (y2 - y1) / d * r_clear,
                  xm  = (x1 + x2) / 2, ym = (y1 + y2) / 2,
                  lab = sprintf(edge_label, rho))
  faint <- e %>% dplyr::filter(!converges)

  p <- ggplot() +
    geom_segment(data = faint, aes(x = x1, y = y1, xend = x2, yend = y2),
                 color = pal[["sand"]], linewidth = 0.8, alpha = 0.55) +
    geom_point(data = xy, aes(x, y), size = 48, shape = 21,
               fill = "#AEC6CC", color = "grey30", show.legend = FALSE) +
    geom_text(data = xy, aes(x, y, label = lab), fontface = "bold", size = 4.5)

  if (nrow(bold) > 0) {
    p <- p +
      geom_segment(data = bold, aes(x = x1c, y = y1c, xend = x2c, yend = y2c),
                   color = pal[["red"]], linewidth = 1.6,
                   arrow = arrow(ends = "both", length = unit(0.22, "cm"), type = "closed")) +
      geom_text(data = bold, aes(x = xm, y = ym + 0.34, label = lab),
                color = pal[["red"]], fontface = "bold", size = 6)
  }

  p +
    coord_cartesian(xlim = c(-2.4, 2.4), ylim = c(-2.3, 2.3)) +
    labs(title = title, caption = caption) +
    theme_void() +
    theme(plot.title   = element_text(face = "bold", size = 18, hjust = 0.5),
          plot.caption = element_text(size = 10, hjust = 0.5, colour = "grey35"))
}

plot_posteriors <- function(d, title, var_levels,
                            fill_pos = unname(pal["red"]),
                            fill_neg = unname(pal["blue"])) {
  d %>%
    mutate(pretty_var = factor(pretty_var, levels = var_levels)) %>%
    group_by(pretty_var) %>%
    mutate(
      pd          = max(mean(value > 0), mean(value < 0)),
      median_est  = median(value),
      direction   = if_else(median_est > 0, "+", "-"),
      significant = pd > 0.975
    ) %>%
    ungroup() %>%
    ggplot(aes(x = value, y = pretty_var)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    stat_halfeye(
      aes(fill = direction, alpha = significant),
      .width = c(0.66, 0.95),
      slab_color = NA, point_size = 2,
      interval_size_range = c(0.6, 1.6)
    ) +
    scale_fill_manual(values = c("+" = fill_pos, "-" = fill_neg), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 0.85, `FALSE` = 0.4), guide = "none") +
    labs(x = NULL, y = NULL, title = title)
}

plot_stratum <- function(draws, sm, title,
                         fill_pos = unname(pal["red"]),
                         fill_neg = unname(pal["blue"])) {
  draws %>%
    filter(submodel == sm) %>%
    mutate(pretty_var = factor(pretty_var, levels = stratum_var_order)) %>%
    group_by(pretty_var, stratum) %>%
    mutate(
      median_est  = median(value),
      direction   = if_else(median_est > 0, "+", "-"),
      pd          = max(mean(value > 0), mean(value < 0)),
      significant = pd > 0.975
    ) %>%
    ungroup() %>%
    ggplot(aes(x = value, y = pretty_var)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    stat_halfeye(
      aes(fill = direction, alpha = significant),
      .width = c(0.66, 0.95),
      slab_color = NA, point_size = 2,
      interval_size_range = c(0.6, 1.6)
    ) +
    scale_fill_manual(values = c("+" = fill_pos, "-" = fill_neg), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 0.85, `FALSE` = 0.4), guide = "none") +
    facet_wrap(~ stratum, scales = "free_x", nrow = 1) +
    labs(x = NULL, y = NULL, title = title) +
    theme(strip.text = element_text(face = "bold", size = 14),
          plot.title = element_text(size = 16), 
          axis.text.x = element_text(size = 16), 
          axis.text.y = element_text(size = 16))
}


# =============================================================================
# HELPERS: noise / sampling design
# 
# =============================================================================

beta_post_sd <- function(k, n) {
  a <- 1 + k; b <- 1 + n - k
  sqrt(a * b / ((a + b)^2 * (a + b + 1)))
}

expected_sd_at_depth <- function(n, p) {
  ks <- 0:n
  sum(dbinom(ks, n, p) * beta_post_sd(ks, n))
}


# =============================================================================
# HELPERS: dynamics (Markov, EDM, CCM)
# 
# =============================================================================

calculate_transitions <- function(data, prev_var, curr_var, drug_name = "Overall") {
  data %>%
    filter(!is.na(.data[[prev_var]]), !is.na(.data[[curr_var]])) %>%
    mutate(
      prev_state = if_else(.data[[prev_var]] == 1, "R", "S"),
      current_state = if_else(.data[[curr_var]] == 1, "R", "S")
    ) %>%
    count(prev_state, current_state) %>%
    group_by(prev_state) %>%
    mutate(
      total = sum(n),
      prob = n / total,
      se = sqrt(prob * (1 - prob) / total),
      ci_lower = pmax(0, prob - 1.96 * se),
      ci_upper = pmin(1, prob + 1.96 * se)
    ) %>%
    ungroup() %>%
    mutate(drug = drug_name)
}

markov_one <- function(class, prev, cur) {
  d  <- pa_visit_full %>% filter(!is.na(.data[[prev]]), !is.na(.data[[cur]]))
  pv <- d[[prev]]; cu <- d[[cur]]
  tibble(class = class,
         p_SR = mean(cu[pv == 0] == 1),   # acquire S -> R
         p_RS = mean(cu[pv == 1] == 0),   # lose    R -> S
         p_RR = 1 - mean(cu[pv == 1] == 0),
         p_SS = 1 - mean(cu[pv == 0] == 1),
         steady_R = mean(cu[pv == 0] == 1) /
           (mean(cu[pv == 0] == 1) + mean(cu[pv == 1] == 0)),
         observed = mean(cu == 1), n_trans = nrow(d))
}

find_optimal_params <- function(df, target_col, max_E = 12, lib_frac = 0.7) {
  
  df_clean <- df %>% filter(!is.na(.data[[target_col]]))
  n_obs <- nrow(df_clean)
  
  if (n_obs < 100) {
    cat("  Insufficient data for", target_col, "(N =", n_obs, ")\n")
    return(list(E = NA, theta = NA, rho_E = NA, rho_theta = NA, N = n_obs))
  }
  
  lib_end <- floor(n_obs * lib_frac)
  lib_str <- paste("1", lib_end)
  pred_str <- paste(lib_end + 1, n_obs)
  
  # Step 1: Find optimal E
  embed_out <- tryCatch({
    EmbedDimension(
      dataFrame = as.data.frame(df_clean),
      columns = target_col,
      target = target_col,
      lib = lib_str,
      pred = pred_str,
      maxE = max_E,
      # showPlot defaults to TRUE, which opens an implicit device. Under
      # Rscript that is Rplots.pdf in the working directory - and if a second
      # Rscript is running in the same directory both processes write to it,
      # the device errors with "write failed", and every later plot in the
      # session dies with it. Never let rEDM draw.
      showPlot = FALSE
    )
  }, error = function(e) {
    cat("  EmbedDimension error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(embed_out)) {
    return(list(E = NA, theta = NA, rho_E = NA, rho_theta = NA, N = n_obs))
  }
  
  best_E <- embed_out$E[which.max(embed_out$rho)]
  best_rho_E <- max(embed_out$rho, na.rm = TRUE)
  
  # Step 2: Find optimal theta using SMap
  theta_seq <- seq(0, 8, by = 0.5)
  
  theta_results <- map_dfr(theta_seq, function(th) {
    smap_result <- tryCatch({
      SMap(
        dataFrame = as.data.frame(df_clean),
        columns = target_col,
        target = target_col,
        lib = lib_str,
        pred = pred_str,
        E = best_E,
        theta = th
      )
    }, error = function(e) NULL)
    
    if (!is.null(smap_result)) {
      rho <- NA
      if ("stats" %in% names(smap_result) && "rho" %in% names(smap_result$stats)) {
        rho <- smap_result$stats$rho
      } else if ("predictions" %in% names(smap_result)) {
        preds <- smap_result$predictions
        pred_col <- names(preds)[grep("Pred|pred", names(preds))[1]]
        obs_col <- names(preds)[grep("Obs|obs", names(preds))[1]]
        if (!is.na(pred_col) && !is.na(obs_col)) {
          valid <- !is.na(preds[[pred_col]]) & !is.na(preds[[obs_col]])
          if (sum(valid) > 10) {
            rho <- cor(preds[[pred_col]][valid], preds[[obs_col]][valid])
          }
        }
      }
      tibble(theta = th, rho = rho)
    } else {
      tibble(theta = th, rho = NA)
    }
  })
  
  best_theta <- theta_results$theta[which.max(theta_results$rho)]
  best_rho_theta <- max(theta_results$rho, na.rm = TRUE)
  
  return(list(
    E = best_E,
    theta = best_theta,
    rho_E = best_rho_E,
    rho_theta = best_rho_theta,
    embed_data = embed_out,
    theta_data = theta_results,
    N = n_obs
  ))
}

find_E <- function(df, v, maxE = 8) {
  d <- df %>% filter(!is.na(.data[[v]]))
  if (nrow(d) < 60) return(3L)
  rng <- paste(1, nrow(d))   # leave-one-out over the whole series (not a
  # patient-splitting train/test, which returns NA)
  out <- tryCatch(EmbedDimension(dataFrame = as.data.frame(d), columns = v, target = v,
                                 lib = rng, pred = rng, maxE = maxE, showPlot = FALSE),
                  error = function(e) NULL)
  if (is.null(out) || all(is.na(out$rho))) return(3L)
  E <- out$E[which.max(out$rho)]
  if (length(E) == 0 || is.na(E)) return(3L)
  as.integer(E)
}

# Separator-aware E selection. find_E() filters out NA rows, which would strip
# the very NA separators that stop an embedding vector spanning two patients -
# so E ends up chosen on a boundary-crossing frame even when the CCM that uses
# it is boundary-safe. This keeps the separators and lets rEDM drop the
# crossing vectors itself. E propagates into every CCM call, including the
# surrogate tests, so it has to be selected on the same frame geometry.
find_E_sep <- function(df_sep, v, maxE = 8) {
  n_ok <- sum(!is.na(df_sep[[v]]))
  if (n_ok < 60) return(3L)
  rng <- paste(1, nrow(df_sep))
  out <- tryCatch(suppressWarnings(
    EmbedDimension(dataFrame = as.data.frame(df_sep), columns = v, target = v,
                   lib = rng, pred = rng, maxE = maxE, showPlot = FALSE)),
    error = function(e) NULL)
  if (is.null(out) || all(is.na(out$rho))) return(3L)
  E <- out$E[which.max(out$rho)]
  if (length(E) == 0 || is.na(E)) return(3L)
  as.integer(E)
}

run_ccm_pair <- function(df, a, b, E_map) {
  d <- df %>% filter(!is.na(.data[[a]]), !is.na(.data[[b]])) %>% as.data.frame()
  if (nrow(d) < 80) { cat("  skip", a, "/", b, "- N =", nrow(d), "\n"); return(NULL) }
  lib_max   <- nrow(d) - 1
  lib_sizes <- unique(round(seq(25, lib_max, length.out = 8)))
  lib_sizes <- lib_sizes[lib_sizes >= 10 & lib_sizes <= lib_max]
  Ea <- E_map[[a]]; Eb <- E_map[[b]]
  if (length(Ea) == 0 || is.na(Ea)) Ea <- 3L
  if (length(Eb) == 0 || is.na(Eb)) Eb <- 3L
  c_ab <- tryCatch(CCM(dataFrame = d, columns = b, target = a, E = Eb,
                       libSizes = lib_sizes, sample = 100), error = function(e) NULL)  # a -> b
  c_ba <- tryCatch(CCM(dataFrame = d, columns = a, target = b, E = Ea,
                       libSizes = lib_sizes, sample = 100), error = function(e) NULL)  # b -> a
  if (is.null(c_ab) || is.null(c_ba)) return(NULL)
  bind_rows(
    tibble(LibSize = c_ab$LibSize, rho = c_ab[[paste0(b, ":", a)]], driver = a, affected = b),
    tibble(LibSize = c_ba$LibSize, rho = c_ba[[paste0(a, ":", b)]], driver = b, affected = a)
  )
}

classify_ccm <- function(causal_flag, confounded_flag, glmm_sig) {
  # Both flags missing means no CCM row exists for this variable - the CCM was
  # never run for it. That is NOT the same as the CCM running and finding
  # nothing, and it must not be folded into "confounded"/"null": doing so
  # manufactures a causal verdict out of absent evidence. This fires whenever
  # ccm_summary_overall is absent, which is the normal state of this pipeline
  # now that the drug-level placeholder CCM has been dropped.
  if (is.na(causal_flag) && is.na(confounded_flag)) return("not_assessed")

  forward <- !is.na(causal_flag)     & causal_flag == "✓"
  reverse <- !is.na(confounded_flag) & confounded_flag == "⚠"
  # extract_effects() now reports probability of direction rather than stars.
  # The bar here is pd > 0.95; the pd > 0.90 tier is deliberately NOT counted as
  # evidence, since a one-sided p of 0.10 should not drive a causal label.
  sig_glmm <- glmm_sig %in% c("pd>0.95", "pd>0.975", "pd>0.99")

  if (forward && reverse)          return("bidirectional")
  if (!forward && reverse)         return("reverse")
  if (forward && !reverse)         return("one_way")
  if (!forward && !reverse && sig_glmm) return("confounded")
  "null"
}


# =============================================================================
# HELPERS: external deep-sampling data (Vanderwoude et al. 2024)
# 
# =============================================================================

read_zoi <- function(sheet, id_label) {
  read_excel(supp_path, sheet = sheet) %>%
    filter(!is.na(Isolate)) %>%
    transmute(isolate = Isolate,
              `Aminoglycoside (amikacin)`   = `AK AVG`,
              `Aminoglycoside (tobramycin)` = `TOB AVG`,
              `Carbapenem (meropenem)`      = `MEM AVG`,
              `Cephalosporin (ceftazidime)` = `CAZ AVG`,
              `Fluoroquinolone (cipro)`     = `CIP AVG`,
              `Penicillin (pip-tazo)`       = `TZP AVG`) %>%
    pivot_longer(-isolate, names_to = "drug", values_to = "zoi") %>%
    mutate(patient = id_label)
}

growth_mu <- function(sheet, id_label) {
  raw    <- read_excel(supp_path, sheet = sheet, col_names = FALSE)
  labels <- trimws(as.character(unlist(raw[1, ])))
  tvec   <- suppressWarnings(as.numeric(unlist(raw[-1, 1])))
  M      <- suppressWarnings(matrix(as.numeric(as.matrix(raw[-1, ])), nrow = length(tvec)))
  blank  <- rowMeans(M[, which(tolower(labels) == "blank"), drop = FALSE], na.rm = TRUE)
  isos   <- setdiff(unique(labels[-1]), c("blank", "NA", ""))
  sm3 <- function(x) { n <- length(x); y <- x
  for (i in 2:(n - 1)) y[i] <- mean(x[(i - 1):(i + 1)]); y }
  mu <- vapply(isos, function(lab) {
    od <- rowMeans(M[, which(labels == lab), drop = FALSE], na.rm = TRUE) - blank
    od <- pmax(sm3(od), 1e-3); ln <- log(od); w <- 4; best <- NA_real_
    for (i in seq_len(length(tvec) - w)) {
      if (od[i + w] < 0.02) next
      s <- unname(coef(lm(ln[i:(i + w)] ~ tvec[i:(i + w)]))[2])
      if (is.na(best) || s > best) best <- s
    }
    best
  }, numeric(1))
  tibble(isolate = isos, mu_max = mu, patient = id_label)
}


# ---- end of AMR_common.R ----------------------------------------------------
