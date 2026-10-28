# =============================================================================
# AMR_pipeline.R  -  P. aeruginosa resistance dynamics in cystic fibrosis
# Emory pilot cohort | conservative (R + I) coding
#
# Author: Canan Karakoc
# Assembled: August 25, 2026
#
# Upstream: analytical_data_prep.R cleans the raw csv and writes
#           data/pa_visit_full.* and data/pa_isolates.*; data_inspection_summary.R
#           produces the descriptive pass. This file is the ANALYSIS layer and
#           never touches the raw data. Run order: ./run_pipeline.sh all
# Sideways: AMR_common.R holds packages, palette, themes, priors and every
#           helper function. Edit that file rarely; edit this one freely.
#
# -----------------------------------------------------------------------------
# CONTENTS
# -----------------------------------------------------------------------------
#   0  Setup and analysis frames
#   1  Core GLMMs           M1/M2 summary, M3-M6 stratified, M7 cross-resistance,
#                           effect extraction, forests, colistin, diagnostics
#   2  Exposure noise       route intensity vs binary prescribed; confounding
#   3  Sampling noise       variance decomposition (2- and 4-class), precision
#                           by sampling depth, down-sampling propagation
#   4  Heterogeneity        mu vs phi, naive vs measurement-error, staph strata
#   5  Per-class models     standard, deep dive, clinical-status-aware
#   6  Ecology              exacerbation vs route, co-colonizers, time trend
#   7  Dynamics             Markov chains, EDM/S-map, convergent cross-mapping
#   8  Deep-sampling study  Vanderwoude et al. 2024 external validation
#   9  Summary              cohort panel, integrated comparison, key findings
#  10  Design sensitivity   colony denominators vs drug-test denominators, an
#                           ESTIMATED between-colony variance in place of the
#                           pi^2/3 convention, and sampling depth as a covariate.
#                           Added after review; Sections 1-9 are unchanged, so
#                           the original design stays quotable alongside it.
#  11  Review responses II  CCM across patient boundaries, calendar era, pooled
#                           cross-resistance, continuous-time Markov, and
#                           measurement error on the logit scale
#  12  Review responses III CCM against a surrogate null (base vs surrogate),
#                           and me() on the class-specific prior-resistance terms
#  13  Multispatial CCM     the published Clark et al. 2015 implementation of
#                           the stitch-the-patients design, with bootstrap CIs
#  14  Surveillance design  depth -> reliability, interval -> autocorrelation,
#                           length x cohort -> power; what it would take
#
# -----------------------------------------------------------------------------
# NOISE SOURCES (the narrative backbone, carried over from Talk_figures_v2)
# -----------------------------------------------------------------------------
#   1 EXPOSURE        prescribed != taken/absorbed      -> Section 2 (null)
#   2 SAMPLING        few colonies != population         -> Section 3 (dominant)
#   3 CLASSIFICATION  AST breakpoints / intermediate     -> R + I coding
#   4 TIMING          visits cluster at exacerbations    -> Section 6 caveat
#
# -----------------------------------------------------------------------------
# HOW TO RUN
# -----------------------------------------------------------------------------
# Source top to bottom, or run Section 0 then any later section on its own.
# brms fits are cached via file=/file_refit = "on_change", so the first run is
# slow and re-runs are fast. NOTE: the fixes applied to analytical_data_prep.R
# changed the cleaned data, so the first run after re-prepping WILL refit every
# model.
#
# TWO EXPENSIVE BLOCKS, and they are the only ones that are expensive once the
# brms caches are warm: 12.1 (CCM surrogate null, N_SURR) and 14.3 (design power
# grid, DESIGN_REPS). Together they are roughly a 9-hour run at publication
# settings, and 12.1 runs BEFORE Sections 13 and 14 - so a plain error anywhere
# downstream of it would surface five hours in.
#
# To check the whole file end to end first, run it with AMR_QUICK set:
#     AMR_QUICK=1 ./run_pipeline.sh pipeline
# That drops both constants to a structural-check size (~20 min) and exercises
# every code path. It is NOT a result: both caches are stamped with the constant
# that produced them, so the publication run detects the quick cache as stale
# and recomputes it. A quick number can never silently end up in a figure.
#
# =============================================================================


# #############################################################################
# SECTION 0: SETUP AND ANALYSIS FRAMES
# AMR_common.R brings the packages, palette, themes, priors and helpers.
# The filtering below is analysis-specific subsetting, not data cleaning -
# cleaning lives in analytical_data_prep.R.
# #############################################################################

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
source("AMR_common.R")
dir.create("figures/glmm", showWarnings = FALSE)
dir.create("figures/markov", showWarnings = FALSE)

# Structural-check mode. See HOW TO RUN at the top: this only shrinks the two
# Monte Carlo constants so the file can be exercised end to end in minutes. The
# 12.1 and 14.3 caches are stamped with those constants, so a quick run cannot
# leave a low-precision result behind for a later run to pick up.
AMR_QUICK <- nzchar(Sys.getenv("AMR_QUICK"))
if (AMR_QUICK) {
  cat("\n*** AMR_QUICK: structural check only. N_SURR and DESIGN_REPS reduced.\n")
  cat("*** Numbers from this run are NOT publication precision. Re-run without\n")
  cat("*** AMR_QUICK to recompute both caches at full size.\n\n")
}

.dat          <- load_amr_data()
pa_visit_full <- .dat$pa_visit_full
pa_isolates   <- .dat$pa_isolates
rm(.dat)

# The four classes with enough testing depth to support a variance or design
# calculation. Defined HERE rather than in Section 14 because 10.3c uses it too,
# and it was previously declared ~1800 lines after its first use - which meant
# Section 10.3c only ran in a session that already had the object lying around
# from a previous run, and died in a clean one.
design_classes <- c(ag = "Aminoglycoside", carb = "Carbapenem",
                    ceph = "Cephalosporin", fq = "Fluoroquinolone")

pa_analysis <- pa_visit_full %>%
  filter(
    !is.na(prev_any_resistance),
    !is.na(n_tested),
    n_tested > 0
  ) %>%
  mutate(
    # Current-visit exacerbation (kept only for backwards compatibility;
    # do NOT use as a predictor — see prev_has_exacerbation below)
    has_exacerbation = as.integer(clinical_status %in% c("APE Outpatient", "APE IV Antibiotics")),
    
    # Lagged exacerbation at t-1 — used as a predictor throughout the pipeline
    # to keep the t-1 predictors → t outcome temporal structure consistent with
    # every other covariate. NA (first visit per patient) coerced to 0.
    prev_has_exacerbation = as.integer(
      prev_clinical_status %in% c("APE Outpatient", "APE IV Antibiotics")
    ),
    prev_has_exacerbation = replace_na(prev_has_exacerbation, 0L),
    
    has_any_coinfection = as.integer(n_coinfections > 0),
    
    # Binary outcomes for 4 classes (EXCLUDING colistin as outcome)
    carb_resistant = as.integer(carb_any_resistant_conservative == 1),
    ag_resistant = as.integer(ag_any_resistant_conservative == 1),
    fq_resistant = as.integer(fq_any_resistant_conservative == 1),
    ceph_resistant = as.integer(ceph_any_resistant_conservative == 1)
    # NOTE: coli_resistant NOT created - insufficient outcome data
  )

pa_ecological <- pa_analysis
pa_prescription <- pa_analysis %>% filter(!is.na(prev_aminoglycoside_prescribed))

cat("Ecological analysis N:", nrow(pa_ecological), "\n")
cat("Prescription analysis N:", nrow(pa_prescription), "\n\n")

# -----------------------------------------------------------------------------
# SAMPLE-SIZE SANITY CHECK
# -----------------------------------------------------------------------------
# The prescription analysis subset is substantially smaller than the ecological
# subset because prescription columns (prev_*_prescribed) are only recorded at
# visits where the patient had a pulmonary exacerbation (APE) at t-1. Stable
# visits are structurally NA on the prescription columns, so
# `filter(!is.na(prev_aminoglycoside_prescribed))` implicitly drops all
# t-1-stable visits.
#
# Consequence: within pa_prescription, prev_has_exacerbation is CONSTANT (= 1)
# and therefore drops out of any GLMM fit on this subset. We include
# prev_has_exacerbation in M1 (fit on pa_ecological, where it varies) but
# NOT in M2, M3-M7, m2_enhanced, or the EDM data frame (all built from
# pa_prescription).
#
# When reporting sample sizes, always distinguish:
#   M1 (ecological)     -> nrow(pa_ecological)
#   M2..M7 (drug Rx)    -> nrow(pa_prescription)
#   Stratified (M3-M6)  -> further reduced by prior-resistance split
# -----------------------------------------------------------------------------

stopifnot(nrow(pa_ecological) >= nrow(pa_prescription))

# Detect and warn if the APE-only prescription pattern shifts in future data.
exac_variation_in_rx <- pa_prescription %>%
  mutate(exac = as.integer(prev_clinical_status %in% 
                             c("APE Outpatient", "APE IV Antibiotics"))) %>%
  pull(exac) %>%
  unique() %>%
  length()

if (exac_variation_in_rx > 1) {
  message("NOTE: prev_has_exacerbation now varies within pa_prescription ",
          "(unique values: ", exac_variation_in_rx, "). The rationale for ",
          "excluding it from drug-Rx models may no longer apply.")
} else {
  message("Confirmed: prev_has_exacerbation is constant within pa_prescription ",
          "(all rows exacerbation = 1). Excluded from drug-Rx models by design.")
}

# Cohort strapline (patients | isolates | visits | span). Talk_figures_v2 had
# these numbers hard-coded into a figure caption; recomputed here instead.
# Must come after the frames above - it reads pa_analysis and pa_prescription.
cohort_facts()


# #############################################################################
# SECTION 1: CORE GLMMs
# Asymmetric colistin design: 4 classes as OUTCOMES, 5 as PREDICTORS.
# #############################################################################


# -----------------------------------------------------------------------------
# 1.1  Summary models M1-M2
# 
# -----------------------------------------------------------------------------

# The conflicts_prefer() block that used to sit here is gone, along with the
# library(extraDistr) call further down that made it necessary. extraDistr is
# only a brms *Suggests* - brms reaches its distributions through :: and never
# needed it attached - and attaching it collided with brms on ddirichlet,
# dfrechet, rdirichlet and friends. Each collision then had to be silenced one
# at a time, because rstan calls utils::apropos() while compiling a model,
# apropos() calls exists() across every attached package, and conflicted aborts
# the fit on the first ambiguous name it meets. Dropping extraDistr removes the
# whole class of failure; Talk_figures*.R never loaded either package.

# M1: Ecological (unchanged)
m1_ecological <- brm(
  n_resistant | trials(n_tested) ~ 
    prev_any_resistance +
    prev_has_mrsa +
    prev_has_mssa +
    prev_has_any_mucoid +
    prev_has_exacerbation +
    (1 | individual_id),
  data = pa_ecological,
  family = beta_binomial(),
  prior = priors_beta_random,
  chains = 4, cores = 4, 
  iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95),
  seed = 42,
  file = "models/m1_ecological_betabin_re",
  file_refit = "on_change"
)

# M2: Overall - KEEPS all 5 drug classes as PREDICTORS (including colistin)
# prev_has_exacerbation intentionally OMITTED: it is constant (= 1) in
# pa_prescription because prescription columns are only recorded for APE
# visits (see sanity-check block above). Including it causes rank deficiency;
# the posterior collapses to the prior (~±3 log-odds CI).
m2_overall <- brm(
  n_resistant | trials(n_tested) ~ 
    prev_any_resistance +
    prev_has_mrsa +
    prev_has_mssa +
    prev_has_any_mucoid +
    prev_carbapenem_prescribed +
    prev_aminoglycoside_prescribed +
    prev_fluoroquinolones_prescribed +
    prev_cephalosporin_prescribed +
    prev_colistin_prescribed +    # KEPT: colistin as predictor
    (1 | individual_id),
  family = beta_binomial(),
  prior = priors_beta_random,
  data = pa_prescription,
  chains = 4, cores = 4,
  iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95),
  seed = 42,
  file = "models/m2_overall_5drugs",
  file_refit = "on_change"
)

# -----------------------------------------------------------------------------
# 1.2  Stratified models M3-M6
# 
# -----------------------------------------------------------------------------

# =============================================================================

# M3: Carbapenem
carb_models <- fit_stratified_models(
  data = pa_prescription,
  class_name = "CARBAPENEM",
  outcome_var = "carb_resistant",
  prev_resist_var = "prev_carb_any_resistant",
  prev_rx_var = "prev_carbapenem_prescribed",
  file_prefix = "m3_carb"
)

# M4: Aminoglycoside
ag_models <- fit_stratified_models(
  data = pa_prescription,
  class_name = "AMINOGLYCOSIDE",
  outcome_var = "ag_resistant",
  prev_resist_var = "prev_ag_any_resistant",
  prev_rx_var = "prev_aminoglycoside_prescribed",
  file_prefix = "m4_ag"
)

# M5: Fluoroquinolone
fq_models <- fit_stratified_models(
  data = pa_prescription,
  class_name = "FLUOROQUINOLONE",
  outcome_var = "fq_resistant",
  prev_resist_var = "prev_fq_any_resistant",
  prev_rx_var = "prev_fluoroquinolones_prescribed",
  file_prefix = "m5_fq"
)

# M6: Cephalosporin
ceph_models <- fit_stratified_models(
  data = pa_prescription,
  class_name = "CEPHALOSPORIN",
  outcome_var = "ceph_resistant",
  prev_resist_var = "prev_ceph_any_resistant",
  prev_rx_var = "prev_cephalosporin_prescribed",
  file_prefix = "m6_ceph"
)

cat("NOTE: Colistin stratified models (M7) skipped\n")
cat("Reason: Insufficient colistin resistance testing data (N=14 resistant cases)\n")
cat("Colistin RETAINED as predictor in M2 and cross-resistance models\n")


# -----------------------------------------------------------------------------
# 1.3  Cross-resistance models M7
# 
# -----------------------------------------------------------------------------

# =============================================================================

pa_cross <- pa_prescription %>%
  filter(
    !is.na(carb_resistant),
    !is.na(ag_resistant),
    !is.na(fq_resistant),
    !is.na(ceph_resistant)
    # NOTE: coli_resistant filter REMOVED - not used as outcome
  )

cat("Cross-resistance analysis N:", nrow(pa_cross), "\n\n")

# Helper function for cross-resistance models

# -----------------------------------------------------------------------------
# Carbapenem exposure → Other resistances (3 models, excludes Coli outcome)
# -----------------------------------------------------------------------------

m7a_carb_to_ag <- fit_cross_model(
  "ag_resistant", "prev_ag_any_resistant",
  "prev_carbapenem_prescribed", "prev_aminoglycoside_prescribed",
  "M7A: Carb→AG", "m7a_carb_to_ag"
)

m7b_carb_to_fq <- fit_cross_model(
  "fq_resistant", "prev_fq_any_resistant",
  "prev_carbapenem_prescribed", "prev_fluoroquinolones_prescribed",
  "M7B: Carb→FQ", "m7b_carb_to_fq"
)

m7c_carb_to_ceph <- fit_cross_model(
  "ceph_resistant", "prev_ceph_any_resistant",
  "prev_carbapenem_prescribed", "prev_cephalosporin_prescribed",
  "M7C: Carb→Ceph", "m7c_carb_to_ceph"
)

# M7D: Carb→Coli REMOVED - insufficient coli_resistant outcome data

# -----------------------------------------------------------------------------
# Aminoglycoside exposure → Other resistances (3 models)
# -----------------------------------------------------------------------------

m7d_ag_to_carb <- fit_cross_model(
  "carb_resistant", "prev_carb_any_resistant",
  "prev_aminoglycoside_prescribed", "prev_carbapenem_prescribed",
  "M7D: AG→Carb", "m7d_ag_to_carb"
)

m7e_ag_to_fq <- fit_cross_model(
  "fq_resistant", "prev_fq_any_resistant",
  "prev_aminoglycoside_prescribed", "prev_fluoroquinolones_prescribed",
  "M7E: AG→FQ", "m7e_ag_to_fq"
)

m7f_ag_to_ceph <- fit_cross_model(
  "ceph_resistant", "prev_ceph_any_resistant",
  "prev_aminoglycoside_prescribed", "prev_cephalosporin_prescribed",
  "M7F: AG→Ceph", "m7f_ag_to_ceph"
)

# AG→Coli REMOVED

# -----------------------------------------------------------------------------
# Fluoroquinolone exposure → Other resistances (3 models)
# -----------------------------------------------------------------------------

m7g_fq_to_carb <- fit_cross_model(
  "carb_resistant", "prev_carb_any_resistant",
  "prev_fluoroquinolones_prescribed", "prev_carbapenem_prescribed",
  "M7G: FQ→Carb", "m7g_fq_to_carb"
)

m7h_fq_to_ag <- fit_cross_model(
  "ag_resistant", "prev_ag_any_resistant",
  "prev_fluoroquinolones_prescribed", "prev_aminoglycoside_prescribed",
  "M7H: FQ→AG", "m7h_fq_to_ag"
)

m7i_fq_to_ceph <- fit_cross_model(
  "ceph_resistant", "prev_ceph_any_resistant",
  "prev_fluoroquinolones_prescribed", "prev_cephalosporin_prescribed",
  "M7I: FQ→Ceph", "m7i_fq_to_ceph"
)

# FQ→Coli REMOVED

# -----------------------------------------------------------------------------
# Cephalosporin exposure → Other resistances (3 models)
# -----------------------------------------------------------------------------

m7j_ceph_to_carb <- fit_cross_model(
  "carb_resistant", "prev_carb_any_resistant",
  "prev_cephalosporin_prescribed", "prev_carbapenem_prescribed",
  "M7J: Ceph→Carb", "m7j_ceph_to_carb"
)

m7k_ceph_to_ag <- fit_cross_model(
  "ag_resistant", "prev_ag_any_resistant",
  "prev_cephalosporin_prescribed", "prev_aminoglycoside_prescribed",
  "M7K: Ceph→AG", "m7k_ceph_to_ag"
)

m7l_ceph_to_fq <- fit_cross_model(
  "fq_resistant", "prev_fq_any_resistant",
  "prev_cephalosporin_prescribed", "prev_fluoroquinolones_prescribed",
  "M7L: Ceph→FQ", "m7l_ceph_to_fq"
)

# Ceph→Coli REMOVED

# -----------------------------------------------------------------------------
# Colistin exposure → Other resistances (4 models) - KEPT!
# Key question: Does colistin use promote resistance to other classes?
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# M7A-strict: is the cross-resistance COEFFICIENT coding-invariant too?
# -----------------------------------------------------------------------------
# AMR_sir_ordinal.R shows the VARIANCE STRUCTURE does not care whether S/I/R is
# collapsed as R+I or strict R. That result comes from a drug-only model, so it
# licenses a claim about variance components and nothing more - and the
# proportional-odds check in that file shows the two thresholds respond to drug
# differently, so coefficients need not be invariant.
#
# This refits ONE cross-resistance model under strict coding (R only, I counted
# as susceptible) on the same rows. If the coefficient holds, the invariance
# claim can be stated broadly; if it moves, it stays a statement about variance.
# The strict LAG does not exist upstream, and it cannot be built inside
# pa_cross - that frame is already filtered, so lag() there would pair visits
# that are not actually consecutive. Build it on the full visit series the way
# analytical_data_prep.R builds the conservative lags, then join.
strict_lags <- pa_visit_full %>%
  arrange(individual_id, date_collected_from_subject) %>%
  group_by(individual_id) %>%
  mutate(prev_ag_any_res_strict   = lag(as.integer(ag_n_resistant_strict   > 0)),
         prev_carb_any_res_strict = lag(as.integer(carb_n_resistant_strict > 0))) %>%
  ungroup() %>%
  select(individual_id, date_collected_from_subject,
         prev_ag_any_res_strict, prev_carb_any_res_strict)

pa_cross_strict <- pa_cross %>%
  left_join(strict_lags, by = c("individual_id", "date_collected_from_subject")) %>%
  mutate(ag_resistant_strict = as.integer(ag_n_resistant_strict > 0)) %>%
  filter(!is.na(ag_resistant_strict), !is.na(prev_ag_any_res_strict))

cat("  M7A-strict frame: N =", nrow(pa_cross_strict), "of", nrow(pa_cross), "cross rows\n")

m7a_strict <- brm(
  ag_resistant_strict ~ prev_ag_any_res_strict + prev_carbapenem_prescribed +
    prev_aminoglycoside_prescribed + prev_has_mrsa + prev_has_mssa +
    prev_has_any_mucoid + (1 | individual_id),
  data = pa_cross_strict, family = bernoulli(), prior = priors_bernoulli,
  chains = 4, cores = 4, iter = 2000, warmup = 1000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/m7a_carb_to_ag_strict", file_refit = "on_change"
)

cat("\n=== M7A coding sensitivity: conservative (R+I) vs strict (R only) ===\n")
.cmp <- bind_rows(
  fixef(m7a_carb_to_ag) %>% as.data.frame() %>% rownames_to_column("term") %>%
    mutate(coding = "conservative (R + I)"),
  fixef(m7a_strict) %>% as.data.frame() %>% rownames_to_column("term") %>%
    mutate(coding = "strict (R only)")) %>%
  mutate(term = str_remove(term, "_strict$"),
         term = str_replace(term, "prev_ag_any_res", "prev_ag_any_resistant")) %>%
  select(term, coding, Estimate, Q2.5, Q97.5)
print(.cmp %>% mutate(across(where(is.numeric), ~round(., 2))) %>% as.data.frame(),
      row.names = FALSE)
write_csv(.cmp, "results/m7a_coding_sensitivity.csv")
cat("  Coefficients agreeing within their intervals -> the invariance claim\n")
cat("  extends past variance components. Diverging -> keep it to variance.\n\n")

m7m_coli_to_carb <- fit_cross_model(
  "carb_resistant", "prev_carb_any_resistant",
  "prev_colistin_prescribed", "prev_carbapenem_prescribed",
  "M7M: Coli→Carb", "m7m_coli_to_carb"
)

m7n_coli_to_ag <- fit_cross_model(
  "ag_resistant", "prev_ag_any_resistant",
  "prev_colistin_prescribed", "prev_aminoglycoside_prescribed",
  "M7N: Coli→AG", "m7n_coli_to_ag"
)

m7o_coli_to_fq <- fit_cross_model(
  "fq_resistant", "prev_fq_any_resistant",
  "prev_colistin_prescribed", "prev_fluoroquinolones_prescribed",
  "M7O: Coli→FQ", "m7o_coli_to_fq"
)

m7p_coli_to_ceph <- fit_cross_model(
  "ceph_resistant", "prev_ceph_any_resistant",
  "prev_colistin_prescribed", "prev_cephalosporin_prescribed",
  "M7P: Coli→Ceph", "m7p_coli_to_ceph"
)


# -----------------------------------------------------------------------------
# 1.4  Effect extraction
# 
# -----------------------------------------------------------------------------



# Parameter mappings
# NOTE: Only M1 (fit on pa_ecological) retains prev_has_exacerbation because
# the variable is constant within pa_prescription. If you later extend M1 to
# use prev_clinical_status as a 3-level factor (Stable / APE Outpatient /
# APE IV), add the two treatment-contrast entries here:
#   "b_prev_clinical_statusAPEOutpatient"   = "APE Outpatient (vs Stable)"
#   "b_prev_clinical_statusAPEIVAntibiotics" = "APE IV (vs Stable)"
params_m1 <- c(
  "b_prev_any_resistance" = "Prior Resistance",
  "b_prev_has_mrsa" = "MRSA Co-colonization",
  "b_prev_has_mssa" = "MSSA Co-colonization",
  "b_prev_has_any_mucoid" = "Mucoid Phenotype",
  "b_prev_has_exacerbation" = "Exacerbation"
)

params_m2 <- c(
  "b_prev_any_resistance" = "Prior Resistance",
  "b_prev_has_mrsa" = "MRSA Co-colonization",
  "b_prev_has_mssa" = "MSSA Co-colonization",
  "b_prev_has_any_mucoid" = "Mucoid Phenotype",
  "b_prev_carbapenem_prescribed" = "Carbapenem Rx",
  "b_prev_aminoglycoside_prescribed" = "Aminoglycoside Rx",
  "b_prev_fluoroquinolones_prescribed" = "Fluoroquinolone Rx",
  "b_prev_cephalosporin_prescribed" = "Cephalosporin Rx",
  "b_prev_colistin_prescribed" = "Colistin Rx"
)

params_acquisition <- c(
  "b_prev_any_resistance" = "Prior Overall Resistance",
  "b_prev_carbapenem_prescribed" = "Carbapenem Rx",
  "b_prev_aminoglycoside_prescribed" = "Aminoglycoside Rx",
  "b_prev_fluoroquinolones_prescribed" = "Fluoroquinolone Rx",
  "b_prev_cephalosporin_prescribed" = "Cephalosporin Rx",
  "b_prev_colistin_prescribed" = "Colistin Rx",
  "b_prev_has_mrsa" = "MRSA Co-colonization",
  "b_prev_has_mssa" = "MSSA Co-colonization",
  "b_prev_has_any_mucoid" = "Mucoid Phenotype"
)

params_persistence <- c(
  "b_prev_carbapenem_prescribed" = "Carbapenem Rx",
  "b_prev_aminoglycoside_prescribed" = "Aminoglycoside Rx",
  "b_prev_fluoroquinolones_prescribed" = "Fluoroquinolone Rx",
  "b_prev_cephalosporin_prescribed" = "Cephalosporin Rx",
  "b_prev_colistin_prescribed" = "Colistin Rx",
  "b_prev_has_mrsa" = "MRSA Co-colonization",
  "b_prev_has_mssa" = "MSSA Co-colonization",
  "b_prev_has_any_mucoid" = "Mucoid Phenotype"
)

params_cross <- c(
  "b_prev_carb_any_resistant" = "Prior Carb Resistance",
  "b_prev_ag_any_resistant" = "Prior AG Resistance",
  "b_prev_fq_any_resistant" = "Prior FQ Resistance",
  "b_prev_ceph_any_resistant" = "Prior Ceph Resistance",
  "b_prev_carbapenem_prescribed" = "Carbapenem Rx",
  "b_prev_aminoglycoside_prescribed" = "Aminoglycoside Rx",
  "b_prev_fluoroquinolones_prescribed" = "Fluoroquinolone Rx",
  "b_prev_cephalosporin_prescribed" = "Cephalosporin Rx",
  "b_prev_colistin_prescribed" = "Colistin Rx",
  "b_prev_has_mrsa" = "MRSA Co-colonization",
  "b_prev_has_mssa" = "MSSA Co-colonization",
  "b_prev_has_any_mucoid" = "Mucoid Phenotype"
)

# Extract all effects
all_effects <- bind_rows(
  # Summary models
  extract_effects(m1_ecological, "M1: Ecological", params_m1),
  extract_effects(m2_overall, "M2: Overall Rx", params_m2),
  
  # Stratified models (4 classes - no colistin)
  extract_effects(carb_models$acquisition, "M3A: Carb Acquisition", params_acquisition),
  extract_effects(carb_models$persistence, "M3B: Carb Persistence", params_persistence),
  extract_effects(ag_models$acquisition, "M4A: AG Acquisition", params_acquisition),
  extract_effects(ag_models$persistence, "M4B: AG Persistence", params_persistence),
  extract_effects(fq_models$acquisition, "M5A: FQ Acquisition", params_acquisition),
  extract_effects(fq_models$persistence, "M5B: FQ Persistence", params_persistence),
  extract_effects(ceph_models$acquisition, "M6A: Ceph Acquisition", params_acquisition),
  extract_effects(ceph_models$persistence, "M6B: Ceph Persistence", params_persistence),
  # NOTE: M7 (Colistin stratified) removed
  
  # Cross-resistance models (16 pairs - asymmetric for colistin)
  extract_effects(m7a_carb_to_ag, "M7A: Carb→AG", params_cross),
  extract_effects(m7b_carb_to_fq, "M7B: Carb→FQ", params_cross),
  extract_effects(m7c_carb_to_ceph, "M7C: Carb→Ceph", params_cross),
  extract_effects(m7d_ag_to_carb, "M7D: AG→Carb", params_cross),
  extract_effects(m7e_ag_to_fq, "M7E: AG→FQ", params_cross),
  extract_effects(m7f_ag_to_ceph, "M7F: AG→Ceph", params_cross),
  extract_effects(m7g_fq_to_carb, "M7G: FQ→Carb", params_cross),
  extract_effects(m7h_fq_to_ag, "M7H: FQ→AG", params_cross),
  extract_effects(m7i_fq_to_ceph, "M7I: FQ→Ceph", params_cross),
  extract_effects(m7j_ceph_to_carb, "M7J: Ceph→Carb", params_cross),
  extract_effects(m7k_ceph_to_ag, "M7K: Ceph→AG", params_cross),
  extract_effects(m7l_ceph_to_fq, "M7L: Ceph→FQ", params_cross),
  # Colistin as EXPOSURE only (4 models)
  extract_effects(m7m_coli_to_carb, "M7M: Coli→Carb", params_cross),
  extract_effects(m7n_coli_to_ag, "M7N: Coli→AG", params_cross),
  extract_effects(m7o_coli_to_fq, "M7O: Coli→FQ", params_cross),
  extract_effects(m7p_coli_to_ceph, "M7P: Coli→Ceph", params_cross)
)


# -----------------------------------------------------------------------------
# 1.5  Cross-resistance visualization
# 
# -----------------------------------------------------------------------------

# =============================================================================

# Explicit mapping for 16 cross-resistance models
cross_model_mapping <- tibble(
  Model = c(
    "M7A: Carb→AG", "M7B: Carb→FQ", "M7C: Carb→Ceph",
    "M7D: AG→Carb", "M7E: AG→FQ", "M7F: AG→Ceph",
    "M7G: FQ→Carb", "M7H: FQ→AG", "M7I: FQ→Ceph",
    "M7J: Ceph→Carb", "M7K: Ceph→AG", "M7L: Ceph→FQ",
    "M7M: Coli→Carb", "M7N: Coli→AG", "M7O: Coli→FQ", "M7P: Coli→Ceph"
  ),
  Cross_Exposure = c(
    rep("Carbapenem", 3), rep("Aminoglycoside", 3), rep("Fluoroquinolone", 3),
    rep("Cephalosporin", 3), rep("Colistin", 4)
  ),
  Outcome = c(
    "AG Resistance", "FQ Resistance", "Ceph Resistance",
    "Carb Resistance", "FQ Resistance", "Ceph Resistance",
    "Carb Resistance", "AG Resistance", "Ceph Resistance",
    "Carb Resistance", "AG Resistance", "FQ Resistance",
    "Carb Resistance", "AG Resistance", "FQ Resistance", "Ceph Resistance"
  ),
  Cross_Variable = c(
    rep("Carbapenem Rx", 3), rep("Aminoglycoside Rx", 3), rep("Fluoroquinolone Rx", 3),
    rep("Cephalosporin Rx", 3), rep("Colistin Rx", 4)
  )
)

# Extract cross-class effects only
cross_effects <- all_effects %>%
  filter(str_detect(Model, "M7")) %>%
  inner_join(cross_model_mapping, by = "Model") %>%
  filter(Variable == Cross_Variable)

cat("\nCross-effects extracted (should be 16):\n")
print(cross_effects %>% select(Model, Cross_Exposure, Outcome, Estimate, Significance))

# ASYMMETRIC Heatmap (5 exposures × 4 outcomes)
fig_cross_matrix <- cross_effects %>%
  mutate(
    Exposure_short = recode(Cross_Exposure,
                            "Carbapenem" = "Carb",
                            "Aminoglycoside" = "AG", 
                            "Fluoroquinolone" = "FQ",
                            "Cephalosporin" = "Ceph",
                            "Colistin" = "Coli"),
    Outcome_short = recode(Outcome,
                           "Carb Resistance" = "Carb R",
                           "AG Resistance" = "AG R",
                           "FQ Resistance" = "FQ R",
                           "Ceph Resistance" = "Ceph R"),
    label = paste0(sprintf("%.2f", Estimate), 
                   ifelse(Significance != "", paste0("\n", Significance), ""))
  ) %>%
  ggplot(aes(x = Exposure_short, y = Outcome_short)) +
  geom_tile(aes(fill = Estimate), color = "white", linewidth = 1.5) +
  geom_text(aes(label = label), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-1.5, 1.5),
    name = "Log-odds"
  ) +
  scale_x_discrete(limits = c("Carb", "AG", "FQ", "Ceph", "Coli")) +
  scale_y_discrete(limits = c("Carb R", "AG R", "FQ R", "Ceph R")) +
  labs(
    title = "Cross-Resistance Matrix (Asymmetric: 5 Exposures × 4 Outcomes)",
    subtitle = "Colistin as exposure only (insufficient outcome data)\n*** p>97.5% | ** p>95% | * p>90%",
    x = "Prior Drug Exposure",
    y = "Resistance Outcome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank()
  )

ggsave("figures/glmm/cross_resistance_matrix_asymmetric.pdf", fig_cross_matrix,
       width = 10, height = 7, dpi = 300)


# -----------------------------------------------------------------------------
# 1.6  Colistin-specific analysis
# 
# -----------------------------------------------------------------------------

# =============================================================================

colistin_effects <- cross_effects %>%
  filter(Cross_Exposure == "Colistin") %>%
  select(Outcome, Estimate, CI_lower, CI_upper, OR, Significance, Direction)

print(colistin_effects)

# Forest plot for colistin effects
fig_coli_effects <- colistin_effects %>%
  mutate(Outcome = fct_reorder(Outcome, Estimate)) %>%
  ggplot(aes(x = Estimate, y = Outcome)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2, color = "grey40") +
  geom_point(aes(color = Direction), size = 4) +
  geom_text(aes(label = Significance), hjust = -0.5, size = 5) +
  scale_color_manual(values = c("↑ Resistance" = "#E41A1C", 
                                "↓ Resistance" = "#377EB8",
                                "Uncertain" = "grey50")) +
  labs(
    title = "Colistin Exposure -> Resistance to Other Drug Classes",
    subtitle = "Does colistin use promote cross-resistance?",
    x = "Log-odds (effect on resistance)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("figures/glmm/forest_colistin_cross_effects.pdf", fig_coli_effects,
       width = 8, height = 4, dpi = 300)

write_csv(colistin_effects, "results/colistin_cross_resistance.csv")


# -----------------------------------------------------------------------------
# 1.7  Model summary and forest plots
# 
# -----------------------------------------------------------------------------


model_summary <- tibble(
  `Model Group` = c(
    "M1: Ecological",
    "M2: Overall Rx",
    "M3-M6: Stratified (4 classes)",
    "M7: Cross-resistance (16 pairs)"
  ),
  `N Models` = c(1, 1, 8, 16),
  `Colistin Role` = c(
    "Not included",
    "Predictor",
    "Predictor only (no outcome)",
    "Exposure only (4 models)"
  ),
  Notes = c(
    "Ecological factors only",
    "All 5 drugs as predictors",
    "Carb, AG, FQ, Ceph outcomes; Colistin outcome excluded",
    "5 exposures × 4 outcomes (excluding Coli→Coli)"
  )
)

print(model_summary)

# Save all effects
write_csv(all_effects, "results/glmm_all_effects_asymmetric.csv")

significant_effects <- all_effects %>%
  filter(Significance != "") %>%
  select(Model, Variable, Estimate, CI_lower, CI_upper, OR, Significance, Direction) %>%
  arrange(Model, desc(abs(Estimate)))

write_csv(significant_effects, "results/glmm_significant_effects_asymmetric.csv")

write_csv(cross_effects %>%
            select(Cross_Exposure, Outcome, Estimate, CI_lower, CI_upper, Significance, Direction),
          "results/cross_resistance_summary_asymmetric.csv")


#===========================================================================
# FOREST PLOTS FOR SUMMARY MODELS (M1, M2)
# =============================================================================

# -----------------------------------------------------------------------------
# Forest plot function
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# M1: Ecological model forest plot
# -----------------------------------------------------------------------------

fig_m1 <- create_forest_plot(
  
  all_effects %>% filter(Model == "M1: Ecological"),
  "M1: Ecological Factors",
  paste0("Beta-binomial GLMM | N = ", nrow(pa_ecological), " visits")
)

ggsave("figures/glmm/forest_m1_ecological.pdf", fig_m1, 
       width = 8, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# M2: Overall Rx model forest plot
# -----------------------------------------------------------------------------

fig_m2 <- create_forest_plot(
  all_effects %>% filter(Model == "M2: Overall Rx"),
  "M2: Overall Antibiotic Exposure",
  paste0("Beta-binomial GLMM | N = ", nrow(pa_prescription), " visits with Rx data")
)

ggsave("figures/glmm/forest_m2_overall.pdf", fig_m2, 
       width = 8, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Combined M1 + M2 figure
# -----------------------------------------------------------------------------


fig_summary_combined <- fig_m1 + fig_m2 +
  
  plot_annotation(
    title = "Summary Models: Predictors of Overall PA Resistance",
    subtitle = "Ecological factors (M1) vs. with antibiotic exposures (M2)",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "grey30")
    )
  )

ggsave("figures/glmm/forest_summary_m1_m2.pdf", fig_summary_combined, 
       width = 14, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# Two companion models to M2:
#   m2_drugs_only  — strips ecological covariates; useful as a contrast to see
#                    how much the drug coefficients change when MRSA/MSSA/
#                    mucoid/exacerbation are removed (confounding check).
#   m2_enhanced    — adds temporal/utilization covariates on top of M2 to
#                    check robustness: days since previous visit, visit
#                    number, and prior count of drug classes prescribed.
# Note: prev_total_antibiotic_pressure is intentionally OMITTED because it is
# the sum of the 5 per-drug exposure scores and is therefore collinear with
# the individual prev_<drug>_prescribed terms already in the model.
# n_isolates is OMITTED because it is measured at time t (testing intensity
# at the current visit) and would violate the lagged-predictor design.
# -----------------------------------------------------------------------------

# ---- m2_drugs_only: M2 stripped of ecological + clinical covariates ---------
m2_drugs_only <- brm(
  n_resistant | trials(n_tested) ~ 
    prev_any_resistance + 
    prev_carbapenem_prescribed +
    prev_aminoglycoside_prescribed +
    prev_fluoroquinolones_prescribed +
    prev_cephalosporin_prescribed +
    prev_colistin_prescribed +
    (1 | individual_id),
  family = beta_binomial(),
  prior = priors_beta_random,
  data = pa_prescription,
  chains = 4, cores = 4,
  iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95),
  seed = 42,
  file = "models/m2_drugs_only",
  file_refit = "on_change"
)

# ---- m2_enhanced: M2 + temporal/utilization covariates ----------------------
# Does NOT include prev_has_exacerbation (constant in pa_prescription).
m2_enhanced <- brm(
  n_resistant | trials(n_tested) ~ 
    prev_any_resistance +
    prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_carbapenem_prescribed +
    prev_aminoglycoside_prescribed +
    prev_fluoroquinolones_prescribed +
    prev_cephalosporin_prescribed +
    prev_colistin_prescribed +
    scale(days_since_prev_visit) +
    scale(visit_number) +
    scale(prev_n_drug_classes) +
    (1 | individual_id),
  family = beta_binomial(),
  prior = priors_beta_random,
  data = pa_prescription %>% filter(!is.na(days_since_prev_visit),
                                    !is.na(prev_n_drug_classes)),
  chains = 4, cores = 4,
  iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95),
  seed = 42,
  file = "models/m2_enhanced",
  file_refit = "on_change"
)

# -----------------------------------------------------------------------------
# Helper functions — these are ALREADY defined earlier in this file
# (extract_effects, create_forest_plot). The duplicate definitions that used
# to live here have been removed to avoid silent drift if one copy is updated.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Parameter map for m2_enhanced
# brms strips non-alphanumerics from parameter names, so scale(x) appears as
# b_scalex. Both raw and scaled variants are listed so the map is robust to
# either coding.
# -----------------------------------------------------------------------------

params_m2_enhanced <- c(
  # Shared with m2_overall
  "b_prev_any_resistance"               = "Prior Resistance",
  "b_prev_has_mrsa"                     = "MRSA Co-colonization",
  "b_prev_has_mssa"                     = "MSSA Co-colonization",
  "b_prev_has_any_mucoid"               = "Mucoid Phenotype",
  "b_prev_carbapenem_prescribed"        = "Carbapenem Rx",
  "b_prev_aminoglycoside_prescribed"    = "Aminoglycoside Rx",
  "b_prev_fluoroquinolones_prescribed"  = "Fluoroquinolone Rx",
  "b_prev_cephalosporin_prescribed"     = "Cephalosporin Rx",
  "b_prev_colistin_prescribed"          = "Colistin Rx",
  
  # New in m2_enhanced (scaled continuous covariates)
  "b_scaledays_since_prev_visit"        = "Days Since Prev Visit (z)",
  "b_scalevisit_number"                 = "Visit Number (z)",
  "b_scaleprev_n_drug_classes"          = "N Drug Classes, prior (z)"
)

# -----------------------------------------------------------------------------
# Extract effects + print sanity check
# -----------------------------------------------------------------------------

effects_m2_drugs_only <- extract_effects(
  m2_drugs_only, "M2 drugs only", params_m2
)

effects_m2_enhanced <- extract_effects(
  m2_enhanced, "M2 enhanced", params_m2_enhanced
)

cat("Extracted effects (m2_enhanced):\n")
effects_m2_enhanced %>%
  arrange(desc(abs(Estimate))) %>%
  dplyr::select(Variable, Estimate, CI_lower, CI_upper, OR, Prob_positive, Significance) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  print(n = 30)

# Flag any b_* params in the fitted model that aren't in the map — so we
# don't silently drop a variable from the forest plot.
all_b_params_enh <- as_draws_df(m2_enhanced) %>%
  names() %>%
  keep(~ startsWith(.x, "b_") & .x != "b_Intercept")

missing_from_map <- setdiff(all_b_params_enh, names(params_m2_enhanced))
if (length(missing_from_map) > 0) {
  cat("\n⚠  Parameters in m2_enhanced but not in params_m2_enhanced:\n")
  cat(paste("  ", missing_from_map, collapse = "\n"), "\n")
}

# -----------------------------------------------------------------------------
# Forest plot
# -----------------------------------------------------------------------------

fig_m2_enhanced <- create_forest_plot(
  effects_m2_enhanced,
  "Effect of Antibiotic Exposure on Any PA Resistance (M2 enhanced)",
  paste0("Beta-binomial GLMM | N = ", nobs(m2_enhanced), " visits")
)

print(fig_m2_enhanced)

ggsave("figures/glmm/forest_m2_enhanced.pdf", fig_m2_enhanced,
       width = 9, height = 7, dpi = 300)
ggsave("figures/glmm/forest_m2_enhanced.pdf", fig_m2_enhanced,
       width = 9, height = 7)

# Append m2_drugs_only and m2_enhanced into all_effects so downstream
# consumers (causal-network figure, comparison CSV) can see them.
all_effects <- bind_rows(
  all_effects,
  effects_m2_drugs_only,
  effects_m2_enhanced
)



# -----------------------------------------------------------------------------
# 1.8  Model diagnostics
# 
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Diagnostic extraction function
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Extract diagnostics for all models
# -----------------------------------------------------------------------------

diagnostics_list <- list(
  # Summary models
  extract_diagnostics(m1_ecological, "M1: Ecological"),
  extract_diagnostics(m2_overall, "M2: Overall Rx"),
  
  # Stratified models
  extract_diagnostics(carb_models$acquisition, "M3A: Carb Acquisition"),
  extract_diagnostics(carb_models$persistence, "M3B: Carb Persistence"),
  extract_diagnostics(ag_models$acquisition, "M4A: AG Acquisition"),
  extract_diagnostics(ag_models$persistence, "M4B: AG Persistence"),
  extract_diagnostics(fq_models$acquisition, "M5A: FQ Acquisition"),
  extract_diagnostics(fq_models$persistence, "M5B: FQ Persistence"),
  extract_diagnostics(ceph_models$acquisition, "M6A: Ceph Acquisition"),
  extract_diagnostics(ceph_models$persistence, "M6B: Ceph Persistence"),
  
  # Cross-resistance models
  extract_diagnostics(m7a_carb_to_ag, "M7A: Carb→AG"),
  extract_diagnostics(m7b_carb_to_fq, "M7B: Carb→FQ"),
  extract_diagnostics(m7c_carb_to_ceph, "M7C: Carb→Ceph"),
  extract_diagnostics(m7d_ag_to_carb, "M7D: AG→Carb"),
  extract_diagnostics(m7e_ag_to_fq, "M7E: AG→FQ"),
  extract_diagnostics(m7f_ag_to_ceph, "M7F: AG→Ceph"),
  extract_diagnostics(m7g_fq_to_carb, "M7G: FQ→Carb"),
  extract_diagnostics(m7h_fq_to_ag, "M7H: FQ→AG"),
  extract_diagnostics(m7i_fq_to_ceph, "M7I: FQ→Ceph"),
  extract_diagnostics(m7j_ceph_to_carb, "M7J: Ceph→Carb"),
  extract_diagnostics(m7k_ceph_to_ag, "M7K: Ceph→AG"),
  extract_diagnostics(m7l_ceph_to_fq, "M7L: Ceph→FQ"),
  extract_diagnostics(m7m_coli_to_carb, "M7M: Coli→Carb"),
  extract_diagnostics(m7n_coli_to_ag, "M7N: Coli→AG"),
  extract_diagnostics(m7o_coli_to_fq, "M7O: Coli→FQ"),
  extract_diagnostics(m7p_coli_to_ceph, "M7P: Coli→Ceph")
)

# Combine
model_diagnostics <- bind_rows(diagnostics_list)

# -----------------------------------------------------------------------------
# Print diagnostic summary
# -----------------------------------------------------------------------------

# Summary table
cat("CONVERGENCE & SAMPLING:\n")
print(
  model_diagnostics %>%
    select(Model, N_obs, Rhat_max, ESS_bulk_min, N_divergent, Quality) %>%
    as.data.frame(),
  row.names = FALSE
)

# Count issues
n_good <- sum(model_diagnostics$Quality == "✓ Good", na.rm = TRUE)
n_total <- nrow(model_diagnostics)

cat("\n")
cat("OVERALL: ", n_good, "/", n_total, " models passed all diagnostics\n", sep = "")

# Flag problematic models
problematic <- model_diagnostics %>% filter(Quality != "✓ Good")
if (nrow(problematic) > 0) {
  cat("\nMODELS REQUIRING ATTENTION:\n")
  print(problematic %>% select(Model, Rhat_max, ESS_bulk_min, N_divergent, Quality))
}

# Save diagnostics
write_csv(model_diagnostics, "results/model_diagnostics.csv")

# -----------------------------------------------------------------------------
# Diagnostic visualization
# -----------------------------------------------------------------------------

# R-hat plot
fig_rhat <- model_diagnostics %>%
  mutate(
    Model_short = str_remove(Model, ": .*"),
    Model = fct_reorder(Model, Rhat_max)
  ) %>%
  ggplot(aes(x = Rhat_max, y = Model, fill = Rhat_ok)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 1.01, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C"),
                    labels = c("TRUE" = "OK", "FALSE" = "Issue"),
                    name = "R-hat < 1.01") +
  labs(
    title = "Convergence Diagnostic: R-hat",
    subtitle = "Values should be < 1.01 (red line)",
    x = "Maximum R-hat across parameters"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ESS plot
fig_ess <- model_diagnostics %>%
  mutate(Model = fct_reorder(Model, ESS_bulk_min)) %>%
  ggplot(aes(x = ESS_bulk_min, y = Model, fill = ESS_ok)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 400, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C"),
                    labels = c("TRUE" = "OK", "FALSE" = "Low"),
                    name = "ESS > 400") +
  labs(
    title = "Effective Sample Size (Bulk ESS)",
    subtitle = "Values should be > 400 (red line)",
    x = "Minimum Bulk ESS across parameters"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Combined diagnostic figure
fig_diagnostics_combined <- fig_rhat + fig_ess +
  plot_annotation(
    title = "Model Quality Diagnostics",
    tag_levels = "A",
    theme = theme(plot.title = element_text(face = "bold", size = 16))
  )

ggsave("figures/glmm/model_diagnostics.pdf", fig_diagnostics_combined, 
       width = 14, height = 10, dpi = 300)

# -----------------------------------------------------------------------------
# Posterior predictive checks for key models
# -----------------------------------------------------------------------------

# M1: Ecological
pp_m1 <- pp_check(m1_ecological, type = "dens_overlay", ndraws = 100) +
  labs(title = "M1: Ecological", subtitle = "Posterior predictive check") +
  theme_minimal()

# M2: Overall
pp_m2 <- pp_check(m2_overall, type = "dens_overlay", ndraws = 100) +
  labs(title = "M2: Overall Rx", subtitle = "Posterior predictive check") +
  theme_minimal()

# Combine PP checks
fig_pp_checks <- pp_m1 + pp_m2 +
  plot_annotation(
    title = "Posterior Predictive Checks: Summary Models",
    subtitle = "Dark line = observed data | Light lines = posterior predictions",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave("figures/glmm/posterior_predictive_checks.pdf", fig_pp_checks, 
       width = 12, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# Trace plots for M1 and M2 (key parameters)
# -----------------------------------------------------------------------------

# M1 trace plots
fig_trace_m1 <- mcmc_trace(m1_ecological, 
                           pars = c("b_prev_any_resistance", 
                                    "b_prev_has_mrsa", 
                                    "b_prev_has_mssa",
                                    "b_prev_has_any_mucoid")) +
  labs(title = "M1: Ecological - MCMC Trace Plots") +
  theme_minimal()

ggsave("figures/glmm/trace_m1_ecological.pdf", fig_trace_m1, 
       width = 10, height = 8, dpi = 300)

# M2 trace plots (drug effects)
fig_trace_m2 <- mcmc_trace(m2_overall, 
                           pars = c("b_prev_carbapenem_prescribed",
                                    "b_prev_aminoglycoside_prescribed",
                                    "b_prev_fluoroquinolones_prescribed",
                                    "b_prev_cephalosporin_prescribed",
                                    "b_prev_colistin_prescribed")) +
  labs(title = "M2: Overall Rx - MCMC Trace Plots (Drug Effects)") +
  theme_minimal()

ggsave("figures/glmm/trace_m2_drug_effects.pdf", fig_trace_m2, 
       width = 10, height = 8, dpi = 300)


check_m2_drugs_only <- extract_diagnostics(m2_drugs_only, "m2_drugs_only")
check_m2_enhanced   <- extract_diagnostics(m2_enhanced,   "m2_enhanced")


# #############################################################################
# SECTION 2: EXPOSURE NOISE  -  prescribed != taken
# Models fit once here; the analysis and talk figures below both read them.
# #############################################################################


# -----------------------------------------------------------------------------
# 2.0  Conservative-coding guard
# 
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# 7.0 CONSERVATIVE-CODING GUARD
# These models inherit their coding from the CSVs. Conservative (R+I) requires
# the v2 prep. Warn loudly if the data look like the old (strict) prep.
# -----------------------------------------------------------------------------
if (!"n_resistant_strict" %in% names(pa_visit_full)) {
  warning(
    "pa_visit_full has no 'n_resistant_strict' column -> looks like PRE-v2 prep.\n",
    "  Resistance counts may be STRICT (R only), not conservative (R+I).\n",
    "  Re-run analytical_data_prep.R and reload the CSVs before fitting."
  )
} else {
  cat("v2 prep detected: n_resistant is conservative (R + I). Good.\n")
}


# -----------------------------------------------------------------------------
# 2.1  Route intensity vs binary prescribed  [Slide 4]
# 
# -----------------------------------------------------------------------------

# Question: does route-weighted exposure intensity (0-3) explain resistance
# better than mere prescribed-yes/no, and does it shrink between-patient
# heterogeneity? Both models are fit on the SAME rows so LOO is comparable.
#
# Caveats baked in:
#   - Route can only touch BETWEEN-patient / BETWEEN-visit variance, NOT the
#     within-sample (between-isolate) slice from 7.1 - that needs deep sampling.
#   - "IV = supervised dose" is imperfect (home IV / OPAT is common in CF).
#   - Colistin keeps its binary term (no lagged exposure score exists); it is
#     identical in both models, so the comparison isolates the 4-class swap.
# =============================================================================

route_predictors <- c("prev_aminoglycoside_exposure", "prev_carbapenem_exposure",
                      "prev_cephalosporin_exposure", "prev_fluoroquinolone_exposure")
bin_predictors   <- c("prev_aminoglycoside_prescribed", "prev_carbapenem_prescribed",
                      "prev_cephalosporin_prescribed", "prev_fluoroquinolones_prescribed")
shared_predictors <- c("prev_any_resistance", "prev_has_mrsa", "prev_has_mssa",
                       "prev_has_any_mucoid", "prev_colistin_prescribed")

route_df <- pa_prescription %>%
  filter(if_all(all_of(c("n_resistant", "n_tested", "individual_id",
                         shared_predictors, route_predictors, bin_predictors)),
                ~ !is.na(.))) %>%
  # z-score the exposure scores so the normal(0,1.5) prior is sensible and the
  # route coefficients are per-SD (comparable to the binary 0/1 coefficients)
  mutate(across(all_of(route_predictors), ~ as.numeric(scale(.)), .names = "{.col}_z"))

cat("\n7.2 route comparison frame: N =", nrow(route_df), "visits\n")

# Binary M2 refit on the common frame (apples-to-apples with m2_route)
m2_binary <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance + prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_aminoglycoside_prescribed + prev_carbapenem_prescribed +
    prev_cephalosporin_prescribed + prev_fluoroquinolones_prescribed +
    prev_colistin_prescribed + (1 | individual_id),
  data = route_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/m2_binary_common", file_refit = "on_change"
)

m2_route <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance + prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_aminoglycoside_exposure_z + prev_carbapenem_exposure_z +
    prev_cephalosporin_exposure_z + prev_fluoroquinolone_exposure_z +
    prev_colistin_prescribed + (1 | individual_id),
  data = route_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/m2_route", file_refit = "on_change"
)

# (i) Predictive comparison
loo_route_cmp <- loo_compare(loo(m2_binary), loo(m2_route))
cat("\nLOO: binary vs route (more negative elpd_diff = worse than best)\n")
print(loo_route_cmp)

# (ii) Does route absorb between-patient heterogeneity?
sd_bin   <- as_draws_df(m2_binary)$sd_individual_id__Intercept
sd_route <- as_draws_df(m2_route)$sd_individual_id__Intercept
patient_sd_compare <- tibble(
  model = c("binary", "route"),
  patient_sd_median = c(median(sd_bin), median(sd_route)),
  lo = c(quantile(sd_bin, .025), quantile(sd_route, .025)),
  hi = c(quantile(sd_bin, .975), quantile(sd_route, .975))
)
cat("\nBetween-patient SD - if route < binary, route resolves some heterogeneity:\n")
print(patient_sd_compare)
write_csv(patient_sd_compare, "results/noise/route_patient_sd_compare.csv")


# -----------------------------------------------------------------------------
# 2.2  Talk figures: main forest and route comparison
# From Talk_figures_v2 Section 3. Consumes m2_binary / m2_route above.
# -----------------------------------------------------------------------------

# MAIN FOREST: the standard analysis, all predictors in one model  [Slide 4]
# m2_binary already carries prior resistance + ecology + drug exposure. Order
# drugs first, then prior state & ecology; shade by whether the 95% CI excludes
# zero. NOTE: conditioned on exacerbation visits (drugs only recorded then), so
# a drug effect here partly reflects how sick the patient was - see Section 3B.
# -----------------------------------------------------------------------------
drug_terms <- c("Aminoglycoside", "Carbapenem", "Cephalosporin",
                "Fluoroquinolone", "Colistin")
main_recodes <- c(prev_any_resistance              = "Prior resistance",
                  prev_has_mrsa                    = "MRSA",
                  prev_has_mssa                    = "MSSA",
                  prev_has_any_mucoid              = "Mucoid",
                  prev_aminoglycoside_prescribed   = "Aminoglycoside",
                  prev_carbapenem_prescribed       = "Carbapenem",
                  prev_cephalosporin_prescribed    = "Cephalosporin",
                  prev_fluoroquinolones_prescribed = "Fluoroquinolone",
                  prev_colistin_prescribed         = "Colistin")

# Full posteriors (not point estimates) - same halfeye idiom as fig2/3/4, so
# the audience is in the distribution mindset from the first model slide.
# Significance = probability of direction > 0.975 (solid vs faded), same rule
# used throughout. Self-contained: does not need the Section 4 draw helpers.
main_draws <- as_draws_df(m2_binary) %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  mutate(variable = str_remove(param, "^b_")) %>%
  filter(variable != "Intercept") %>%
  mutate(pretty_var = recode(variable, !!!main_recodes)) %>%
  group_by(pretty_var) %>%
  mutate(pd          = max(mean(value > 0), mean(value < 0)),
         direction   = if_else(median(value) > 0, "+", "-"),
         significant = pd > 0.975) %>%
  ungroup() %>%
  mutate(group = if_else(pretty_var %in% drug_terms,
                         "Prior state & ecology", "Antibiotic exposure (prior)"),
         group = factor(group, levels = c("Prior state & ecology",
                                          "Antibiotic exposure (prior)"
         )),
         pretty_var = fct_reorder(pretty_var, value, .fun = median))

fig_main_forest <- ggplot(main_draws, aes(value, pretty_var)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  stat_halfeye(aes(fill = direction, alpha = significant),
               .width = c(0.66, 0.95), slab_color = NA, point_size = 2,
               interval_size_range = c(0.6, 1.6)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("+" = pal[["red"]], "-" = pal[["blue"]]), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.35), guide = "none") +
  labs(x = "Effect on resistance (log-odds)", y = NULL,
       title = "What a standard analysis says",
       caption = lab_n(route_df)) +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        plot.title = element_text(face = "bold", size = 16))

ggsave("figures/talk/fig_main_forest.pdf", fig_main_forest,
       width = 9, height = 6.5, dpi = 300, bg = "transparent")

route_forest_df <- bind_rows(
  get_effects(m2_binary, "Binary (prescribed y/n)"),
  get_effects(m2_route,  "Route intensity")
) %>%
  mutate(term = fct_reorder(term, est))

fig_route <- ggplot(route_forest_df,
                    aes(est, term, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.5), linewidth = 0.8) +
  scale_color_manual(values = c("Binary (prescribed y/n)" = "#AB9968",
                                "Route intensity" = "#207E9D"), name = NULL) +
  labs(x = "Effect on resistance (log-odds)", y = NULL,
       title = "Does route intensity sharpen the drug signal?",
       subtitle = "Same patients, same rows; binary vs route-weighted exposure",
       caption = lab_n(route_df)) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16))

ggsave("figures/talk/fig_route_vs_binary.pdf", fig_route, width = 9, height = 6, dpi = 300, bg = "transparent")


# =============================================================================

# -----------------------------------------------------------------------------
# 2.3  Evidence of confounding by indication
# 
# -----------------------------------------------------------------------------

# Prescriptions are logged almost only at exacerbations (99% of APE-IV visits,
# 0% when stable), and the strong drugs (carbapenem, colistin, IV) concentrate
# in the sickest visits. So a "drug prescribed" predictor partly encodes HOW
# SICK the patient was. We do not adjust this away - it is collinear with
# clinical status by design (the prev_has_exacerbation = 1 problem). We SHOW it
# and explain it. On-theme: the sample is not random; it is driven by how sick
# the patient is when someone decides to culture and treat.
# -----------------------------------------------------------------------------

# (1) Treatment responds to prior resistance (indication / feedback) ----------
# Avoidance of failing first-line drugs, escalation to last-resort. The point is
# that exposure and prior resistance are entangled, not that direction is fixed.

feedback <- bind_rows(
  rx_by_prior(pa_visit_full, "aminoglycoside_prescribed",   "prev_ag_any_resistant",   "Aminoglycoside"),
  rx_by_prior(pa_visit_full, "carbapenem_prescribed",       "prev_carb_any_resistant", "Carbapenem"),
  rx_by_prior(pa_visit_full, "fluoroquinolones_prescribed", "prev_fq_any_resistant",   "Fluoroquinolone"),
  rx_by_prior(pa_visit_full, "cephalosporin_prescribed",    "prev_ceph_any_resistant", "Cephalosporin")
) %>%
  mutate(prior = factor(prior, levels = c("Prior S", "Prior R")))

fig_conf_rx <- ggplot(feedback, aes(class, rate, fill = prior)) +
  geom_col(position = position_dodge(0.7), width = 0.62) +
  scale_fill_manual(values = c("Prior S" = pal[["blue"]], "Prior R" = pal[["red"]])) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.05))) +
  coord_flip() +
  labs(x = NULL, y = "Prescribed this class this visit",
       title = "Treatment responds to prior resistance",
       caption = lab_n(pa_prescription)) +
  theme(legend.position = "top")

ggsave("figures/talk/fig_confound_rx_by_prior.pdf", fig_conf_rx,
       width = 9, height = 6, dpi = 300, bg = "transparent")

# (2) Strong drugs concentrate in the sickest visits (severity) ---------------
# Prescription is only recorded at acute visits, so at stable visits we treat a
# missing record as "not prescribed" (rate contribution 0) - clinically these
# patients are not on acute antipseudomonal courses.
sev_long <- pa_visit_full %>%
  filter(!is.na(clinical_status), clinical_status != "") %>%
  mutate(status = factor(clinical_status,
                         levels = c("Clinically Stable", "APE Outpatient", "APE IV Antibiotics"),
                         labels = c("Stable", "APE outpt", "APE IV"))) %>%
  transmute(status,
            Aminoglycoside  = coalesce(aminoglycoside_prescribed, 0),
            Fluoroquinolone = coalesce(fluoroquinolones_prescribed, 0),
            Carbapenem      = coalesce(carbapenem_prescribed, 0),
            Colistin        = coalesce(colistin_prescribed, 0)) %>%
  filter(!is.na(status)) %>%
  pivot_longer(-status, names_to = "drug", values_to = "v") %>%
  group_by(status, drug) %>%
  summarise(rate = mean(v > 0), .groups = "drop") %>%
  mutate(drug = factor(drug, levels = c("Aminoglycoside", "Fluoroquinolone",
                                        "Carbapenem", "Colistin")))

fig_conf_sev <- ggplot(sev_long, aes(status, rate, fill = drug)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  scale_fill_manual(values = c("Aminoglycoside" = pal[["blue"]],
                               "Fluoroquinolone" = pal[["gold"]],
                               "Carbapenem" = pal[["red"]],
                               "Colistin" = pal[["plum"]])) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.05))) +
  labs(x = NULL, y = "Visits prescribed",
       title = "The strong drugs are markers of how sick the patient is",
       caption = lab_n(filter(pa_visit_full, !is.na(clinical_status)))) +
  theme(legend.position = "top")

ggsave("figures/talk/fig_confound_severity.pdf", fig_conf_sev,
       width = 9, height = 6, dpi = 300, bg = "transparent")

# =============================================================================

# #############################################################################
# SECTION 3: SAMPLING NOISE  -  few colonies != the population
# The dominant noise source: where the variance actually lives, and how
# precision scales with the number of colonies tested.
# #############################################################################


# -----------------------------------------------------------------------------
# 3.1  Variance decomposition, AG + carbapenem  [Slide 3]
# Fits vardecomp_ag / vardecomp_carb and binds them into vardecomp_2class.
# -----------------------------------------------------------------------------

# Isolate-level logistic with patient and visit random effects. Conservative
# (R + I) resistance call. Three components:
#   between-patient                 = var(individual_id)
#   between-visit, within-patient   = var(individual_id:visit_number)
#   between-isolate, within-visit   = residual (pi^2/3, latent scale)
# The last is what single-colony sampling discards.
# =============================================================================


vardecomp_ag   <- fit_variance_decomp(pa_isolates,
                                      c("amikacin", "tobramycin", "gentamicin"),
                                      "Aminoglycoside", "ag")
vardecomp_carb <- fit_variance_decomp(pa_isolates,
                                      c("meropenem", "imipenem"),
                                      "Carbapenem", "carb")

vardecomp_2class <- bind_rows(vardecomp_ag, vardecomp_carb)
write_csv(vardecomp_2class, "results/noise/variance_decomposition_2class.csv")
print(vardecomp_2class)


# -----------------------------------------------------------------------------
# 3.2  Variance decomposition, 4 well-sampled classes
# Adds FQ and cephalosporin. The AG and carbapenem fits that opened this
# block in Talk_figures_v2 are dropped - 3.1 above already ran them with
# identical arguments and cache tags. Both the 4-class and the 2-class
# figures are drawn here; the 2-class one needs vardecomp_2class from 3.1.
# -----------------------------------------------------------------------------

# 2A. Where the variance lives (AG + carbapenem)  [Slide 3]
# 7.1 VARIANCE DECOMPOSITION - AG + CARBAPENEM  [Slide 3]
# Isolate-level logistic with patient and visit random effects. Conservative
# (R + I) resistance call. Three components:
#   between-patient                 = var(individual_id)
#   between-visit, within-patient   = var(individual_id:visit_number)
#   between-isolate, within-visit   = residual (pi^2/3, latent scale)
# The last is what single-colony sampling discards.
# =============================================================================


vardecomp_fq   <- fit_variance_decomp(pa_isolates,
                                      c("ciprofloxacin", "levofloxacin"),
                                      "Fluoroquinolone", "fq")
vardecomp_ceph <- fit_variance_decomp(pa_isolates,
                                      c("cefepime", "ceftazidime"),
                                      "Cephalosporin", "ceph")

vardecomp_4class <- bind_rows(vardecomp_ag, vardecomp_carb,
                              vardecomp_fq, vardecomp_ceph) %>%
  mutate(class = factor(class, levels = c("Aminoglycoside", "Carbapenem",
                                          "Fluoroquinolone", "Cephalosporin")))
write_csv(vardecomp_4class, "results/noise/variance_decomposition_4class.csv")
print(vardecomp_4class)

# =============================================================================

# 2A figure: variance decomposition bars
# --- 3.3 Variance decomposition, AG + carbapenem [Slide 3] -------------------
# 4 well-sampled classes (AG, carbapenem, FQ, cephalosporin)
vd <- vardecomp_4class %>%
  mutate(component = recode(component,
                            "between_patient" = "Between patients",
                            "between_visit_within_patient" = "Between visits\n(within patient)",
                            "between_isolate_within_visit" = "Between isolates\n(within visit)"),
         component = factor(component, levels = c("Between patients",
                                                  "Between isolates\n(within visit)",
                                                  "Between visits\n(within patient)")))

fig_variance <- ggplot(vd, aes(component, median, fill = class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = scales::percent(median, accuracy = 1)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 6) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.12))) +
  scale_fill_manual(values = c("Aminoglycoside"  = pal[["blue"]],
                               "Carbapenem"      = pal[["red"]],
                               "Fluoroquinolone" = pal[["gold"]],
                               "Cephalosporin"   = pal[["plum"]]),
                    name = NULL) +
  labs(x = NULL, y = "Share of total variance",
       #title = "Where does resistance variance live?",
       #subtitle = "The within-visit slice is what single-colony sampling discards",
       caption = lab_n(pa_isolates)) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16))

ggsave("figures/talk/fig_variance_4class.pdf", fig_variance, width = 9, height = 6, dpi = 300, bg = "transparent")

# =============================================================================

# 2A figure: variance decomposition bars
# --- 3.3 Variance decomposition, AG + carbapenem [Slide 3] -------------------
# vardecomp_2class comes from Synthesis.R Section 7.1
vd <- vardecomp_2class %>%
  mutate(component = recode(component,
                            "between_patient" = "Between patients",
                            "between_visit_within_patient" = "Between visits\n(within patient)",
                            "between_isolate_within_visit" = "Between isolates\n(within visit)"),
         component = factor(component, levels = c("Between patients",
                                                  "Between isolates\n(within visit)",
                                                  "Between visits\n(within patient)")))

fig_variance <- ggplot(vd, aes(component, median, fill = class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = scales::percent(median, accuracy = 0.1)),
            position = position_dodge(width = 0.8), vjust = -0.4, size = 4) +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.12))) +
  scale_fill_manual(values = c("Aminoglycoside" = "#207E9D", "Carbapenem" = "#90606C"),
                    name = NULL) +
  labs(x = NULL, y = "Share of total variance",
       title = "Where does resistance variance live?",
       subtitle = "The within-visit slice is what single-colony sampling discards",
       caption = lab_n(pa_isolates)) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16))

ggsave("figures/talk/fig_variance_2class.pdf", fig_variance, width = 9, height = 6, dpi = 300, bg = "transparent")


# 2B. One culture barely measures the within-host frequency  [Slide 5]
# 6A. PRECISION vs SAMPLING DEPTH  [TALK]
#
# How precisely can a single culture pin down the within-host resistance
# frequency, as a function of how many colonies you test?
#
# Model the per-visit count as k resistant of n tested, latent frequency theta
# with a Beta(1,1) prior. Posterior is Beta(1+k, 1+n-k); its SD measures how
# tightly the data constrain theta. We average that SD over plausible k for a
# given true frequency p (k ~ Binomial(n, p)) to get the EXPECTED precision at
# depth n. p = 0.5 is the hardest (widest) case; p = 0.1 / 0.9 are easier.
#
# The payoff: overlay the depth distribution that surveillance ACTUALLY uses
# (n_isolates = colonies per culture) to show where current practice sits.
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# 3.3  Precision vs sampling depth  [Slide 5]
# Restored from Talk_figures.R: the empirical-depth bottom panel (fig6a_bot)
# was computed but never saved in Talk_figures_v2. Both the stacked png and
# the slide-sized top-panel pdf are written.
# -----------------------------------------------------------------------------

# Section 5 shows noise EXISTS and that ignoring it changes inferences.
# Section 6 turns that into (6A) a surveillance design ask and (6B) a
# demonstration that the three methods of Sections 2–4 are not equally robust
# to sampling noise.
#
# 6A is analytic and cheap; it always runs.
# 6B refits models on thinned data and is slower / newly written; it is fenced
#    behind RUN_SECTION_6B and should be run interactively and sanity-checked
#    before it goes anywhere near a slide.
# =============================================================================
# =============================================================================

# -----------------------------------------------------------------------------
# 6A. PRECISION vs SAMPLING DEPTH  [TALK]
#
# How precisely can a single culture pin down the within-host resistance
# frequency, as a function of how many colonies you test?
#
# Model the per-visit count as k resistant of n tested, latent frequency theta
# with a Beta(1,1) prior. Posterior is Beta(1+k, 1+n-k); its SD measures how
# tightly the data constrain theta. We average that SD over plausible k for a
# given true frequency p (k ~ Binomial(n, p)) to get the EXPECTED precision at
# depth n. p = 0.5 is the hardest (widest) case; p = 0.1 / 0.9 are easier.
#
# The payoff: overlay the depth distribution that surveillance ACTUALLY uses
# (ag_n_tested) to show where current practice sits on the curve.
# -----------------------------------------------------------------------------


depth_grid <- expand_grid(n = 1:12, p = c(0.1, 0.3, 0.5)) %>%
  rowwise() %>%
  mutate(
    exp_post_sd    = expected_sd_at_depth(n, p),
    ci95_halfwidth = 1.96 * exp_post_sd          # approx; Beta is skewed at small n
  ) %>%
  ungroup() %>%
  mutate(p_lab = factor(sprintf("true freq = %.1f", p)))

write_csv(depth_grid, "results/noise/precision_by_depth.csv")

# Empirical sampling depth currently in use.
#
# BUG FIX: this used ag_n_tested, which counts DRUG-TESTS, not colonies. AG has
# three agents, so one colony tested against all three scores n = 3 - and the
# distribution gave itself away as multiples of three (3 -> 611 visits,
# 6 -> 350, 9 -> 66, i.e. 1, 2 and 3 colonies). Every one of the 611 visits at
# "n = 3" contributed exactly ONE colony. The panel was therefore labelled
# "colonies" while plotting drug-tests, and the "current practice" reference
# line sat at 3 when the true median depth is 1. That made current sampling
# look three times deeper than it is, against a precision curve whose x axis
# really is colonies. Counting colonies directly fixes both.
ag_drugs_depth <- c("amikacin", "tobramycin", "gentamicin")
depth_emp <- pa_isolates %>%
  mutate(.vd = as.Date(date_collected_from_subject),
         # if_any(), not sapply(.data[[x]]) - the data-mask pronoun does not
         # survive inside sapply's closure ("object 'x' not found").
         .tested = as.integer(if_any(all_of(ag_drugs_depth),
                                     ~ .x %in% c("Susceptible", "Intermediate", "Resistant")))) %>%
  group_by(individual_id, .vd) %>%
  summarise(n_colonies = sum(.tested), .groups = "drop") %>%
  filter(n_colonies >= 1) %>%
  transmute(n = pmin(n_colonies, 12L))

med_depth <- median(depth_emp$n)
pct_shallow <- round(100 * mean(depth_emp$n <= 2), 1)

cat("\n6A precision-by-depth:\n")
cat("  Median AG sampling depth in pilot:", med_depth, "COLONIES/visit\n")
cat("  Visits sampled at n <= 2 colonies:", pct_shallow, "%\n")
cat("  (for contrast, median ag_n_tested = ",
    median(pa_visit_full$ag_n_tested[pa_visit_full$ag_n_tested >= 1], na.rm = TRUE),
    " DRUG-TESTS, which is what this panel used to plot)\n", sep = "")

# Top panel: precision curves (the load-bearing element)
fig6a_top <- ggplot(depth_grid, aes(x = n, y = ci95_halfwidth, color = p_lab)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_vline(xintercept = med_depth, linetype = "dashed", color = "grey40") +
  annotate("text", x = med_depth + 0.2, y = max(depth_grid$ci95_halfwidth) * 0.95,
           hjust = 0, size = 4, color = "grey30",
           label = sprintf("current median depth = %.0f", med_depth)) +
  scale_x_continuous(breaks = 1:12) +
  scale_color_manual(values = c("#1A9641", "#FC8D59", "#D7301F"), name = NULL) +
  labs(x = NULL,
       y = "Approx. 95% CI half-width\non resistance frequency",
       title = "How many colonies does a culture need?",
       subtitle = "Precision improves steeply then plateaus; dashed line = current practice") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        legend.position = c(0.78, 0.8),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13))

# Bottom panel: where current sampling actually sits
fig6a_bot <- ggplot(depth_emp, aes(x = n)) +
  geom_bar(fill = "grey55", width = 0.8) +
  geom_vline(xintercept = med_depth, linetype = "dashed", color = "grey40") +
  scale_x_continuous(breaks = 1:12, limits = c(0.4, 12.6)) +
  labs(x = "Colonies tested per culture (sampling depth n)",
       y = "Visits",
       subtitle = sprintf("Current practice: %s%% of visits sampled at n <= 2", pct_shallow)) +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 13))

fig6a <- fig6a_top / fig6a_bot + plot_layout(heights = c(2, 1))

ggsave("figures/noise/fig6a_precision_by_depth.pdf", fig6a,
       width = 9, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# 6B. NOISE-PROPAGATION DOWN-SAMPLING EXPERIMENT  (NEW — run interactively)
#
# Idea: treat each visit's FULL set of tested isolates as ground truth, then
# artificially thin to m = 1, 2, 3, 5 colonies and re-estimate the same
# quantity with each of the three method families. The three should degrade at
# DIFFERENT rates as m shrinks — that is the "which model do you trust under
# noise" figure.
#
#   GLMM lens   : persistence coefficient (effect of prior resistance), with the
#                 PREDICTOR measured from only m colonies at t-1 (mirrors fig3,
#                 now as a function of depth). Expect attenuation toward 0.
#   Markov lens : steady-state resistant fraction implied by S<->R transitions
#                 built from m-colony state calls. Expect pull toward 0.5.
#   CCM lens    : cross-map skill (rho) for a chosen coupling on m-colony series.
#                 Expect the steepest decay (CCM is the most noise-fragile).
#
# This refits brms models and is newly written — it has NOT been executed here.
# Run it interactively, check each piece, then assemble the combined figure.
# -----------------------------------------------------------------------------

RUN_SECTION_6B <- FALSE   # flip to TRUE to run the down-sampling experiment

if (RUN_SECTION_6B) {
  
  set.seed(42)
  ag_drugs <- c("amikacin", "gentamicin", "tobramycin")
  
  # ---- Per-isolate AG resistance call (same rule as fig1 variance decomp) ----
  iso_ag <- pa_isolates %>%
    rowwise() %>%
    mutate(
      ag_test_count = sum(c_across(all_of(ag_drugs)) %in%
                            c("Resistant", "Intermediate", "Susceptible")),
      ag_R = if_else(
        ag_test_count == 0, NA_integer_,
        as.integer(any(c_across(all_of(ag_drugs)) %in% c("Resistant", "Intermediate")))
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(ag_R)) %>%
    transmute(individual_id, date_collected_from_subject, ag_R)
  
  # List of per-visit isolate calls
  visit_calls <- iso_ag %>%
    group_by(individual_id, date_collected_from_subject) %>%
    summarise(calls = list(ag_R), n_iso = n(), .groups = "drop")
  
  # Subsample each visit to depth m; return visit-level summary
  subsample_to_depth <- function(m) {
    visit_calls %>%
      mutate(
        draw = map(calls, ~ if (length(.x) <= m) .x else sample(.x, m)),
        n_tested_sub    = map_int(draw, length),
        n_resistant_sub = map_int(draw, sum),
        any_resistant_sub = as.integer(n_resistant_sub > 0)
      ) %>%
      arrange(individual_id, date_collected_from_subject) %>%
      group_by(individual_id) %>%
      mutate(prev_any_resistant_sub = lag(any_resistant_sub)) %>%
      ungroup() %>%
      select(individual_id, date_collected_from_subject,
             n_tested_sub, n_resistant_sub, any_resistant_sub,
             prev_any_resistant_sub)
  }
  
  depths <- c(1, 2, 3, 5, 99)   # 99 == "use all isolates" (full sample)
  

# Slide-sized top panel on its own (Talk_figures_v2 shipped only this one).
ggsave("figures/noise/fig6a_precision_by_depth.pdf", fig6a_top,
       width = 5, height = 4, dpi = 300, bg = "transparent")


# -----------------------------------------------------------------------------
# 3.4  Down-sampling noise propagation
# 
# -----------------------------------------------------------------------------

  # ---- Markov lens: cheap, repeat subsampling to get CIs -------------------
  markov_steady_at_depth <- function(m, reps = 100) {
    vals <- map_dbl(seq_len(reps), function(i) {
      d <- subsample_to_depth(m) %>% filter(!is.na(prev_any_resistant_sub))
      p_sr <- mean(d$any_resistant_sub[d$prev_any_resistant_sub == 0])
      p_rs <- mean(d$any_resistant_sub[d$prev_any_resistant_sub == 1] == 0)
      if (is.na(p_sr) || is.na(p_rs) || (p_sr + p_rs) == 0) return(NA_real_)
      p_sr / (p_sr + p_rs)
    })
    tibble(m = m,
           steady_R = mean(vals, na.rm = TRUE),
           lo = quantile(vals, 0.025, na.rm = TRUE),
           hi = quantile(vals, 0.975, na.rm = TRUE))
  }
  markov_curve <- map_dfr(depths, markov_steady_at_depth)
  write_csv(markov_curve, "results/noise/downsample_markov.csv")
  
  # ---- GLMM lens: one beta-binomial fit per depth (cached) -----------------
  # Outcome = full-sample current resistance; predictor = m-colony prior state.
  covars <- pa_visit_full %>%
    select(individual_id, date_collected_from_subject,
           ag_n_resistant, ag_n_tested,
           prev_has_any_mucoid, prev_has_mrsa, prev_has_mssa,
           prev_aminoglycoside_prescribed)
  
  glmm_coef_at_depth <- function(m) {
    d <- subsample_to_depth(m) %>%
      inner_join(covars, by = c("individual_id", "date_collected_from_subject")) %>%
      filter(ag_n_tested >= 2,
             !is.na(prev_any_resistant_sub),
             !is.na(prev_has_any_mucoid), !is.na(prev_has_mrsa),
             !is.na(prev_has_mssa), !is.na(prev_aminoglycoside_prescribed))
    
    fit <- brm(
      ag_n_resistant | trials(ag_n_tested) ~
        prev_any_resistant_sub + prev_has_any_mucoid +
        prev_has_mrsa + prev_has_mssa + prev_aminoglycoside_prescribed +
        (1 | individual_id),
      data = d, family = beta_binomial(),
      prior = c(prior(normal(0, 1.5), class = "b"),
                prior(normal(0, 2),   class = "Intercept"),
                prior(exponential(1), class = "sd"),
                prior(gamma(2, 0.1),  class = "phi")),
      chains = 4, cores = 4, iter = 4000, warmup = 2000,
      control = list(adapt_delta = 0.95), seed = 42,
      file = paste0("models/noise/glmm_downsample_m", m),
      file_refit = "on_change"
    )
    dr <- as_draws_df(fit)$b_prev_any_resistant_sub
    tibble(m = m, estimate = median(dr),
           lo = quantile(dr, 0.025), hi = quantile(dr, 0.975))
  }
  glmm_curve <- map_dfr(depths, glmm_coef_at_depth)
  write_csv(glmm_curve, "results/noise/downsample_glmm.csv")
  
  # ---- CCM lens: optional, slowest -----------------------------------------
  # Left as a stub: build a per-patient subsampled resistance series, run CCM
  # for one coupling (e.g. AG prescription -> AG resistance) at each depth, take
  # the terminal rho. Use run_ccm_pair() from AMR_common.R (the old helper went with the placeholder CCM).
  
  # ---- Combined "degradation" figure ---------------------------------------
  glmm_full   <- glmm_curve$estimate[glmm_curve$m == 99]
  markov_full <- markov_curve$steady_R[markov_curve$m == 99]
  
  combined <- bind_rows(
    glmm_curve   %>% transmute(m, method = "GLMM persistence coef",
                               signal = estimate / glmm_full),
    markov_curve %>% transmute(m, method = "Markov steady-state (dist. from 0.5)",
                               signal = abs(steady_R - 0.5) / abs(markov_full - 0.5))
  ) %>% filter(m < 99)
  
  fig6b <- ggplot(combined, aes(x = m, y = signal, color = method)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey60") +
    geom_line(linewidth = 1.2) + geom_point(size = 2.5) +
    scale_x_continuous(breaks = c(1, 2, 3, 5)) +
    labs(x = "Colonies tested per culture (sampling depth)",
         y = "Signal recovered (fraction of full-sample value)",
         title = "Different models break at different sampling depths",
         subtitle = "Each method's estimate as a fraction of its full-sample value",
         color = NULL) +
    theme(plot.title = element_text(size = 16, face = "bold"),
          legend.position = "bottom")
  
  ggsave("figures/noise/fig6b_downsample_propagation.pdf", fig6b,
         width = 9, height = 6, dpi = 300)
}




# #############################################################################
# SECTION 4: HETEROGENEITY  -  modelling the spread, not just the mean
# mu vs phi, naive vs measurement-error-aware, and staph-stratified dispersion.
# Talk_figures_v2 labelled this block 'POSTER figures only / superseded'. That
# comment is wrong - its own slide map (v2 lines 14-15) sends fig2 to Slide 6
# and fig3 to Slide 7, and priors_dist_betabin feeds Section 5. Kept in full.
# #############################################################################

#
# This is the conceptual pivot of the talk. Three axes of inference:
#   AXIS 1  what is TRUE in the lung      (biological clonal heterogeneity)
#   AXIS 2  how well we can MEASURE it    (sampling precision; fixable)
#   AXIS 3  what it REPRESENTS clinically (AST/in-vivo validity gap; not fixable
#                                          by more sampling)
# Section 5 builds the figures; Section 6 turns Axis 2 into a design ask.
# =============================================================================
# =============================================================================
# AG subset: ag_n_tested >= 2 with all 5 predictors non-missing.
# CAUTION - ag_n_tested counts DRUG-TESTS, not colonies. AG has three agents,
# so a single colony tested on amikacin + tobramycin already satisfies >= 2:
# 54% of the visits this filter keeps contributed only ONE colony. It is
# therefore NOT a filter for "visits where heterogeneity is definable", which
# is what this comment used to claim. Section 10.2 refits the same model with
# a true colony denominator (ag_k | trials(ag_n_iso)) and n_isolates >= 2.
het_data <- pa_visit_full %>%
  filter(
    ag_n_tested >= 2,
    !is.na(prev_ag_any_resistant),
    !is.na(prev_has_any_mucoid),
    !is.na(prev_has_mrsa),
    !is.na(prev_has_mssa),
    !is.na(prev_aminoglycoside_prescribed)
  )

cat("v5 AG heterogeneity dataset:\n")
cat("  N visits:", nrow(het_data), "\n")
cat("  Patients:", n_distinct(het_data$individual_id), "\n")
cat("  Total AG isolates tested:", sum(het_data$ag_n_tested), "\n")
cat("  Predictor distributions:\n")
cat("    prev_ag_any_resistant         :", round(100*mean(het_data$prev_ag_any_resistant),1), "%\n")
cat("    prev_has_any_mucoid           :", round(100*mean(het_data$prev_has_any_mucoid),1), "%\n")
cat("    prev_has_mrsa                 :", round(100*mean(het_data$prev_has_mrsa),1), "%\n")
cat("    prev_has_mssa                 :", round(100*mean(het_data$prev_has_mssa),1), "%\n")
cat("    prev_aminoglycoside_prescribed:", round(100*mean(het_data$prev_aminoglycoside_prescribed),1), "%\n\n")


# Halfeye plot

# =============================================================================
# STEP 1: THREE-LEVEL VARIANCE DECOMPOSITION (Figure 1, unchanged)
# =============================================================================

ag_drugs <- c("amikacin", "gentamicin", "tobramycin")

iso_ag <- pa_isolates %>%
  rowwise() %>%
  mutate(
    ag_test_count = sum(c_across(all_of(ag_drugs)) %in%
                          c("Resistant", "Intermediate", "Susceptible")),
    ag_R          = if_else(
      ag_test_count == 0,
      NA_integer_,
      as.integer(any(c_across(all_of(ag_drugs)) %in%
                       c("Resistant", "Intermediate")))
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(ag_R)) %>%
  mutate(visit_key = paste(individual_id, date_collected_from_subject, sep = "_"))

m_var <- brm(
  ag_R ~ 1 + (1 | individual_id) + (1 | visit_key),
  data = iso_ag,
  family = bernoulli(),
  prior = c(
    prior(normal(0, 2), class = "Intercept"),
    prior(exponential(1), class = "sd")
  ),
  chains = 4, cores = 4,
  iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95),
  seed = 42,
  file = "models/heterogeneity/m_var_decomp",
  file_refit = "on_change"
)

vc <- VarCorr(m_var)
sigma_pat   <- vc$individual_id$sd[, "Estimate"]^2
sigma_visit <- vc$visit_key$sd[, "Estimate"]^2
sigma_resid <- pi^2 / 3
total       <- sigma_pat + sigma_visit + sigma_resid

var_decomp <- tibble(
  level = factor(
    c("Between patients",
      "Between isolates (within visit)",
      "Between visits (within patient)"),
    levels = c("Between patients",
               "Between isolates (within visit)",
               "Between visits (within patient)")
  ),
  proportion = c(sigma_pat, sigma_resid, sigma_visit) / total
)
write_csv(var_decomp, "results/heterogeneity/variance_decomposition.csv")

fig1 <- ggplot(var_decomp,
               aes(x = level, y = proportion * 100, fill = level)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", proportion * 100)),
            vjust = -0.5, size = 5.5, fontface = "bold") +
  scale_fill_manual(values = c("#1B7837", "#5AAE61", "#A6DBA0"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Share of total variance (%)",
       title = "Variance partitioning of aminoglycoside resistance",
       subtitle = "Within-visit slice (22.6%) is what single-clone sampling discards") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14), 
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 13))

ggsave("figures/heterogeneity/fig1_variance_decomposition.pdf",
       fig1, width = 10, height = 6, dpi = 300)

# =============================================================================
# STEP 2: DISTRIBUTIONAL BETA-BINOMIAL - AG with FOCUSED predictor set
# =============================================================================

formula_main <- bf(
  ag_n_resistant | trials(ag_n_tested) ~
    prev_ag_any_resistant +
    prev_has_any_mucoid +
    prev_has_mrsa +
    prev_has_mssa +
    prev_aminoglycoside_prescribed +
    (1 | individual_id),
  phi ~
    prev_ag_any_resistant +
    prev_has_any_mucoid +
    prev_has_mrsa +
    prev_has_mssa +
    prev_aminoglycoside_prescribed +
    (1 | individual_id),
  family = beta_binomial()
)

m_main_v5 <- brm(
  formula_main,
  data = het_data,
  prior = priors_dist_betabin,
  chains = 4, cores = 4,
  iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 42,
  file = "models/heterogeneity/m_main_v5",
  file_refit = "on_change"
)

print(summary(m_main_v5))

main_draws <- extract_draws_long(m_main_v5, predictor_recodes)
main_summary <- extract_disp_effects(m_main_v5, predictor_recodes)
write_csv(main_summary, "results/heterogeneity/main_effects_v5.csv")

# Lock factor order using mu medians (theme-grouped: phenotype, species, drug)
# Manual ordering reads more naturally for the audience than data-driven sort:
main_var_order <- c(
  "Prior AG prescription",
  "Prior MSSA",
  "Prior MRSA",
  "Prior mucoid phenotype",
  "Prior AG resistance"
)

fig2_mu  <- plot_posteriors(
  main_draws %>% filter(submodel == "mu"),
  "Where the distribution sits: logit(mu)",
  main_var_order
) + theme(plot.title = element_text(size = 16), axis.text.x = element_text(size = 16), 
          axis.text.y = element_text(size = 16))

fig2_phi <- plot_posteriors(
  main_draws %>% filter(submodel == "phi"),
  "How spread out it is: log(phi)   (+ = homogeneous)",
  main_var_order
) + theme(plot.title = element_text(size = 16), axis.text.y = element_text("none"), 
          axis.text.x = element_text(size = 16)
)

fig2 <- (fig2_mu | fig2_phi) +
  plot_annotation(
    title = "Predictor effects on resistance: location vs spread",
    subtitle = sprintf(
      "AG resistance, N = %d visits, %d patients | densities = full posteriors | bars = 66%% / 95%% CIs",
      nrow(het_data), n_distinct(het_data$individual_id)
    ),
    theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave("figures/heterogeneity/fig2_mu_vs_phi_v5.pdf",
       fig2, width = 14, height = 6, dpi = 300)

# Posterior predictive check
pp_main <- pp_check(m_main_v5, ndraws = 100, type = "hist") +
  ggtitle("Posterior predictive check: distributional beta-binomial")
ggsave("figures/heterogeneity/fig2b_ppc_v5.pdf",
       pp_main, width = 10, height = 6, dpi = 300)

# =============================================================================
# STEP 3: NAIVE vs MEASUREMENT-ERROR-AWARE (Figure 3, redesigned)
#
# Two contrasts shown together:
#   (A) Top panel: same FIVE predictors, but the prior-resistance row is
#       different between models - naive uses binary, ME uses latent theta.
#       The OTHER FOUR predictors are identical between models, so they
#       act as a control/sanity check (their coefficients should NOT shift).
#   (B) Bottom panel: scatter showing what the predictor MEANS in each model.
#       Naive: the binary value (0 or 1). ME: the posterior mean of latent
#       theta, which is shrunken toward 0.5 for low-n visits.
#       This makes the "sampling precision matters" point visceral.
# =============================================================================

# Build the latent-theta predictor for the ME model
me_data <- pa_visit_full %>%
  arrange(individual_id, date_collected_from_subject) %>%
  group_by(individual_id) %>%
  mutate(
    prev_n_resistant_lag = lag(ag_n_resistant),
    prev_n_tested_lag    = lag(ag_n_tested)
  ) %>%
  ungroup() %>%
  filter(
    ag_n_tested >= 2,
    !is.na(prev_n_resistant_lag),
    !is.na(prev_n_tested_lag),
    prev_n_tested_lag >= 1,
    !is.na(prev_has_any_mucoid),
    !is.na(prev_has_mrsa),
    !is.na(prev_has_mssa),
    !is.na(prev_aminoglycoside_prescribed)
  ) %>%
  mutate(
    # DESIGN-EFFECT CORRECTION.
    # ag_n_tested counts DRUG-tests, not independent tests: one colony assayed
    # against amikacin / tobramycin / gentamicin contributes 3. Those three
    # calls are correlated within a colony (rho measured in Section 14), so
    # treating them as 3 independent draws overstates the precision of theta
    # and understates the shrinkage - which is the whole point of the ME model.
    #
    #   n_eff = n / (1 + (k - 1) * rho)      k = 3 AG drugs
    #
    # k_eff rescales the numerator so the observed proportion is preserved and
    # only the precision changes. rho comes from ag_within_colony_rho(), the
    # same helper Section 14 uses, so the two cannot drift apart.
    .rho_w  = ag_within_colony_rho(pa_isolates),
    .k_drug = 3,
    .deff   = 1 + (.k_drug - 1) * .rho_w,
    n_eff   = prev_n_tested_lag / .deff,
    k_eff   = prev_n_resistant_lag * n_eff / prev_n_tested_lag,
    # Posterior mean of theta (latent prev resistance frequency, Beta(1,1) prior)
    prev_theta_pm = (1 + k_eff) / (2 + n_eff),
    # Posterior SD of theta
    prev_theta_psd = sqrt(
      (1 + k_eff) * (1 + n_eff - k_eff) /
        ((2 + n_eff)^2 * (3 + n_eff))
    ),
    # keep the uncorrected version so the change is auditable, not silent
    prev_theta_pm_drugtests  = (1 + prev_n_resistant_lag) / (2 + prev_n_tested_lag),
    prev_theta_psd_drugtests = sqrt(
      (1 + prev_n_resistant_lag) *
        (1 + prev_n_tested_lag - prev_n_resistant_lag) /
        ((2 + prev_n_tested_lag)^2 * (3 + prev_n_tested_lag))
    )
  )

cat(sprintf("\n  Design-effect correction: rho = %.2f, deff = %.2f\n",
            me_data$.rho_w[1], me_data$.deff[1]))
cat(sprintf("  median prev_theta_psd: %.3f (drug-tests) -> %.3f (effective tests)\n",
            median(me_data$prev_theta_psd_drugtests), median(me_data$prev_theta_psd)))
cat("  NOTE (measured 2026-08-27): the correction does NOT strengthen the\n")
cat("  contrast here. prior-resistance went 3.99 [1.68, 5.85] -> 2.40 [-0.38,\n")
cat("  5.09]: the point estimate fell and the interval now spans zero. A less\n")
cat("  precisely known predictor is also a less identifiable one, and on 232\n")
cat("  visits that dominates the attenuation correction. Report the corrected\n")
cat("  number; do not claim the correction helped the case.\n")

cat("\n=== Step 3: ME comparison data ===\n")
cat("N visits:", nrow(me_data),
    "| Patients:", n_distinct(me_data$individual_id), "\n")
cat("prev_n_tested distribution (sampling depth at t-1):\n")
print(summary(me_data$prev_n_tested_lag))

# Naive model: binary prior R
m_naive_v5 <- brm(
  bf(
    ag_n_resistant | trials(ag_n_tested) ~
      prev_ag_any_resistant +
      prev_has_any_mucoid +
      prev_has_mrsa +
      prev_has_mssa +
      prev_aminoglycoside_prescribed +
      (1 | individual_id),
    family = beta_binomial()
  ),
  data = me_data,
  prior = c(
    prior(normal(0, 1.5), class = "b"),
    prior(normal(0, 2),   class = "Intercept"),
    prior(exponential(1), class = "sd"),
    prior(gamma(2, 0.1),  class = "phi")
  ),
  chains = 4, cores = 4,
  iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95),
  seed = 42,
  file = "models/heterogeneity/m_naive_v5",
  file_refit = "on_change"
)

# ME-aware model: latent theta replacing binary prior R
m_me_v5 <- brm(
  bf(
    ag_n_resistant | trials(ag_n_tested) ~
      me(prev_theta_pm, prev_theta_psd) +
      prev_has_any_mucoid +
      prev_has_mrsa +
      prev_has_mssa +
      prev_aminoglycoside_prescribed +
      (1 | individual_id),
    family = beta_binomial()
  ),
  data = me_data,
  prior = c(
    prior(normal(0, 1.5), class = "b"),
    prior(normal(0, 2),   class = "Intercept"),
    prior(exponential(1), class = "sd"),
    prior(gamma(2, 0.1),  class = "phi"),
    prior(normal(0.4, 0.3), class = "meanme"),
    prior(exponential(2),   class = "sdme")
  ),
  chains = 4, cores = 4,
  iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 13),
  seed = 42,
  file = "models/heterogeneity/m_me_v5",
  file_refit = "on_change"
)

print(summary(m_naive_v5))
print(summary(m_me_v5))

# Tidy: pull effects with consistent labeling for both models
naive_eff <- as_draws_df(m_naive_v5) %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  filter(param != "b_Intercept") %>%
  mutate(
    model = "Naive\n(prior R as binary)",
    variable = case_when(
      param == "b_prev_ag_any_resistant" ~ "Prior resistance",
      TRUE ~ str_remove(param, "^b_")
    )
  )

me_eff <- as_draws_df(m_me_v5) %>%
  select(starts_with("b_"), starts_with("bsp_")) %>%
  pivot_longer(everything(), names_to = "param", values_to = "value") %>%
  filter(!str_detect(param, "Intercept")) %>%
  mutate(
    model = "ME-aware\n(prior R as latent theta)",
    variable = case_when(
      str_detect(param, "meprev_theta_pm") ~ "Prior resistance",
      TRUE ~ str_remove(param, "^b_")
    )
  )

compare_me <- bind_rows(naive_eff, me_eff) %>%
  mutate(
    pretty_var = recode(variable,
                        "Prior resistance"               = "Prior resistance",
                        "prev_has_any_mucoid"            = "Prior mucoid phenotype",
                        "prev_has_mrsa"                  = "Prior MRSA",
                        "prev_has_mssa"                  = "Prior MSSA",
                        "prev_aminoglycoside_prescribed" = "Prior AG prescription"
    )
  )

# Variable order: prior resistance at top, then theme groupings
me_var_order <- c(
  "Prior AG prescription",
  "Prior MSSA",
  "Prior MRSA",
  "Prior mucoid phenotype",
  "Prior resistance"
)

# === Panel A: same five rows, two posteriors per row, side-by-side
# Critical visual: only the "Prior resistance" row should differ between
# models. The other four rows are identical predictors so their posteriors
# should be visually similar (this is the control/sanity check the audience
# needs to see).

fig3a <- compare_me %>%
  mutate(
    pretty_var = factor(pretty_var, levels = me_var_order),
    model = factor(model, levels = c("Naive\n(prior R as binary)",
                                     "ME-aware\n(prior R as latent theta)"))
  ) %>%
  group_by(pretty_var, model) %>%
  mutate(
    median_est  = median(value),
    direction   = if_else(median_est > 0, "+", "-"),
    pd          = max(mean(value > 0), mean(value < 0)),
    significant = pd > 0.975,
    is_focal    = pretty_var == "Prior resistance"
  ) %>%
  ungroup() %>%
  ggplot(aes(x = value, y = pretty_var)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  stat_halfeye(
    aes(fill = direction, alpha = is_focal),
    .width = c(0.66, 0.95),
    slab_color = NA, point_size = 2,
    interval_size_range = c(0.6, 1.6)
  ) +
  scale_fill_manual(values = c("+" = "#B2182B", "-" = "#2166AC"), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.35), guide = "none") +
  facet_wrap(~ model, scales = "free_x", nrow = 1) +
  # Greek spelled out, NOT "\u03bc" and NOT expression(mu). The Google Fonts
  # build of Lato registered by showtext has no Greek block, so mu/theta/beta
  # and the arrow render as tofu boxes. showtext routes plotmath through that
  # same font, so expression() does not rescue it either - verified by test
  # render. Only ASCII is safe here; x, +/- and <= do survive.
  labs(x = "Effect on logit(mu)", y = NULL,
       title = "Same predictors, two ways of encoding prior resistance",
       subtitle = "Highlighted rows differ between models; faded rows are identical predictors (sanity check)") +
  theme(strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 14, face = "bold"), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 16))

# === Panel B: scatter showing what the predictor LOOKS LIKE under each model
# This is the visceral "sampling precision matters" panel.
# X-axis: observed prior proportion (k/n) - what naive treats as the predictor
# Y-axis: posterior mean of latent theta - what ME treats as the predictor
# Color: prev_n_tested (sampling depth)
# Diagonal: y = x (where they would coincide if there were no shrinkage)

# The old version binned prev_n_tested_lag as 1-2 / 3-5 / 6+ and called it
# "sampling depth". Two problems: it is DRUG-TESTS, not colonies (AG has three
# agents, so one colony scores 3), and against this data the "1-2" bin is EMPTY
# and 217 of 232 points fall in "6+" - one colour, no visible divergence.
#
# Rebuilt two ways. Points are coloured by COLONIES (the honest unit, and the
# cohort only reaches 4). The fan of grey reference lines shows where the
# ME predictor would sit at depths the cohort does not have: shrinkage is
# deterministic given n, pm = (1 + k)/(2 + n), which in terms of the observed
# proportion is a straight line of slope n/(n+2) and intercept 1/(n+2). That is
# what makes the divergence legible - you can see a 1-colony visit pulled almost
# to 0.5 while a 20-test visit sits nearly on the diagonal, even though this
# cohort never samples that deep.
# Reference lines labelled in COLONIES, not drug-tests. Two reasons the old
# "60 tests (20 colonies)" labelling was wrong:
#   - "tests" and "colonies tested" are different quantities and putting both in
#     one label invites reading them as interchangeable;
#   - there are only THREE aminoglycosides, so "60 tests" reads as 60 drugs. It
#     meant 20 colonies x 3 drugs, which is a realistic deep-sampling target
#     (Vanderwoude et al. sampled ~75 colonies) but not 60 distinct agents.
#
# The model's n is drug-tests: m_me_v5 uses prev_n_tested_lag, so c colonies
# enter as n = 3c. That assumes the three AG drugs on one colony are three
# independent draws, and they are not - within a colony they agree 63% of the
# time and correlate 0.45 (ami-tob) to 0.70 (ami-gen). The model therefore
# credits each colony with about three times the information it carries, and
# shrinks too little. The dashed line shows where one colony would sit if the
# COLONY were the independent unit (n = 1), which is the denominator Section
# 10.2 refits on. The gap between it and the "1 colony" solid line is the
# overstatement.
colony_levels <- c(1, 2, 4, 10, 20)
shrink_lines <- tibble::tibble(n_colonies = colony_levels) %>%
  mutate(n_model   = 3 * n_colonies,
         slope     = n_model / (n_model + 2),
         intercept = 1 / (n_model + 2),
         lab = sprintf("%d %s", n_colonies,
                       ifelse(n_colonies == 1, "colony", "colonies")))
# same thing, but counting the colony (not the drug-test) as the unit
shrink_conservative <- tibble::tibble(n_colonies = 1) %>%
  mutate(slope = n_colonies / (n_colonies + 2), intercept = 1 / (n_colonies + 2))

fig3b <- me_data %>%
  mutate(
    obs_prop = prev_n_resistant_lag / prev_n_tested_lag,
    n_colonies = prev_n_tested_lag / 3,
    n_bin = case_when(
      n_colonies <= 1 ~ "1 colony",
      n_colonies <= 2 ~ "2 colonies",
      TRUE            ~ "3-4 colonies"
    ) %>% factor(levels = c("1 colony", "2 colonies", "3-4 colonies"))
  ) %>%
  ggplot(aes(x = obs_prop, y = prev_theta_pm)) +
  geom_abline(data = shrink_lines, aes(slope = slope, intercept = intercept),
              colour = "grey65", linewidth = 0.5) +
  geom_abline(data = shrink_conservative, aes(slope = slope, intercept = intercept),
              colour = pal[["red"]], linewidth = 0.6, linetype = "dotted") +
  geom_text(data = shrink_conservative, aes(x = 1.02, y = slope + intercept),
            label = "1 colony, counted as 1 draw", hjust = 0, size = 3.0,
            colour = pal[["red"]]) +
  geom_text(data = shrink_lines, aes(x = 1.02, y = slope + intercept, label = lab),
            hjust = 0, size = 3.1, colour = "grey40") +
  geom_abline(slope = 1, intercept = 0, color = "grey20", linetype = "dashed") +
  annotate("text", x = 1.02, y = 1.06, hjust = 0, size = 3.1, fontface = "bold",
           label = "no shrinkage (perfect measurement)") +
  geom_jitter(aes(color = n_bin), width = 0.02, height = 0,
              size = 2.4, alpha = 0.7) +
  scale_color_manual(
    values = c("1 colony"     = "#D7301F",
               "2 colonies"   = "#FC8D59",
               "3-4 colonies" = "#1A9641"),
    name = "Prior visit\ncolonies tested"
  ) +
  scale_x_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, 0.25)) +
  coord_cartesian(clip = "off") +
  labs(x = "Observed prior proportion (k/n)",
       y = "Posterior mean of latent theta",
       title = "What the predictor means in each model",
       subtitle = paste("Grey: colonies as the model counts them (3 drug-tests each). Dotted red: one colony as one draw.",
                        "\nThis cohort never exceeds 4 colonies.")) +
  theme(plot.title = element_text(size = 14, face = "bold"),
        legend.position = c(0.82, 0.22),
        plot.margin = margin(5.5, 90, 5.5, 5.5),
        axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 16), 
        axis.title.y = element_text(size = 16))

fig3 <- fig3a / fig3b +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Naive vs measurement-error-aware: respecting sampling precision",
    subtitle = sprintf("N = %d visits, %d patients", 
                       nrow(me_data), n_distinct(me_data$individual_id)),
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave("figures/heterogeneity/fig3_naive_vs_me_v5.pdf",
       fig3, width = 13, height = 11, dpi = 300)

# =============================================================================
# STEP 4: STAPH-STRATIFIED MODEL (Figure 4 - ecological context)
#
# Same model, fit separately for Staph+ and Staph- visits.
# Question: does the structure of dispersion differ between PA-only and
# PA+Staph communities? Posterior shape comparison across strata = effect
# modification by ecological context.
# =============================================================================

# For stratification, build a single Staph indicator from MRSA OR MSSA.
# We KEEP MRSA and MSSA separate as covariates within strata? No - within
# each stratum, "has Staph" is constant by definition, so MRSA/MSSA only
# vary in the Staph+ stratum. Cleaner to stratify on EITHER staph and
# include only the cross-stratum-varying predictors in each model.

het_data_strat <- het_data %>%
  mutate(has_any_staph = pmax(prev_has_mrsa, prev_has_mssa, na.rm = TRUE))

cat("\nStaph stratification:\n")
cat("  Staph+ visits:", sum(het_data_strat$has_any_staph == 1), "\n")
cat("  Staph- visits:", sum(het_data_strat$has_any_staph == 0), "\n")


m_staph_pos_v5 <- fit_staph_stratum(het_data_strat, 1, "pos")
m_staph_neg_v5 <- fit_staph_stratum(het_data_strat, 0, "neg")

stratum_recodes <- c(
  prev_ag_any_resistant         = "Prior AG resistance",
  prev_has_any_mucoid           = "Prior mucoid phenotype",
  prev_aminoglycoside_prescribed = "Prior AG prescription"
)

n_pos <- sum(het_data_strat$has_any_staph == 1)
n_neg <- sum(het_data_strat$has_any_staph == 0)

# fit_staph_stratum() returns NULL for a stratum with fewer than 80 visits, and
# nothing downstream used to check - as_draws_df(NULL) then killed the run with
# "Don't know how to transform an object of class 'NULL'". This bites now that
# het_data is complete-case on prev_has_any_mucoid and
# prev_aminoglycoside_prescribed, which are honestly NA since the prep fixes
# (662 and 555 visits respectively): het_data is 232 visits, so Staph+ is 73 and
# falls under the guard. Skip the figure and say so, rather than die here.
if (is.null(m_staph_pos_v5) || is.null(m_staph_neg_v5)) {
  cat("\nSKIPPING fig4_staph_stratified_v5: stratum too small to fit",
      sprintf("(Staph+ n=%d, Staph- n=%d; fit_staph_stratum needs >= 80).\n", n_pos, n_neg))
  cat("  het_data is", nrow(het_data), "visits after complete-case filtering.\n")
  cat("  To restore this figure, relax the het_data filter - dropping\n")
  cat("  prev_has_any_mucoid alone recovers most of the rows - or accept that\n")
  cat("  the stratification is not estimable on the honestly-missing data.\n")
} else {

draws_pos <- extract_draws_long(m_staph_pos_v5, stratum_recodes) %>%
  mutate(stratum = sprintf("Staph+ (n=%d)", n_pos))
draws_neg <- extract_draws_long(m_staph_neg_v5, stratum_recodes) %>%
  mutate(stratum = sprintf("Staph- (n=%d)", n_neg))

stratum_draws <- bind_rows(draws_pos, draws_neg) %>%
  mutate(stratum = factor(stratum, levels = c(sprintf("Staph- (n=%d)", n_neg),
                                              sprintf("Staph+ (n=%d)", n_pos))))

stratum_var_order <- c(
  "Prior AG prescription",
  "Prior mucoid phenotype",
  "Prior AG resistance"
)


fig4_mu  <- plot_stratum(stratum_draws, "mu",
                         "Where the distribution sits: logit(mu)")
fig4_phi <- plot_stratum(stratum_draws, "phi",
                         "How spread out it is: log(phi)  (+ = more homogeneous)")

fig4 <- (fig4_mu / fig4_phi) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title    = "Ecological context: Staph co-colonization changes the structure of resistance dynamics",
    subtitle = "Same model, fit separately for Staph+ and Staph- visits. Compare posterior shapes between strata.",
    theme    = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave("figures/heterogeneity/fig4_staph_stratified_v5.pdf",
       fig4, width = 13, height = 8, dpi = 300)

}  # end if (both staph strata fitted)

# =============================================================================
# DONE
# =============================================================================

cat("\n\n=== HETEROGENEITY ANALYSIS v5 COMPLETE ===\n")
cat("Four figures, consistent predictor set, ready for collaborator review:\n")
cat("  fig1_variance_decomposition.pdf - the motivation\n")
cat("  fig2_mu_vs_phi_v5.pdf           - location vs spread (5 predictors, 3 themes)\n")
cat("  fig3_naive_vs_me_v5.pdf         - sampling-precision correction (2 panels)\n")
cat("  fig4_staph_stratified_v5.pdf    - ecological context modulates dynamics\n")
cat("\nPredictor structure (consistent across all model figures):\n")
cat("  Phenotype:  prior AG resistance, prior mucoid\n")
cat("  Species:    prior MRSA, prior MSSA\n")
cat("  Drug:       prior AG prescription\n")

# =============================================================================

# #############################################################################
# SECTION 5: PER-CLASS MODELS
# Noise-aware composite, per-class standard + deep dive, clinical-status-aware.
# New in Talk_figures_v2; no counterpart in Synthesis.R.
# #############################################################################

# PURPOSE: apples-to-apples naive vs aware on the SAME composite outcome and
# predictor set, so the survivors that feed CCM come from the honest model.
# (fig2/fig3/fig4 stay AG-specific and unchanged; this is the composite analog
#  used only to choose what CCM should test.)
# Depends on: route_df, priors_beta_random, pal, lab_n  (all defined upstream).
# =============================================================================

# Latent theta for prior COMPOSITE resistance, from the already-lagged
# conservative counts (Beta(1,1) conjugate posterior; same formula as fig3).
# Computed outside the mutate on purpose: across()/all_of() inside a nested
# transmute() gets captured by the outer mutate's data mask and fails.
.panel9 <- c("amikacin", "gentamicin", "tobramycin", "meropenem", "cefepime",
             "ceftazidime", "levofloxacin", "piperacillin_tazobactam", "aztreonam")
.rho_c  <- ag_within_colony_rho(pa_isolates, drugs = .panel9)
.k_c    <- mean(rowSums(pa_isolates[.panel9] %>%
                          lapply(function(x) x %in% c("Susceptible", "Intermediate",
                                                      "Resistant")) %>% as.data.frame()))
.deff_c <- 1 + (.k_c - 1) * .rho_c

me_comp_df <- route_df %>%
  filter(!is.na(prev_n_resistant), !is.na(prev_n_tested), prev_n_tested >= 1) %>%
  mutate(
    # DESIGN-EFFECT CORRECTION, composite version.
    # Same conflation as m_me_v5 in Section 4: prev_n_tested counts DRUG-tests,
    # and one colony assayed against the panel contributes ~8.5 of them. Those
    # calls are correlated within a colony, so they are not 8.5 independent
    # draws. The composite uses its OWN rho and k - measured across the
    # well-populated 9-drug panel rather than the 3 aminoglycosides - because
    # correlation across classes is weaker than within one class
    # (rho ~ 0.41 vs 0.57), and k is larger.
    n_eff_c = prev_n_tested / .deff_c,
    k_eff_c = prev_n_resistant * n_eff_c / prev_n_tested,
    prev_theta_pm  = (1 + k_eff_c) / (2 + n_eff_c),
    prev_theta_psd = sqrt(
      (1 + k_eff_c) * (1 + n_eff_c - k_eff_c) /
        ((2 + n_eff_c)^2 * (3 + n_eff_c))
    ),
    prev_theta_pm_drugtests  = (1 + prev_n_resistant) / (2 + prev_n_tested),
    prev_theta_psd_drugtests = sqrt(
      (1 + prev_n_resistant) *
        (1 + prev_n_tested - prev_n_resistant) /
        ((2 + prev_n_tested)^2 * (3 + prev_n_tested))
    )
  )

cat(sprintf("\n  Composite design effect: rho = %.2f, k = %.1f drugs/colony, deff = %.2f\n",
            .rho_c, .k_c, .deff_c))
cat(sprintf("  median prev_theta_psd: %.3f (drug-tests) -> %.3f (effective tests)\n",
            median(me_comp_df$prev_theta_psd_drugtests), median(me_comp_df$prev_theta_psd)))

cat("\n3C composite ME frame: N =", nrow(me_comp_df),
    "visits,", n_distinct(me_comp_df$individual_id), "patients",
    "| dropped from route_df:", nrow(route_df) - nrow(me_comp_df), "\n")

# Coding check: does the binary prior-R used by the naive forest match the
# counts used to build theta? (catches any strict-vs-conservative mismatch)
cat("Coding check  prev_any_resistance == (prev_n_resistant > 0):",
    round(100 * mean(me_comp_df$prev_any_resistance ==
                       as.integer(me_comp_df$prev_n_resistant > 0)), 1), "%\n")

# Naive composite, refit on the SAME rows as the ME model (matched comparison).
m2_naive_comp <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance +
    prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_aminoglycoside_prescribed + prev_carbapenem_prescribed +
    prev_cephalosporin_prescribed + prev_fluoroquinolones_prescribed +
    prev_colistin_prescribed +
    (1 | individual_id),
  data = me_comp_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/m2_naive_comp", file_refit = "on_change"
)

# Noise-aware composite: prior resistance -> latent theta with known SD.
m2_me_comp <- brm(
  n_resistant | trials(n_tested) ~
    me(prev_theta_pm, prev_theta_psd) +
    prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_aminoglycoside_prescribed + prev_carbapenem_prescribed +
    prev_cephalosporin_prescribed + prev_fluoroquinolones_prescribed +
    prev_colistin_prescribed +
    (1 | individual_id),
  data = me_comp_df, family = beta_binomial(),
  prior = c(
    prior(normal(0, 1.5), class = "b"),
    prior(normal(0, 2),   class = "Intercept"),
    prior(exponential(1), class = "sd"),
    prior(gamma(2, 0.1),  class = "phi"),
    prior(normal(0.4, 0.3), class = "meanme"),
    prior(exponential(2),   class = "sdme")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 13), seed = 42,
  file = "models/m2_me_comp", file_refit = "on_change"
)

print(summary(m2_naive_comp))
print(summary(m2_me_comp))

# --- pull both posteriors; map the prior-R param in each model to one label --
comp_lab <- c(
  prev_any_resistance              = "Prior resistance",
  meprev_theta_pm                  = "Prior resistance",
  prev_has_mrsa                    = "MRSA",
  prev_has_mssa                    = "MSSA",
  prev_has_any_mucoid              = "Mucoid",
  prev_aminoglycoside_prescribed   = "Aminoglycoside",
  prev_carbapenem_prescribed       = "Carbapenem",
  prev_cephalosporin_prescribed    = "Cephalosporin",
  prev_fluoroquinolones_prescribed = "Fluoroquinolone",
  prev_colistin_prescribed         = "Colistin"
)

# --- pull both posteriors; map the prior-R param in each model to one label --
comp_lab <- c(
  prev_any_resistance              = "Prior resistance",
  meprev_theta_pm                  = "Prior resistance",
  prev_has_mrsa                    = "MRSA",
  prev_has_mssa                    = "MSSA",
  prev_has_any_mucoid              = "Mucoid",
  prev_aminoglycoside_prescribed   = "Aminoglycoside",
  prev_carbapenem_prescribed       = "Carbapenem",
  prev_cephalosporin_prescribed    = "Cephalosporin",
  prev_fluoroquinolones_prescribed = "Fluoroquinolone",
  prev_colistin_prescribed         = "Colistin"
)


comp_compare <- bind_rows(
  grab_comp(m2_naive_comp, "Naive (binary prior R)"),
  grab_comp(m2_me_comp,    "Aware (latent theta)")
)

# Survivors table (pd > 0.975) under each model = the CCM shortlist.
comp_survivors <- comp_compare %>%
  group_by(model, pretty_var) %>%
  summarise(est = median(value),
            lo  = quantile(value, .025),
            hi  = quantile(value, .975),
            pd  = max(mean(value > 0), mean(value < 0)),
            .groups = "drop") %>%
  mutate(direction = if_else(est > 0, "+", "-"),
         survives  = pd > 0.975) %>%
  arrange(pretty_var, model)

cat("\n=== Composite survivors: naive vs aware (pd > 0.975) ===\n")
print(comp_survivors, n = 40)
write_csv(comp_survivors, "results/noise/composite_naive_vs_aware_survivors.csv")

# --- figure: the full forest under each model, side by side ------------------
drug_terms_c <- c("Aminoglycoside","Carbapenem","Cephalosporin",
                  "Fluoroquinolone","Colistin")
comp_plot_df <- comp_compare %>%
  group_by(model, pretty_var) %>%
  mutate(med = median(value),
         direction   = if_else(med > 0, "+", "-"),
         pd          = max(mean(value > 0), mean(value < 0)),
         significant = pd > 0.975) %>%
  ungroup() %>%
  mutate(group = if_else(pretty_var %in% drug_terms_c,
                         "Antibiotic exposure (prior)", "Prior state & ecology"),
         group = factor(group, levels = c("Prior state & ecology",
                                          "Antibiotic exposure (prior)")),
         pretty_var = fct_reorder(pretty_var, med),
         model = factor(model, levels = c("Naive (binary prior R)",
                                          "Aware (latent theta)")))

fig_aware_forest <- ggplot(comp_plot_df, aes(value, pretty_var)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  stat_halfeye(aes(fill = direction, alpha = significant),
               .width = c(0.66, 0.95), slab_color = NA, point_size = 1.8,
               interval_size_range = c(0.5, 1.4)) +
  facet_grid(group ~ model, scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c("+" = pal[["red"]], "-" = pal[["blue"]]), guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.35), guide = "none") +
  labs(x = "Effect on resistance (log-odds)", y = NULL,
       title = "Same predictors, naive vs noise-aware",
       caption = lab_n(me_comp_df)) +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"),
        plot.title = element_text(face = "bold", size = 16))

ggsave("figures/talk/fig_aware_forest.pdf", fig_aware_forest,
       width = 11, height = 6.5, dpi = 300, bg = "transparent")

cat("\n3C done: fig_aware_forest.pdf + composite_naive_vs_aware_survivors.csv\n")


# =============================================================================
# 5.2  Per-class models  -  standard (constant phi) + deep-dive (mu/phi)
# One template, fit to each outcome class (AG, carbapenem, ceph, FQ):
#   (A) standard beta-binomial  -> cross-resistance matrix from prescriptions
#   (B) distributional mu/phi    -> which class differs in LEVEL and SPREAD,
#                                    to justify CCM / Markov targets.
# Predictors (both): prior SAME-class resistance (persistence) + ecology
#   (MRSA, MSSA, mucoid) + ALL FIVE prescriptions (same-class = selection,
#   cross-class = cross-resistance). Colistin is predictor-only.
# Depends on: pa_visit_full, priors_beta_random, priors_dist_betabin, pal, lab_n
# =============================================================================
dir.create("models/perclass",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/perclass",  recursive = TRUE, showWarnings = FALSE)
dir.create("figures/perclass",  recursive = TRUE, showWarnings = FALSE)

class_cfg <- tibble::tribble(
  ~cls,   ~label,            ~nres,              ~ntest,          ~persist,
  "ag",   "Aminoglycoside",  "ag_n_resistant",   "ag_n_tested",   "prev_ag_any_resistant",
  "carb", "Carbapenem",      "carb_n_resistant", "carb_n_tested", "prev_carb_any_resistant",
  "ceph", "Cephalosporin",   "ceph_n_resistant", "ceph_n_tested", "prev_ceph_any_resistant",
  "fq",   "Fluoroquinolone", "fq_n_resistant",   "fq_n_tested",   "prev_fq_any_resistant"
)

common_rhs <- paste(
  "prev_has_any_mucoid + prev_has_mrsa + prev_has_mssa +",
  "prev_aminoglycoside_prescribed + prev_carbapenem_prescribed +",
  "prev_cephalosporin_prescribed + prev_fluoroquinolones_prescribed +",
  "prev_colistin_prescribed"
)
common_vars <- c("prev_has_any_mucoid","prev_has_mrsa","prev_has_mssa",
                 "prev_aminoglycoside_prescribed","prev_carbapenem_prescribed",
                 "prev_cephalosporin_prescribed","prev_fluoroquinolones_prescribed",
                 "prev_colistin_prescribed")
common_labels <- c("Mucoid","MRSA","MSSA","AG Rx","Carbapenem Rx",
                   "Cephalosporin Rx","Fluoroquinolone Rx","Colistin Rx")


class_fits <- lapply(seq_len(nrow(class_cfg)), function(i) fit_one_class(class_cfg[i, ]))
names(class_fits) <- class_cfg$cls

# --- extract effects (summary) -> long tables --------------------------------
naive_eff <- bind_rows(lapply(class_fits, function(f) eff_from(f$m_naive, f$label, f$persist, "naive")))
dist_eff  <- bind_rows(lapply(class_fits, function(f) eff_from(f$m_dist,  f$label, f$persist, "dist")))
write_csv(naive_eff, "results/perclass/naive_effects.csv")
write_csv(dist_eff,  "results/perclass/dist_effects.csv")

cat("\n=== Per-class mu survivors (pd > 0.975) - CCM shortlist ===\n")
dist_eff %>% filter(submodel == "mu", survives) %>%
  arrange(class, desc(abs(est))) %>% select(class, pred, est, pd) %>%
  print(n = 60)

# --- FIGURE 1: cross-resistance matrix (standard naive mu, Rx block) ----------
rx_levels  <- c("AG Rx","Carbapenem Rx","Cephalosporin Rx","Fluoroquinolone Rx","Colistin Rx")
diag_key   <- c("AG Rx" = "Aminoglycoside", "Carbapenem Rx" = "Carbapenem",
                "Cephalosporin Rx" = "Cephalosporin", "Fluoroquinolone Rx" = "Fluoroquinolone")
xres_mat <- naive_eff %>%
  filter(submodel == "mu", pred %in% rx_levels) %>%
  mutate(class   = factor(class, levels = class_cfg$label),
         pred    = factor(pred,  levels = rev(rx_levels)),
         lab     = sprintf("%.2f%s", est, if_else(survives, "*", "")),
         txtcol  = if_else(abs(est) > 0.6, "white", "grey15"),
         is_diag = mapply(function(p, c) identical(unname(diag_key[as.character(p)]), as.character(c)),
                          pred, class))

fig_xres_matrix <- ggplot(xres_mat, aes(class, pred)) +
  geom_tile(aes(fill = est), color = "grey90") +
  geom_tile(data = ~ filter(.x, is_diag), fill = NA, color = pal[["ink"]], linewidth = 1.1) +
  geom_text(aes(label = lab, color = txtcol), size = 4.2) +
  scale_color_identity() +
  scale_fill_gradient2(low = pal[["blue"]], mid = "white", high = pal[["red"]],
                       midpoint = 0, name = "log-OR", breaks = c(-0.2, 0.2, 0.6)) +
  labs(x = "Resistance outcome (class)", y = "Prior prescription",
       title = "Cross-resistance: prescriptions vs class-specific resistance",
       caption = "Boxed = same-class (direct selection); * = 95% CI excludes 0. Standard beta-binomial, per class.") +
  theme(plot.title = element_text(face = "bold", size = 15), panel.grid = element_blank())
ggsave("figures/perclass/fig_xres_matrix.pdf", fig_xres_matrix,
       width = 11.5, height = 9.5, dpi = 300, bg = "transparent")

# --- FIGURES 2 & 3: per-class deep dive (distributional mu, then phi) ---------
pred_order <- c("Prior same-class R","Mucoid","MRSA","MSSA",
                "AG Rx","Carbapenem Rx","Cephalosporin Rx","Fluoroquinolone Rx","Colistin Rx")
fig_perclass_mu  <- perclass_forest("mu",
                                    "Effect on resistance level (log-odds)",
                                    "Per-class deep dive: which predictors survive, by class")
fig_perclass_phi <- perclass_forest("phi",
                                    "Effect on dispersion (log phi; + = more homogeneous)",
                                    "Per-class deep dive: what drives heterogeneity, by class")
ggsave("figures/perclass/fig_perclass_mu.pdf",  fig_perclass_mu,  width = 13, height = 5, dpi = 300, bg = "transparent")
ggsave("figures/perclass/fig_perclass_phi.pdf", fig_perclass_phi, width = 13, height = 5, dpi = 300, bg = "transparent")

cat("\n3D done: fig_xres_matrix + fig_perclass_mu + fig_perclass_phi; naive/dist_effects.csv\n")


# =============================================================================
# 5.3  Clinical-status-aware per-class mu model
# Adds current-visit clinical_status (Stable / APE Outpatient / APE IV) to the
# mu submodel of each per-class distributional model. Answers "did you adjust
# for how sick / hospitalized the patient was?" - APE IV Antibiotics is the
# sickest tier (IV abx for an exacerbation), the closest proxy to hospitalized.
# Test: do the diffuse prescription effects ATTENUATE (= severity confounding)
# while persistence HOLDS? Parallels the binary-vs-route comparison.
# Depends on Section 3D objects (class_cfg, common_vars/labels/rhs,
#   priors_dist_betabin, dist_eff, pal).
# =============================================================================
sev_fits <- lapply(seq_len(nrow(class_cfg)), function(i) fit_class_sev(class_cfg[i, ]))

# --- extract mu effects incl. clinical_status terms --------------------------
sev_eff <- bind_rows(lapply(sev_fits, function(f) eff_sev(f$m, f$label, f$persist)))
write_csv(sev_eff, "results/perclass/clinical_status_aware_effects.csv")

cat("\n=== Clinical-status severity effects (vs Stable), by class ===\n")
sev_eff %>% filter(str_detect(pred, "^Status:")) %>%
  arrange(class, pred) %>% select(class, pred, est, lo, hi, pd) %>% print(n = 20)

# --- base vs +clinical-status comparison on persistence + prescriptions -------
key_order <- c("Prior same-class R","AG Rx","Carbapenem Rx","Cephalosporin Rx",
               "Fluoroquinolone Rx","Colistin Rx","Status: APE Outpatient","Status: APE IV")
base_cmp <- dist_eff %>% filter(submodel == "mu", pred %in% key_order) %>%
  transmute(class, pred, est, lo, hi, pd, survives, model = "Base")
sev_cmp  <- sev_eff %>% filter(pred %in% key_order) %>%
  transmute(class, pred, est, lo, hi, pd, survives, model = "+ clinical status")
cmp <- bind_rows(base_cmp, sev_cmp) %>%
  mutate(pred  = factor(pred, levels = rev(key_order)),
         class = factor(class, levels = class_cfg$label),
         model = factor(model, levels = c("Base", "+ clinical status")))
write_csv(cmp, "results/perclass/base_vs_clinical_status.csv")

cat("\n=== Persistence + prescription shift when adjusting for clinical status ===\n")
cmp %>% filter(!str_detect(pred, "^Status:")) %>%
  select(class, pred, model, est, pd) %>%
  tidyr::pivot_wider(names_from = model, values_from = c(est, pd)) %>%
  print(n = 40)

fig_status_aware <- ggplot(cmp, aes(est, pred, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.55),
                  linewidth = 0.6, size = 0.35, na.rm = TRUE) +
  facet_wrap(~ class, nrow = 1) +
  scale_color_manual(values = c("Base" = pal[["sand"]], "+ clinical status" = pal[["ink"]]),
                     name = NULL) +
  labs(x = "Effect on resistance level (log-odds)", y = NULL,
       title = "Clinical-status-aware: do prescription effects survive severity adjustment?",
       caption = "Severity tiers vs Clinically Stable; APE IV = sickest (hospitalized-equivalent).") +
  theme(plot.title = element_text(face = "bold", size = 14),
        strip.text = element_text(face = "bold"), legend.position = "top")
ggsave("figures/perclass/fig_status_aware.pdf", fig_status_aware,
       width = 13, height = 5.5, dpi = 300, bg = "transparent")

cat("\n3E done: fig_status_aware + clinical_status_aware_effects.csv + base_vs_clinical_status.csv\n")



# #############################################################################
# SECTION 6: ECOLOGY  -  co-colonizers and clinical context
# 
# #############################################################################


# -----------------------------------------------------------------------------
# 6.1  Exacerbation vs route  [Slide 6]
# 
# -----------------------------------------------------------------------------

# prev_has_exacerbation is a VALID severity/context adjustment on the full
# ecological subset. Route (here: total route-weighted pressure) is a different
# construct. Fit both on the same rows and compare the CO-COLONIZER coefficients:
#   stable across adjustments -> swap is safe (route ok for the noise theme)
#   coefficients move          -> exacerbation was adjusting for severity that
#                                 route does not capture -> KEEP exacerbation.
# =============================================================================

eco_shared <- c("prev_any_resistance", "prev_has_mrsa", "prev_has_mssa",
                "prev_has_any_mucoid")
eco_df <- pa_ecological %>%
  filter(if_all(all_of(c("n_resistant", "n_tested", "individual_id",
                         eco_shared, "prev_has_exacerbation",
                         "prev_total_antibiotic_pressure")), ~ !is.na(.))) %>%
  mutate(prev_total_pressure_z = as.numeric(scale(prev_total_antibiotic_pressure)))

cat("\n7.3 ecology comparison frame: N =", nrow(eco_df), "visits\n")

m1_exac <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance + prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_has_exacerbation + (1 | individual_id),
  data = eco_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/m1_eco_exac_common", file_refit = "on_change"
)

m1_route <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance + prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_total_pressure_z + (1 | individual_id),
  data = eco_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/m1_eco_route", file_refit = "on_change"
)

cocolonizer_compare <- bind_rows(
  as_draws_df(m1_exac)  %>% transmute(model = "exacerbation",
                                      MRSA = b_prev_has_mrsa,
                                      MSSA = b_prev_has_mssa,
                                      mucoid = b_prev_has_any_mucoid),
  as_draws_df(m1_route) %>% transmute(model = "route_total",
                                      MRSA = b_prev_has_mrsa,
                                      MSSA = b_prev_has_mssa,
                                      mucoid = b_prev_has_any_mucoid)
) %>%
  pivot_longer(c(MRSA, MSSA, mucoid), names_to = "cocolonizer") %>%
  group_by(model, cocolonizer) %>%
  summarise(median = median(value),
            lo = quantile(value, .025),
            hi = quantile(value, .975), .groups = "drop")

cat("\nCo-colonizer effects under each adjustment (compare for stability):\n")
print(cocolonizer_compare)
write_csv(cocolonizer_compare, "results/noise/ecology_exac_vs_route.csv")

# =============================================================================

# -----------------------------------------------------------------------------
# 6.2  Talk figures: co-colonizer effects, staph backup, time trend
# The eco_df / m1_exac / m1_route block that opened this section in
# Talk_figures_v2 is dropped - Section 6.1 above already fits those models
# and writes results/noise/ecology_exac_vs_route.csv. Figures only here.
# -----------------------------------------------------------------------------

# Co-colonizer effects stable whether adjusted for exacerbation or route total ->
# the swap is safe. Staph effects are class-specific (see fig4, poster).
# Section 6 figure: ecology exac vs route forest
# --- 3.2 Ecology: exacerbation vs route forest [Slide 6] ---------------------
# cocolonizer_compare comes from Synthesis.R Section 7.3 (median/lo/hi by model)
fig_ecology <- ggplot(cocolonizer_compare,
                      aes(median, cocolonizer, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.5), linewidth = 0.8) +
  scale_color_manual(values = c("exacerbation" = "#B85C44", "route_total" = "#207E9D"),
                     name = "Adjustment") +
  labs(x = "Co-colonizer effect on resistance (log-odds)", y = NULL,
       title = "Do co-colonizer effects survive the route swap?",
       subtitle = "Stable across adjustments = swap is safe; shifted = keep exacerbation",
       caption = lab_n(eco_df)) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", size = 16))

ggsave("figures/talk/fig_ecology_exac_vs_route.pdf", fig_ecology, width = 9, height = 5, dpi = 300, bg = "transparent")

# =============================================================================
# BACKUP (Q&A): STAPH CO-COLONIZATION  -  descriptive, class-specific
# =============================================================================
# Quick-glance association to have ready for questions. The rigorous version is
# fig4 (stratified mu/phi). These are RAW proportions, not adjusted for prior
# resistance, severity, or sampling. FQ stands out (higher with MRSA, lower with
# MSSA); the other classes are roughly flat. Say that out loud if you show it.
# -----------------------------------------------------------------------------
staph_levels <- c("No staph", "MSSA", "MRSA")
res_recode <- c(ag_any_resistant_conservative   = "Aminoglycoside",
                carb_any_resistant_conservative = "Carbapenem",
                fq_any_resistant_conservative   = "Fluoroquinolone",
                ceph_any_resistant_conservative = "Cephalosporin")

staph_df <- pa_visit_full %>%
  filter(!is.na(prev_has_mrsa) | !is.na(prev_has_mssa)) %>%
  mutate(pstaph = case_when(coalesce(prev_has_mrsa, 0) == 1 ~ "MRSA",
                            coalesce(prev_has_mssa, 0) == 1 ~ "MSSA",
                            TRUE ~ "No staph"),
         pstaph = factor(pstaph, levels = staph_levels)) %>%
  pivot_longer(all_of(names(res_recode)), names_to = "class", values_to = "r") %>%
  filter(!is.na(r)) %>%
  mutate(class = recode(class, !!!res_recode)) %>%
  group_by(class, pstaph) %>%
  summarise(rate = mean(r == 1), .groups = "drop")

fig_staph_backup <- ggplot(staph_df, aes(pstaph, rate, fill = pstaph)) +
  geom_col(width = 0.7) +
  facet_wrap(~class) +
  scale_fill_manual(values = c("No staph" = pal[["sand"]],
                               "MSSA" = pal[["blue"]],
                               "MRSA" = pal[["red"]]), guide = "none") +
  scale_y_continuous(labels = scales::percent, expand = expansion(c(0, 0.08))) +
  labs(x = NULL, y = "Resistant (R + I)",
       title = "Before adjustment, FQ resistance tracks staph",
       caption = "Raw, unadjusted. The FQ gap is severity confounding.") +
  theme(panel.spacing = unit(1, "lines"))

ggsave("figures/talk/fig_staph_backup.pdf", fig_staph_backup,
       width = 9, height = 6.5, dpi = 300, bg = "transparent")

# =============================================================================
# OPTIONAL: DYNAMIC RESISTANCE OVER TIME  -  stable population, restless people
# =============================================================================
# Option 3: per-class population mean (+ ribbon) over time, on a grey wash of
# individual patient trajectories. FRAME AS STATIONARITY, not trend - the
# population line is ~flat (matches the variance decomposition: between-visit
# variance is 1-4%), while individuals bounce. That contrast IS the thesis.
# Use this instead of, or alongside, the summary panel for slide 2/3.
# -----------------------------------------------------------------------------
time_classes <- c(ag_any_resistant_conservative   = "Aminoglycoside",
                  carb_any_resistant_conservative = "Carbapenem",
                  fq_any_resistant_conservative   = "Fluoroquinolone",
                  ceph_any_resistant_conservative = "Cephalosporin")

time_long <- pa_visit_full %>%
  filter(!is.na(days_since_first_visit)) %>%
  mutate(yr = days_since_first_visit / 365.25) %>%
  pivot_longer(all_of(names(time_classes)), names_to = "class", values_to = "r") %>%
  filter(!is.na(r)) %>%
  mutate(class = recode(class, !!!time_classes))

# individual trajectories (grey wash): patients with >= 3 visits in a class
traj_bg <- time_long %>%
  add_count(individual_id, class, name = "nv") %>%
  filter(nv >= 3)

# population mean + binomial CI per ~half-year bin
pop_trend <- time_long %>%
  mutate(ybin = floor(yr * 2) / 2 + 0.25) %>%
  group_by(class, ybin) %>%
  summarise(p = mean(r), n = n(), .groups = "drop") %>%
  filter(n >= 15) %>%
  mutate(se = sqrt(p * (1 - p) / n), lo = pmax(0, p - 1.96 * se), hi = pmin(1, p + 1.96 * se))

fig_time <- ggplot() +
  geom_line(data = traj_bg,
            aes(yr, r, group = interaction(individual_id, class)),
            colour = "grey75", linewidth = 0.3, alpha = 0.5) +
  geom_ribbon(data = pop_trend, aes(ybin, ymin = lo, ymax = hi),
              fill = pal[["red"]], alpha = 0.18) +
  geom_line(data = pop_trend, aes(ybin, p), colour = pal[["red"]], linewidth = 1.1) +
  facet_wrap(~class) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Years since first visit", y = "Resistant (R + I)",
       title = "Stable population, restless individuals",
       caption = lab_n(time_long)) +
  theme(panel.spacing = unit(1, "lines"))

ggsave("figures/talk/fig_resistance_over_time.pdf", fig_time,
       width = 11, height = 7, dpi = 300, bg = "transparent")


# =============================================================================

# #############################################################################
# SECTION 7: DYNAMICS  -  Markov chains, EDM/S-map, cross-mapping
# 
# #############################################################################


# -----------------------------------------------------------------------------
# 7.1  Markov chain analysis, 4 outcome classes
# 
# -----------------------------------------------------------------------------

# 5 drug classes as PREDICTORS: includes Colistin
# =============================================================================

library(tidyverse)
library(markovchain)
library(rEDM)
library(patchwork)

for (dir in c("figures", "results")) {
  dir.create(dir, showWarnings = FALSE)
}

# =============================================================================
# PART 1: MARKOV CHAIN ANALYSIS - 4 OUTCOME CLASSES
# =============================================================================

# -----------------------------------------------------------------------------
# Helper function: Calculate transition matrix with CIs
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Overall resistance transitions
# -----------------------------------------------------------------------------

overall_trans <- pa_prescription %>%
  filter(!is.na(prev_any_resistance), !is.na(any_resistance)) %>%
  mutate(
    prev_state = if_else(prev_any_resistance == 1, "R", "S"),
    current_state = if_else(any_resistance == 1, "R", "S")
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
  mutate(drug = "Overall")

cat("Overall Resistance Transitions:\n")
print(overall_trans)

# Create markovchain object
trans_matrix_overall <- overall_trans %>%
  dplyr::select(prev_state, current_state, prob) %>%
  pivot_wider(names_from = current_state, values_from = prob, values_fill = 0) %>%
  column_to_rownames("prev_state") %>%
  as.matrix()

trans_matrix_overall <- trans_matrix_overall[c("S", "R"), c("S", "R")]

mc_overall <- new("markovchain",
                  states = c("S", "R"),
                  transitionMatrix = trans_matrix_overall,
                  name = "Overall")

steady_overall <- steadyStates(mc_overall)
cat("\nSteady-state distribution (Overall):\n")
cat("  Susceptible:", round(steady_overall[1, "S"], 3), "\n")
cat("  Resistant:", round(steady_overall[1, "R"], 3), "\n\n")

# -----------------------------------------------------------------------------
# Drug-specific transitions (4 OUTCOME CLASSES - excludes Colistin)
# -----------------------------------------------------------------------------

# Define drugs - 5 as EXPOSURES, but only 4 as OUTCOMES
drug_specs <- tibble(
  rx_var = c("prev_carbapenem_prescribed", "prev_aminoglycoside_prescribed",
             "prev_fluoroquinolones_prescribed", "prev_cephalosporin_prescribed"),
  resist_prev = c("prev_carb_any_resistant", "prev_ag_any_resistant",
                  "prev_fq_any_resistant", "prev_ceph_any_resistant"),
  resist_curr = c("carb_resistant", "ag_resistant", "fq_resistant", "ceph_resistant"),
  label = c("Carbapenem", "Aminoglycoside", "Fluoroquinolone", "Cephalosporin")
)

# Calculate transitions for each drug (exposed vs unexposed)
drug_transitions <- pmap_dfr(drug_specs, function(rx_var, resist_prev, resist_curr, label) {
  
  if (!rx_var %in% names(pa_prescription)) {
    cat("Skipping", label, "- column not found\n")
    return(NULL)
  }
  
  if (!resist_prev %in% names(pa_prescription) | !resist_curr %in% names(pa_prescription)) {
    cat("Skipping", label, "- resistance columns not found\n")
    return(NULL)
  }
  
  exposed <- pa_prescription %>% filter(.data[[rx_var]] == 1)
  unexposed <- pa_prescription %>% filter(.data[[rx_var]] == 0)
  
  cat(label, ":\n")
  cat("  Exposed: N =", nrow(exposed), "\n")
  cat("  Unexposed: N =", nrow(unexposed), "\n")
  
  if (nrow(exposed) < 20 | nrow(unexposed) < 20) {
    cat("  Skipping - insufficient data\n\n")
    return(NULL)
  }
  
  bind_rows(
    calculate_transitions(exposed, resist_prev, resist_curr, paste0(label, " (Yes)")),
    calculate_transitions(unexposed, resist_prev, resist_curr, paste0(label, " (No)"))
  ) %>%
    mutate(drug_class = label)
})

# Flag low-N results
drug_transitions <- drug_transitions %>%
  mutate(
    low_n_warning = if_else(total < 20, "⚠ Low N", ""),
    reliable = total >= 20
  )

# -----------------------------------------------------------------------------
# Colistin EXPOSURE effect on OTHER class transitions
# -----------------------------------------------------------------------------

colistin_exposure_transitions <- pmap_dfr(drug_specs, function(rx_var, resist_prev, resist_curr, label) {
  
  coli_exposed <- pa_prescription %>% filter(prev_colistin_prescribed == 1)
  coli_unexposed <- pa_prescription %>% filter(prev_colistin_prescribed == 0)
  
  cat("Colistin exposure →", label, "resistance:\n")
  cat("  Colistin exposed: N =", nrow(coli_exposed), "\n")
  cat("  Colistin unexposed: N =", nrow(coli_unexposed), "\n")
  
  if (nrow(coli_exposed) < 20 | nrow(coli_unexposed) < 20) {
    cat("  Skipping - insufficient data\n\n")
    return(NULL)
  }
  
  bind_rows(
    calculate_transitions(coli_exposed, resist_prev, resist_curr, paste0("Coli→", label, " (Yes)")),
    calculate_transitions(coli_unexposed, resist_prev, resist_curr, paste0("Coli→", label, " (No)"))
  ) %>%
    mutate(
      exposure = "Colistin",
      outcome_class = label
    )
})

# -----------------------------------------------------------------------------
# Summary: Acquisition and Reversion rates
# -----------------------------------------------------------------------------

transition_summary <- drug_transitions %>%
  filter(current_state != prev_state) %>%
  mutate(
    transition_type = if_else(prev_state == "S", "Acquisition (S→R)", "Reversion (R→S)"),
    exposed = if_else(str_detect(drug, "Yes"), "Exposed", "Unexposed")
  ) %>%
  dplyr::select(drug_class, exposed, transition_type, n, total, prob, ci_lower, ci_upper)

cat("\nTransition Summary (4 classes):\n")
print(transition_summary, n = 40)

# Calculate risk ratios
risk_ratios <- transition_summary %>%
  dplyr::select(drug_class, exposed, transition_type, prob) %>%
  pivot_wider(names_from = exposed, values_from = prob) %>%
  mutate(
    risk_ratio = Exposed / Unexposed,
    interpretation = case_when(
      risk_ratio > 1.2 & transition_type == "Acquisition (S→R)" ~ "↑ Acquisition",
      risk_ratio < 0.8 & transition_type == "Acquisition (S→R)" ~ "↓ Acquisition",
      risk_ratio > 1.2 & transition_type == "Reversion (R→S)" ~ "↑ Reversion",
      risk_ratio < 0.8 & transition_type == "Reversion (R→S)" ~ "↓ Reversion",
      TRUE ~ "No strong effect"
    )
  )

cat("\nRisk Ratios (Exposed / Unexposed) - 4 classes:\n")
print(risk_ratios)

write_csv(risk_ratios, "results/markov_risk_ratios_4classes.csv")

# Colistin cross-effect risk ratios
if (nrow(colistin_exposure_transitions) > 0) {
  coli_risk_ratios <- colistin_exposure_transitions %>%
    filter(current_state != prev_state) %>%
    mutate(
      transition_type = if_else(prev_state == "S", "Acquisition (S→R)", "Reversion (R→S)"),
      exposed = if_else(str_detect(drug, "Yes"), "Coli Exposed", "Coli Unexposed")
    ) %>%
    dplyr::select(outcome_class, exposed, transition_type, prob) %>%
    pivot_wider(names_from = exposed, values_from = prob) %>%
    mutate(
      risk_ratio = `Coli Exposed` / `Coli Unexposed`,
      interpretation = case_when(
        risk_ratio > 1.2 & transition_type == "Acquisition (S→R)" ~ "↑ Acquisition",
        risk_ratio < 0.8 & transition_type == "Acquisition (S→R)" ~ "↓ Acquisition",
        TRUE ~ "No strong effect"
      )
    )
  
  cat("\nColistin Cross-Effect Risk Ratios:\n")
  print(coli_risk_ratios)
  
  write_csv(coli_risk_ratios, "results/markov_colistin_cross_effects.csv")
}

# -----------------------------------------------------------------------------
# Steady-state analysis by drug exposure (4 classes)
# -----------------------------------------------------------------------------

steady_states <- drug_transitions %>%
  dplyr::select(drug, prev_state, current_state, prob) %>%
  pivot_wider(names_from = current_state, values_from = prob, values_fill = 0) %>%
  group_by(drug) %>%
  summarise(
    P_SS = first(S[prev_state == "S"]),
    P_SR = first(R[prev_state == "S"]),
    P_RS = first(S[prev_state == "R"]),
    P_RR = first(R[prev_state == "R"]),
    .groups = "drop"
  ) %>%
  mutate(
    steady_R = P_SR / (P_SR + P_RS),
    steady_S = 1 - steady_R
  )

cat("\nSteady-state resistant proportion by drug exposure (4 classes):\n")
print(steady_states %>% dplyr::select(drug, P_SR, P_RS, steady_R) %>% arrange(desc(steady_R)))

write_csv(steady_states, "results/markov_steady_states_4classes.csv")

# -----------------------------------------------------------------------------
# Figure: Transition probabilities by drug (4 classes)
# -----------------------------------------------------------------------------

fig_transitions_drug <- drug_transitions %>%
  filter(prev_state != current_state) %>%
  filter(reliable) %>%
  mutate(
    transition = paste0(prev_state, "→", current_state),
    exposed = if_else(str_detect(drug, "Yes"), "Exposed", "Unexposed"),
    drug_class = factor(drug_class, levels = c("Carbapenem", "Aminoglycoside", 
                                               "Fluoroquinolone", "Cephalosporin"))
  ) %>%
  ggplot(aes(x = drug_class, y = prob, fill = exposed)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(0.8), width = 0.2) +
  geom_text(aes(label = paste0(round(prob*100, 1), "%")),
            position = position_dodge(0.8), vjust = -0.5, size = 2.5) +
  facet_wrap(~transition, scales = "free_y") +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Exposed" = "#E41A1C", "Unexposed" = "#377EB8")) +
  labs(
    title = "Resistance State Transitions by Antibiotic Exposure",
    subtitle = "4 outcome classes (Colistin excluded due to insufficient testing data)",
    x = NULL, y = "Transition Probability", fill = "Same-Class Exposure"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave("figures/markov/markov_drug_transitions_4classes.pdf", fig_transitions_drug, 
       width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Figure: Steady-state comparison (4 classes)
# -----------------------------------------------------------------------------

fig_steady <- steady_states %>%
  mutate(
    drug_class = str_remove(drug, " \\(Yes\\)| \\(No\\)"),
    exposed = if_else(str_detect(drug, "Yes"), "Exposed", "Unexposed")
  ) %>%
  filter(!is.na(steady_R)) %>%
  ggplot(aes(x = reorder(drug_class, steady_R), y = steady_R, fill = exposed)) +
  geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = steady_overall[1, "R"], linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 0.5, y = steady_overall[1, "R"] + 0.03, 
           label = paste0("Overall: ", round(steady_overall[1, "R"]*100, 1), "%"),
           hjust = 0, size = 3, color = "red") +
  geom_text(aes(label = paste0(round(steady_R*100, 1), "%")),
            position = position_dodge(0.8), vjust = -0.5, size = 3) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(values = c("Exposed" = "#E41A1C", "Unexposed" = "#377EB8")) +
  labs(
    title = "Steady-State Resistance by Drug Exposure",
    subtitle = "4 outcome classes | Equilibrium if current prescribing patterns continued",
    x = NULL, y = "Steady-State Resistant (%)", fill = "Drug Exposure"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  coord_flip()

ggsave("figures/markov/markov_steady_states_4classes.pdf", fig_steady, 
       width = 10, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# 7.2  Talk figures: transition dynamics
# 
# -----------------------------------------------------------------------------

# Which associations survive once the system is treated as dynamic.
# --- 4-class Markov: visit-to-visit transitions, computed fresh ---------------
# 2-state (Susceptible / Resistant) per class on the conservative (R+I) binary,
# over consecutive visits (prev_* already lag within patient). Recomputed here
# so nothing carries over from earlier single-class runs.
markov_cfg <- tibble::tribble(
  ~class,            ~prev,                     ~cur,
  "Aminoglycoside",  "prev_ag_any_resistant",   "ag_any_resistant_conservative",
  "Carbapenem",      "prev_carb_any_resistant", "carb_any_resistant_conservative",
  "Cephalosporin",   "prev_ceph_any_resistant", "ceph_any_resistant_conservative",
  "Fluoroquinolone", "prev_fq_any_resistant",   "fq_any_resistant_conservative")

markov <- purrr::pmap_dfr(markov_cfg, markov_one) %>%
  mutate(class = factor(class, levels = markov_cfg$class))

cat("\n=== 4-class Markov (conservative R+I, visit-to-visit) ===\n")
markov %>% mutate(across(where(is.numeric), ~round(., 2))) %>% print()
write_csv(markov, "results/noise/markov_4classes.csv")

# Figure 1: per-class 2-state transition diagrams (acquire / lose / persist) ---
mk_nodes <- markov %>% distinct(class) %>%
  mutate(d = list(tibble(state = c("Susceptible", "Resistant"),
                         x = c(0, 4), y = c(0, 0)))) %>% tidyr::unnest(d)
mk_lab <- markov %>% transmute(class,
                               acq   = sprintf("acquire %.0f%%", 100 * p_SR),
                               los   = sprintf("lose %.0f%%",    100 * p_RS),
                               stayS = sprintf("stay %.0f%%",    100 * p_SS),
                               stayR = sprintf("stay %.0f%%",    100 * p_RR))

fig_markov_trans <- ggplot() +
  geom_curve(data = markov, aes(x = 0.7, y = 0.35, xend = 3.3, yend = 0.35),
             curvature = -0.3, arrow = arrow(length = unit(0.18, "cm")), color = pal[["red"]]) +
  geom_curve(data = markov, aes(x = 3.3, y = -0.35, xend = 0.7, yend = -0.35),
             curvature = -0.3, arrow = arrow(length = unit(0.18, "cm")), color = pal[["blue"]]) +
  geom_point(data = mk_nodes, aes(x, y, fill = state), size = 37, shape = 21,
             color = "grey30", show.legend = FALSE) +
  geom_text(data = mk_nodes, aes(x, y, label = state), fontface = "bold", size = 5) +
  geom_text(data = mk_lab, aes(x = 2, y = 0.95,  label = acq),   color = pal[["red"]],  size = 6) +
  geom_text(data = mk_lab, aes(x = 2, y = -0.95, label = los),   color = pal[["blue"]], size = 6) +
  geom_text(data = mk_lab, aes(x = 0, y = 1.15,  label = stayS), color = "grey45", size = 6) +
  geom_text(data = mk_lab, aes(x = 4, y = 1.15,  label = stayR), color = "grey45", size = 6) +
  scale_fill_manual(values = c("Susceptible" = "#A9C7CC", "Resistant" = "#E3B5A8")) +
  facet_wrap(~ class) +
  coord_cartesian(xlim = c(-1, 5), ylim = c(-1.3, 1.45)) +
  labs(title = sprintf("Once resistant, it stays: persistence is %.0f-%.0f%% every class",
                       100 * min(markov$p_RR), 100 * max(markov$p_RR))) +
  theme_void() +
  theme(strip.text = element_text(face = "bold", size = 15),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.caption = element_blank())

ggsave("figures/talk/markov_drug_transitions_4classes.pdf", fig_markov_trans,
       width = 9, height = 6.5, dpi = 300, bg = "transparent")

# Figure 2: steady-state resistance per class (vs observed) --------------------
mk_steady <- markov %>%
  transmute(class, Observed = observed, `Markov steady state` = steady_R) %>%
  pivot_longer(-class, names_to = "type", values_to = "prop") %>%
  mutate(type = factor(type, levels = c("Observed", "Markov steady state")))

fig_markov_steady <- ggplot(mk_steady, aes(class, prop, fill = type)) +
  geom_col(position = position_dodge(0.72), width = 0.66) +
  geom_text(aes(label = scales::percent(prop, accuracy = 1)),
            position = position_dodge(0.72), vjust = -0.4, size = 3.4) +
  scale_fill_manual(values = c("Observed" = pal[["sand"]],
                               "Markov steady state" = pal[["red"]]), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.92)) +
  labs(x = NULL, y = "Resistant (R + I)",
       title = "Patients already sit near the long-run level",
       caption = "Steady state = acquire / (acquire + lose). Class differences are driven by acquisition, not reversion.") +
  theme_talk() +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 14))

ggsave("figures/talk/markov_steady_states_4classes.pdf", fig_markov_steady,
       width = 8.5, height = 5, dpi = 300, bg = "transparent")

# =============================================================================

# -----------------------------------------------------------------------------
# 7.3  EDM data prep and optimal E / theta
# 
# -----------------------------------------------------------------------------

# =============================================================================

edm_data <- pa_prescription %>%
  arrange(individual_id, date_collected_from_subject) %>%
  mutate(
    time = row_number(),
    
    # Outcomes - overall and 4 class-specific (NO colistin)
    resistance = prop_resistant_conservative,
    carb_resist = as.numeric(carb_resistant),
    ag_resist = as.numeric(ag_resistant),
    fq_resist = as.numeric(fq_resistant),
    ceph_resist = as.numeric(ceph_resistant),
    # NOTE: coli_resist NOT included as outcome
    
    # Antibiotics - ALL 5 as predictors (including colistin)
    carbapenem = replace_na(prev_carbapenem_prescribed, 0),
    aminoglycoside = replace_na(prev_aminoglycoside_prescribed, 0),
    fluoroquinolone = replace_na(prev_fluoroquinolones_prescribed, 0),
    cephalosporin = replace_na(prev_cephalosporin_prescribed, 0),
    colistin = replace_na(prev_colistin_prescribed, 0),  # KEPT as predictor
    
    # Ecological factors
    has_mrsa = replace_na(prev_has_mrsa, 0),
    has_mssa = replace_na(prev_has_mssa, 0),
    has_mucoid = replace_na(prev_has_any_mucoid, 0)
    # NOTE: exacerbation intentionally OMITTED here. edm_data is built from
    # pa_prescription, which is entirely APE-at-t-1 visits, so the variable
    # is constant (= 1) and would carry zero information for CCM or S-map.
  ) %>%
  dplyr::select(time, individual_id,
                resistance, carb_resist, ag_resist, fq_resist, ceph_resist,  # 5 outcomes
                carbapenem, aminoglycoside, fluoroquinolone, cephalosporin, colistin,  # 5 predictors
                has_mrsa, has_mssa, has_mucoid) %>%
  filter(!is.na(resistance))

edm_df <- as.data.frame(edm_data)
n <- nrow(edm_df)

cat("  N observations:", n, "\n")
cat("  N patients:", n_distinct(edm_data$individual_id), "\n")
cat("\nOUTCOMES (5): resistance, carb_resist, ag_resist, fq_resist, ceph_resist\n")
cat("  NOTE: Colistin resistance EXCLUDED as outcome\n")
cat("\nPREDICTORS (5 drugs): carbapenem, aminoglycoside, fluoroquinolone, cephalosporin, colistin\n")
cat("  NOTE: Colistin INCLUDED as predictor\n")
cat("Ecological: has_mrsa, has_mssa, has_mucoid\n")
cat("Clinical: (exacerbation omitted — constant in pa_prescription)\n\n")

# Library and prediction sets
lib_end <- floor(n * 0.7)
lib_str <- paste("1", lib_end)
pred_str <- paste(lib_end + 1, n)

# =============================================================================
# PART 3: FIND OPTIMAL E AND THETA FOR EACH TARGET (5 targets)
# =============================================================================


# -----------------------------------------------------------------------------
# Find parameters for 5 targets (Overall + 4 class-specific)
# -----------------------------------------------------------------------------

cat("1. Overall Resistance\n")
params_overall <- find_optimal_params(edm_df, "resistance")
cat("   E =", params_overall$E, ", θ =", params_overall$theta, 
    ", ρ =", round(params_overall$rho_theta, 3), "\n\n")

cat("2. Carbapenem Resistance\n")
params_carb <- find_optimal_params(edm_df, "carb_resist")
cat("   E =", params_carb$E, ", θ =", params_carb$theta,
    ", ρ =", round(params_carb$rho_theta, 3), "\n\n")

cat("3. Aminoglycoside Resistance\n")
params_ag <- find_optimal_params(edm_df, "ag_resist")
cat("   E =", params_ag$E, ", θ =", params_ag$theta,
    ", ρ =", round(params_ag$rho_theta, 3), "\n\n")

cat("4. Fluoroquinolone Resistance\n")
params_fq <- find_optimal_params(edm_df, "fq_resist")
cat("   E =", params_fq$E, ", θ =", params_fq$theta,
    ", ρ =", round(params_fq$rho_theta, 3), "\n\n")

cat("5. Cephalosporin Resistance\n")
params_ceph <- find_optimal_params(edm_df, "ceph_resist")
cat("   E =", params_ceph$E, ", θ =", params_ceph$theta,
    ", ρ =", round(params_ceph$rho_theta, 3), "\n\n")

# NOTE: Colistin resistance parameter search REMOVED

# Summary table
params_summary <- tibble(
  Target = c("Overall", "Carbapenem", "Aminoglycoside", "Fluoroquinolone", "Cephalosporin"),
  N = c(params_overall$N, params_carb$N, params_ag$N, params_fq$N, params_ceph$N),
  E = c(params_overall$E, params_carb$E, params_ag$E, params_fq$E, params_ceph$E),
  theta = c(params_overall$theta, params_carb$theta, params_ag$theta, 
            params_fq$theta, params_ceph$theta),
  rho = c(params_overall$rho_theta, params_carb$rho_theta, params_ag$rho_theta,
          params_fq$rho_theta, params_ceph$rho_theta)
) %>%
  mutate(
    Dynamics = case_when(
      is.na(theta) ~ "Unknown",
      theta == 0 ~ "Linear",
      theta < 2 ~ "Weakly nonlinear",
      theta < 5 ~ "Moderately nonlinear",
      TRUE ~ "Strongly nonlinear"
    )
  )

cat("PARAMETER SUMMARY (5 targets):\n")
print(params_summary)

write_csv(params_summary, "results/edm_optimal_parameters_5targets.csv")

# -----------------------------------------------------------------------------
# Figure: E and theta by target (5 targets)
# -----------------------------------------------------------------------------

# Combine embedding data
embed_all <- bind_rows(
  if(!is.null(params_overall$embed_data)) params_overall$embed_data %>% mutate(Target = "Overall"),
  if(!is.null(params_carb$embed_data)) params_carb$embed_data %>% mutate(Target = "Carbapenem"),
  if(!is.null(params_ag$embed_data)) params_ag$embed_data %>% mutate(Target = "Aminoglycoside"),
  if(!is.null(params_fq$embed_data)) params_fq$embed_data %>% mutate(Target = "Fluoroquinolone"),
  if(!is.null(params_ceph$embed_data)) params_ceph$embed_data %>% mutate(Target = "Cephalosporin")
)

target_colors <- c(
  "Overall" = "#000000",
  "Carbapenem" = "#E41A1C",
  "Aminoglycoside" = "#377EB8",
  "Fluoroquinolone" = "#4DAF4A",
  "Cephalosporin" = "#984EA3"
)

# Combine theta data
theta_all <- bind_rows(
  if(!is.null(params_overall$theta_data)) params_overall$theta_data %>% mutate(Target = "Overall"),
  if(!is.null(params_carb$theta_data)) params_carb$theta_data %>% mutate(Target = "Carbapenem"),
  if(!is.null(params_ag$theta_data)) params_ag$theta_data %>% mutate(Target = "Aminoglycoside"),
  if(!is.null(params_fq$theta_data)) params_fq$theta_data %>% mutate(Target = "Fluoroquinolone"),
  if(!is.null(params_ceph$theta_data)) params_ceph$theta_data %>% mutate(Target = "Cephalosporin")
)

# If every find_optimal_params() call returned NULL, bind_rows() yields a
# zero-column tibble and aes(x = E) dies with "object 'E' not found" - which
# previously took the entire pipeline down at this line. Skip the panel and say
# so loudly instead; the EDM failure is the thing to fix, not the plot.
edm_params_ok <- all(c("E", "rho") %in% names(embed_all)) && nrow(embed_all) > 0 &&
                 all(c("theta", "rho") %in% names(theta_all)) && nrow(theta_all) > 0
if (!edm_params_ok) {
  cat("\n!! EDM parameter search returned no embedding/theta data -",
      "skipping figures/ccm/edm_parameters_5targets.pdf.\n")
  cat("   Check the 'EmbedDimension error:' lines above: if they say",
      "'write failed' or\n   'PDF_endpage', the graphics device broke, not the EDM.\n")
} else {

fig_embed_all <- ggplot(embed_all, aes(x = E, y = rho, color = Target)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = target_colors) +
  labs(title = "A. Optimal Embedding Dimension by Target",
       subtitle = "Peak indicates optimal E | 5 targets (Colistin excluded)",
       x = "Embedding Dimension (E)", y = "Prediction Skill (ρ)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

fig_theta_all <- ggplot(theta_all, aes(x = theta, y = rho, color = Target)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  annotate("text", x = 0.2, y = min(theta_all$rho, na.rm = TRUE), 
           label = "Linear", hjust = 0, size = 3, color = "gray50") +
  scale_color_manual(values = target_colors) +
  labs(title = "B. Nonlinearity (theta) by Target",
       subtitle = "theta = 0 is linear; higher theta = more nonlinear state-dependence",
       x = "S-map Localization (θ)", y = "Prediction Skill (ρ)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

fig_edm_params <- fig_embed_all / fig_theta_all +
  plot_annotation(title = "EDM Parameter Selection: 5 Resistance Targets")

ggsave("figures/ccm/edm_parameters_5targets.pdf", fig_edm_params,
       width = 12, height = 10, dpi = 300)

}  # end edm_params_ok


# -----------------------------------------------------------------------------
# 7.4  S-map with state-dependent effects
# 
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Prepare multivariate embedding data
# -----------------------------------------------------------------------------

# FIX (review #9): lag WITHIN patient. Ungrouped, ~8% of rows took their lag1
# from another patient's last visit and ~16% their lag2, because edm_df is 92
# patient series stacked into one frame. Grouping costs a couple of rows per
# patient and removes the contamination.
edm_df_smap <- edm_df %>%
  group_by(individual_id) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    resistance_lag1 = dplyr::lag(resistance, 1),
    resistance_lag2 = dplyr::lag(resistance, 2)
  ) %>%
  ungroup() %>%
  filter(!is.na(resistance_lag2))

n_smap <- nrow(edm_df_smap)
lib_end_smap <- floor(n_smap * 0.7)
lib_str_smap <- paste("1", lib_end_smap)
pred_str_smap <- paste(lib_end_smap + 1, n_smap)

# -----------------------------------------------------------------------------
# Run S-map with optimal theta
# Includes ALL 5 drugs as predictors (including colistin)
# -----------------------------------------------------------------------------

optimal_theta <- params_overall$theta
cat("Using theta =", optimal_theta, "from parameter optimization\n\n")

# Define columns - includes colistin as predictor
smap_columns <- "resistance resistance_lag1 resistance_lag2 carbapenem aminoglycoside fluoroquinolone cephalosporin colistin has_mrsa has_mssa has_mucoid"
E_smap <- length(strsplit(smap_columns, " ")[[1]])

cat("Running S-map with E =", E_smap, ", theta =", optimal_theta, "\n")
cat("Predictors include all 5 drug classes (including Colistin)\n\n")

smap_out <- tryCatch({
  SMap(
    dataFrame = as.data.frame(edm_df_smap),
    columns = smap_columns,
    target = "resistance",
    lib = lib_str_smap,
    pred = pred_str_smap,
    E = E_smap,
    theta = optimal_theta,
    embedded = TRUE
  )
}, error = function(e) {
  cat("S-map error:", e$message, "\n")
  NULL
})

if (!is.null(smap_out) && "coefficients" %in% names(smap_out)) {
  
  smap_coeffs <- smap_out$coefficients %>%
    as_tibble()
  
  # Rename columns to readable names
  rename_map <- c(
    "∂resistance/∂carbapenem" = "Carbapenem",
    "∂resistance/∂aminoglycoside" = "Aminoglycoside",
    "∂resistance/∂fluoroquinolone" = "Fluoroquinolone",
    "∂resistance/∂cephalosporin" = "Cephalosporin",
    "∂resistance/∂colistin" = "Colistin",  # Included!
    "∂resistance/∂has_mrsa" = "MRSA",
    "∂resistance/∂has_mssa" = "MSSA",
    "∂resistance/∂has_mucoid" = "Mucoid"
    # exacerbation dropped: constant in edm_df
  )
  
  for (old_name in names(rename_map)) {
    if (old_name %in% names(smap_coeffs)) {
      smap_coeffs <- smap_coeffs %>% rename(!!rename_map[old_name] := !!old_name)
    }
  }
  
  # Add time index
  if (!"time" %in% names(smap_coeffs)) {
    smap_coeffs <- smap_coeffs %>% mutate(time_idx = row_number())
  } else {
    smap_coeffs <- smap_coeffs %>% rename(time_idx = time)
  }
  
  # Select key variables (now includes Colistin)
  key_vars <- c("Carbapenem", "Aminoglycoside", "Fluoroquinolone", "Cephalosporin", 
                "Colistin", "MRSA", "MSSA", "Mucoid")
  available_vars <- intersect(key_vars, names(smap_coeffs))
  
  cat("Available S-map coefficients:", paste(available_vars, collapse = ", "), "\n\n")
  
  if (length(available_vars) > 0) {
    
    # Reshape to long format
    coeffs_long <- smap_coeffs %>%
      dplyr::select(time_idx, all_of(available_vars)) %>%
      pivot_longer(-time_idx, names_to = "Variable", values_to = "Coefficient") %>%
      filter(!is.na(Coefficient), !is.nan(Coefficient), is.finite(Coefficient))
    
    # Summary statistics
    coeffs_summary <- coeffs_long %>%
      group_by(Variable) %>%
      summarise(
        Mean = mean(Coefficient, na.rm = TRUE),
        SD = sd(Coefficient, na.rm = TRUE),
        CV = abs(SD / Mean) * 100,
        Pct_Positive = mean(Coefficient > 0, na.rm = TRUE) * 100,
        Min = min(Coefficient, na.rm = TRUE),
        Max = max(Coefficient, na.rm = TRUE),
        Sign_Changes = if_else(Min < 0 & Max > 0, "Yes", "No"),
        .groups = "drop"
      ) %>%
      arrange(desc(CV))
    
    cat("S-MAP COEFFICIENT SUMMARY (sorted by variability):\n")
    cat("High CV = strong state-dependence\n\n")
    print(coeffs_summary)
    
    # Highlight colistin
    cat("\nCOLISTIN COEFFICIENT:\n")
    print(coeffs_summary %>% filter(Variable == "Colistin"))
    
    write_csv(coeffs_summary, "results/smap_coefficients_summary.csv")
    
    # -------------------------------------------------------------------------
    # Figure: Time-varying coefficients (includes Colistin)
    # -------------------------------------------------------------------------
    
    smap_colors <- c(
      "Carbapenem" = "#E41A1C", "Aminoglycoside" = "#377EB8",
      "Fluoroquinolone" = "#4DAF4A", "Cephalosporin" = "#984EA3",
      "Colistin" = "#FF7F00",  # Orange for colistin
      "MRSA" = "#A65628", "MSSA" = "#F781BF",
      "Mucoid" = "#999999"
    )
    
    fig_smap_time <- ggplot(coeffs_long, aes(x = time_idx, y = Coefficient, color = Variable)) +
      geom_line(alpha = 0.3, linewidth = 0.5) +
      geom_smooth(se = TRUE, linewidth = 1.2, method = "loess", span = 0.3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      facet_wrap(~Variable, scales = "free_y", ncol = 3) +
      scale_color_manual(values = smap_colors, guide = "none") +
      labs(title = paste0("S-map Time-Varying Coefficients (theta = ", optimal_theta, ")"),
           subtitle = "All 5 drug classes as predictors (including Colistin) | High variability = state-dependent",
           x = "Time Index", y = "Local Coefficient") +
      theme_minimal(base_size = 11) +
      theme(strip.text = element_text(face = "bold"))
    
    ggsave("figures/ccm/smap_timevarying_coefficients.pdf", fig_smap_time,
           width = 12, height = 10, dpi = 300)
    
    # -------------------------------------------------------------------------
    # S-map coefficients by resistance state
    # -------------------------------------------------------------------------
    
    smap_coeffs_state <- smap_coeffs %>%
      mutate(
        resistance_state = case_when(
          edm_df_smap$resistance[time_idx] < 0.3 ~ "Low (<30%)",
          edm_df_smap$resistance[time_idx] < 0.6 ~ "Medium (30-60%)",
          TRUE ~ "High (>60%)"
        )
      )
    
    coeffs_by_state <- smap_coeffs_state %>%
      dplyr::select(resistance_state, all_of(available_vars)) %>%
      pivot_longer(-resistance_state, names_to = "Variable", values_to = "Coefficient") %>%
      filter(!is.na(Coefficient), !is.nan(Coefficient), is.finite(Coefficient)) %>%
      group_by(resistance_state, Variable) %>%
      summarise(
        Mean = mean(Coefficient, na.rm = TRUE),
        SD = sd(Coefficient, na.rm = TRUE),
        N = n(),
        .groups = "drop"
      )
    
    
    coeffs_wide <- coeffs_by_state %>%
      dplyr::select(Variable, resistance_state, Mean) %>%
      pivot_wider(names_from = resistance_state, values_from = Mean)
    print(coeffs_wide)
    
    # Colistin by state
    print(coeffs_by_state %>% filter(Variable == "Colistin"))
    
    write_csv(coeffs_by_state, "results/smap_coefficients_by_state.csv")
    
    # Figure: Coefficients by state
    fig_smap_state <- coeffs_by_state %>%
      mutate(
        resistance_state = factor(resistance_state, 
                                  levels = c("Low (<30%)", "Medium (30-60%)", "High (>60%)")),
        Variable = factor(Variable, levels = c("Carbapenem", "Aminoglycoside", "Fluoroquinolone",
                                               "Cephalosporin", "Colistin", "MRSA", "MSSA",
                                               "Mucoid"))
      ) %>%
      ggplot(aes(x = Variable, y = Mean, fill = resistance_state)) +
      geom_col(position = position_dodge(0.8), width = 0.7) +
      geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                    position = position_dodge(0.8), width = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("Low (<30%)" = "#377EB8", 
                                   "Medium (30-60%)" = "#FFD92F",
                                   "High (>60%)" = "#E41A1C")) +
      labs(title = "S-map Coefficients by Resistance State",
           subtitle = "All 5 drugs as predictors | State-dependent effects: same variable, different effect by stage",
           x = NULL, y = "Mean Local Coefficient", fill = "Resistance State") +
      theme_minimal(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "bottom")
    
    ggsave("figures/ccm/smap_coefficients_by_state.pdf", fig_smap_state,
           width = 12, height = 7, dpi = 300)
  }
}


# -----------------------------------------------------------------------------
# 7.5  CCM on endogenous state pairs
# Replaces the drug-level placeholder CCM that used to sit in Synthesis.R.
# Restored from Talk_figures.R: mucoid is back in the pair set (v2 narrowed
# vars5 -> vars4 and dropped the four mucoid pairs).
# -----------------------------------------------------------------------------

# host states with their own dynamics - NOT prescriptions (exogenous, binary,
# clinician-driven) and NOT self-persistence (that is S-map / Markov territory).
# Pairs (both directions):
#   mucoid  x {AG, carb, ceph, FQ} resistance      (the one GLMM ecology survivor)
#   resistance x resistance, all 6 class pairs      (dynamical cross-resistance)
# Uses the pilot's rEDM tooling (EmbedDimension for E, bidirectional CCM).
# Continuous class proportions (embed better than binaries); current values.
# CAVEAT: series are concatenated patient-then-date with a single time index
#   (matches the pilot); embeddings can cross patient boundaries, so treat this
#   as secondary / triangulating and read convergence conservatively.
# =============================================================================
suppressMessages(library(rEDM))
dir.create("results/ccm", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/ccm", recursive = TRUE, showWarnings = FALSE)

edm_endo <- pa_ecological %>%
  arrange(individual_id, date_collected_from_subject) %>%
  transmute(individual_id,
            ag_R   = ag_prop_resistant_conservative,
            carb_R = carb_prop_resistant_conservative,
            ceph_R = ceph_prop_resistant_conservative,
            fq_R   = fq_prop_resistant_conservative,
            mucoid = has_any_mucoid) %>%
  mutate(time = row_number()) %>%
  as.data.frame()

vars5 <- c(ag_R = "Aminoglycoside", carb_R = "Carbapenem",
           ceph_R = "Cephalosporin", fq_R = "Fluoroquinolone", mucoid = "Mucoid")

# --- optimal E per variable (simplex; capped for short series) ----------------
E_map <- sapply(names(vars5), function(v) find_E(edm_endo, v))
cat("\nCCM optimal E per variable:\n"); print(E_map)

# --- bidirectional CCM for one pair ------------------------------------------
# CCM(columns = M, target = X) cross-maps X from M's manifold => converging
#   means X drives M.  So a->b is tested with columns=b, target=a (col "b:a").

pairs <- list(
  c("mucoid","ag_R"), c("mucoid","carb_R"), c("mucoid","ceph_R"), c("mucoid","fq_R"),
  c("ag_R","carb_R"), c("ag_R","ceph_R"), c("ag_R","fq_R"),
  c("carb_R","ceph_R"), c("carb_R","fq_R"), c("ceph_R","fq_R")
)
cat("\nRunning CCM on", length(pairs), "endogenous pairs (both directions)...\n")
ccm_all <- map_dfr(pairs, function(p) run_ccm_pair(edm_endo, p[[1]], p[[2]], E_map))

# --- convergence summary: rho at smallest vs largest library -----------------
ccm_sum <- ccm_all %>%
  filter(!is.na(rho)) %>%
  group_by(driver, affected) %>%
  summarise(rho_start = rho[which.min(LibSize)],
            rho_end   = rho[which.max(LibSize)],
            delta     = rho_end - rho_start, .groups = "drop") %>%
  mutate(converges    = delta > 0.05 & rho_end > 0.1,
         driver_lab   = unname(vars5[driver]),
         affected_lab = unname(vars5[affected]))
write_csv(ccm_all, "results/ccm/ccm_curves_endogenous.csv")
write_csv(ccm_sum, "results/ccm/ccm_summary_endogenous.csv")

cat("\n=== CCM convergence (driver -> affected): converging links ===\n")
ccm_sum %>% filter(converges) %>% arrange(desc(rho_end)) %>%
  select(driver_lab, affected_lab, rho_start, rho_end, delta) %>% print(n = 30)
cat("\n(Full table in ccm_summary_endogenous.csv)\n")

# --- FIGURE: directed coupling matrix (driver -> affected) -------------------
lvl <- unname(vars5)
ccm_mat <- ccm_sum %>%
  mutate(driver_lab   = factor(driver_lab,   levels = lvl),
         affected_lab = factor(affected_lab, levels = rev(lvl)),
         txtcol = if_else(abs(rho_end) > 0.4, "white", "grey15"),
         lab    = sprintf("%.2f%s", rho_end, if_else(converges, "*", "")))
fig_ccm_endo <- ggplot(ccm_mat, aes(driver_lab, affected_lab)) +
  geom_tile(aes(fill = rho_end), color = "grey90") +
  geom_text(aes(label = lab, color = txtcol), size = 6) +
  scale_color_identity() +
  scale_fill_gradient2(low = pal[["blue"]], mid = "white", high = pal[["red"]],
                       midpoint = 0, name = "cross-map skill", breaks = c(0.2, 0.4, 0.6)) +
  labs(x = "Driver", y = "Affected",
       title = "Dynamical coupling among resistance classes") +
  theme(plot.title = element_text(face = "bold", size = 18), panel.grid = element_blank())
ggsave("figures/ccm/fig_ccm_endogenous.pdf", fig_ccm_endo,
       width = 8, height = 5.5, dpi = 300, bg = "transparent")

cat("\n5B done: fig_ccm_endogenous + ccm_summary_endogenous.csv + ccm_curves_endogenous.csv\n")


# --- FIGURE (alt): Markov-style coupling diagram -----------------------------
# Same content as the matrix, read at a glance: nodes = classes, ONE bold
# bidirectional edge for the converging pair, faint edges for tested-but-flat
# pairs. Drop-in alternative to fig_ccm_endo for slide 9. Data-driven: the bold
# link is whatever `converges` flags, so a rerun/surrogate test updates it.
# NOTE: built for the 4 resistance classes only - mucoid is dropped here
#       (binary, does not embed for CCM; also matches the deck matrix).
ccm_classes <- c("Aminoglycoside", "Carbapenem", "Cephalosporin", "Fluoroquinolone")

# node coords: carbapenem & cephalosporin centred (bold link dead-centre),
# aminoglycoside / fluoroquinolone pushed to top / bottom.
ccm_xy <- tibble::tribble(
  ~lab,               ~x,    ~y,
  "Carbapenem",       -1.5,   0,
  "Cephalosporin",     1.5,   0,
  "Aminoglycoside",    0,     1.7,
  "Fluoroquinolone",   0,    -1.7
)

# collapse the two CCM directions per pair -> one undirected edge
ccm_edges <- ccm_sum %>%
  filter(driver_lab %in% ccm_classes, affected_lab %in% ccm_classes) %>%
  mutate(a = ifelse(driver_lab < affected_lab, driver_lab, affected_lab),
         b = ifelse(driver_lab < affected_lab, affected_lab, driver_lab)) %>%
  group_by(a, b) %>%
  summarise(rho = max(rho_end), converges = any(converges), .groups = "drop") %>%
  left_join(ccm_xy, by = c("a" = "lab")) %>% rename(x1 = x, y1 = y) %>%
  left_join(ccm_xy, by = c("b" = "lab")) %>% rename(x2 = x, y2 = y)

# NOTE: this run is the UNCORRECTED one - embeddings may span patient
# boundaries and `converges` is a bare threshold with no null. Section 12.1
# rebuilds this same diagram from the surrogate-tested verdict and writes it to
# the slide path figures/ccm/fig_ccm_coupling.pdf. This version is kept for
# comparison under an explicit _uncorrected name.
fig_ccm_coupling <- ccm_coupling_diagram(
  edges   = ccm_edges,
  xy      = ccm_xy,
  title   = "Dynamical coupling among resistance classes (uncorrected)",
  caption = paste("Boundary-crossing embeddings, no surrogate null.",
                  "Bold = converging pair. Faint = tested, flat.")
)

ggsave("figures/ccm/fig_ccm_coupling_uncorrected.pdf", fig_ccm_coupling,
       width = 6.5, height = 5, dpi = 300, bg = "transparent")

# #############################################################################
# SECTION 8: DEEP-SAMPLING CASE STUDY
# External ground truth: Vanderwoude et al. 2024, two of our patients
# sampled at ~75 colonies. Illustration, not model input.
# #############################################################################

# Vanderwoude et al. 2024 (mBio; CC BY 4.0 - reusable with attribution):
#   their P1 = our patient 7  - hypermutator, heterogeneous
#   their P3 = our patient 38 - normomutator, tight
# External ground truth for within-host heterogeneity. NOT model input
# (continuous ZOI vs clinical R/I/S; depth 75 vs ~1; single timepoint; n=2).
# =============================================================================
library(readxl)
# Relative to the working directory set at the top of this file. It was an
# absolute Dropbox path, which broke when the project moved to GitHub.
supp_path <- "mbio.03109-23-s0008.xlsx"
stopifnot(file.exists(supp_path))


deep <- bind_rows(read_zoi("S1", "P1  (heterogeneous)"), #id 7
                  read_zoi("S3", "P3  (tight)")) %>%     #id 38
  filter(!is.na(zoi))

# CLSI disc-diffusion breakpoints (ZOI mm): R if <= R, S if >= S, else I
thr <- tribble(
  ~drug,                          ~R,  ~S,
  "Aminoglycoside (amikacin)",    14,  17,
  "Aminoglycoside (tobramycin)",  12,  15,
  "Carbapenem (meropenem)",       15,  19,
  "Cephalosporin (ceftazidime)",  14,  18,
  "Fluoroquinolone (cipro)",      18,  25,
  "Penicillin (pip-tazo)",        14,  21)

xmax_axis <- 52
bands <- thr %>%
  transmute(drug, r_xmin = 0, r_xmax = R, i_xmin = R, i_xmax = S,
            s_xmin = S, s_xmax = xmax_axis) %>%
  pivot_longer(-drug, names_to = c("zone", ".value"), names_sep = "_") %>%
  mutate(zone = recode(zone, r = "Resistant", i = "Intermediate", s = "Susceptible"),
         zone = factor(zone, levels = c("Resistant", "Intermediate", "Susceptible")))
band_cols <- c(Resistant = pal[["red"]], Intermediate = pal[["gold"]],
               Susceptible = pal[["blue"]])

fig_deep_zoi <- ggplot() +
  geom_rect(data = bands, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
                              fill = zone), alpha = 0.16) +
  geom_jitter(data = deep, aes(zoi, patient), height = 0.18, size = 0.7,
              alpha = 0.55, colour = pal[["ink"]]) +
  scale_fill_manual(values = band_cols, name = NULL) +
  facet_wrap(~ drug, ncol = 2) +
  labs(x = "Zone of inhibition (mm)   -   smaller = more resistant", y = NULL,
       title = "One sputum, 75 colonies: what a single colony misses",
       # attribution restored - Talk_figures_v2 dropped this caption, and the
       # Vanderwoude data is CC BY 4.0, which requires the credit to travel
       # with the figure.
       caption = "Each point = 1 of ~75 isolates from a single sample. Bands = CLSI R/I/S. Two of our patients, deep-sampled by Vanderwoude et al. 2024 (CC BY 4.0). Illustration, not model input.") +
  theme_talk() +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 18))

ggsave("figures/talk/fig_deep_zoi.pdf", fig_deep_zoi, width = 8.5, height = 7,
       dpi = 300, bg = "transparent")

# --- the bridge: deep population split vs what the clinic would call ---------
deep_pct <- deep %>% left_join(thr, by = "drug") %>%
  group_by(patient, drug) %>%
  summarise(n = n(), pctR = round(mean(zoi <= R) * 100),
            pctI = round(mean(zoi > R & zoi < S) * 100),
            pctS = round(mean(zoi >= S) * 100), .groups = "drop")
cat("\n=== Deep-sample R/I/S split (n=75) - the population a single colony is drawn from ===\n")
print(deep_pct, n = 40)
write_csv(deep_pct, "results/deep/deep_sample_RIS_split.csv")
cat("\nHeadline: id 7 amikacin = 11% R / 47% I / 43% S, but the clinic's nearest\n",
    "culture (2 colonies) called it Resistant - caught the minority tail.\n")

cat("\nSection 7 done: fig_deep_zoi.pdf + deep_sample_RIS_split.csv\n")


# =============================================================================
# 8.2  Within-population growth-rate spread (Vanderwoude S14-S17) [poster]
# Second heterogeneity axis: per-isolate max specific growth rate from the raw
# OD600 curves. The slow-grower tail is what standard culture under-samples.
# mu_max = max log-phase slope of ln(net OD); a proxy, not GrowthCurver's r.
# =============================================================================

growth <- bind_rows(
  growth_mu("S14", "P1 (hypermutator)"),
  growth_mu("S16", "P3 (normomutator)"),
  growth_mu("S15", "P2 (context)"),
  growth_mu("S17", "P4 (context)")) %>%
  filter(is.finite(mu_max)) %>%
  mutate(patient = fct_reorder(patient, mu_max, median),
         cohort  = if_else(str_detect(patient, "mutator"), "our cohort", "context"))

fig_growth_spread <- ggplot(growth, aes(mu_max, patient)) +
  geom_jitter(aes(colour = cohort), height = 0.16, size = 1.1, alpha = 0.6) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 7, colour = pal[["ink"]]) +
  scale_colour_manual(values = c("our cohort" = pal[["red"]], "context" = "grey65"),
                      name = NULL) +
  labs(x = "Max specific growth rate (/h)   -   smaller = slower grower", y = NULL,
       title = "One population, many speeds: the slow-grower tail",
       caption = "Each point = 1 of ~75 isolates (mean of triplicate OD600 curves). | = median. Log-phase slope proxy; Vanderwoude et al. 2024 (CC BY 4.0).") +
  theme_talk() +
  theme(legend.position = "top", plot.title = element_text(face = "bold", size = 15))

ggsave("figures/talk/fig_growth_spread.pdf", fig_growth_spread,
       width = 9, height = 5, dpi = 300, bg = "transparent")

growth %>% group_by(patient) %>%
  summarise(n = n(), min = round(min(mu_max), 2), median = round(median(mu_max), 2),
            max = round(max(mu_max), 2), CV = round(sd(mu_max) / mean(mu_max) * 100),
            .groups = "drop") %>% print()

# #############################################################################
# SECTION 9: SUMMARY
# Cohort panel for the talk, then the integrated comparison and key findings.
# #############################################################################


# -----------------------------------------------------------------------------
# 9.1  Cohort summary panel  [Slide 2]
# 
# -----------------------------------------------------------------------------

# Built defensively: a panel only uses columns that exist in your schema, so
# this runs as-is and you adjust column names if any differ.
# =============================================================================

cat("\n==== COHORT ====\n")
cat("  Patients:", n_distinct(pa_visit_full$individual_id), "\n")
cat("  Visits  :", nrow(pa_visit_full), "\n")
if (exists("pa_isolates")) cat("  Isolates:", nrow(pa_isolates), "\n")

# % of non-NA rows coded 1 (resistant / present), for columns that exist


# (a) OUTCOME over time: class resistance levels + their (lack of) trend ------
# Replaces BOTH the static "resistance by class" bar and the trajectory panel:
# the height of each line is the class level; the flatness is the (non-)trend.
res_cols <- c(ag_any_resistant_conservative   = "Aminoglycoside",
              carb_any_resistant_conservative = "Carbapenem",
              fq_any_resistant_conservative   = "Fluoroquinolone",
              ceph_any_resistant_conservative = "Cephalosporin")
res_time <- pa_visit_full %>%
  filter(!is.na(days_since_first_visit)) %>%
  mutate(yr = days_since_first_visit / 365.25) %>%
  pivot_longer(all_of(names(res_cols)), names_to = "class", values_to = "r") %>%
  filter(!is.na(r)) %>%
  mutate(class = recode(class, !!!res_cols), ybin = floor(yr)) %>%
  group_by(class, ybin) %>%
  summarise(rate = mean(r == 1), n = n(), .groups = "drop") %>%
  filter(n >= 15)
p_res <- ggplot(res_time, aes(ybin, rate, colour = class)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_colour_manual(values = c("Aminoglycoside"  = pal[["blue"]],
                                 "Carbapenem"      = pal[["red"]],
                                 "Fluoroquinolone" = pal[["gold"]],
                                 "Cephalosporin"   = pal[["plum"]])) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Years since first visit", y = "Resistant (R + I)",
       title = "Class resistance over time") +
  theme(legend.position = c(0.75, 0.2), legend.title = element_blank(),
        plot.title = element_text(face = "bold", size = 18))

ggsave("figures/talk/p_res.pdf", p_res, width = 6, height = 4, dpi = 300, bg = "transparent")

# (b) CONTEXT: co-colonizer prevalence (ecological confounding) ---------------
col_cols <- c(prev_has_mrsa       = "MRSA",
              prev_has_mssa       = "MSSA",
              prev_has_any_mucoid = "Mucoid PA")
col_df <- pct_present(pa_visit_full, names(col_cols)) %>%
  mutate(item = recode(item, !!!col_cols))
# subtitle computed, not typed - it was hardcoded "1,809 PA isolates"
p_col <- barpanel(col_df, "Ecological context",
                  sprintf("%s PA isolates", format(nrow(pa_isolates), big.mark = ",")),
                  "#207E9D")

ggsave("figures/talk/p_col.pdf", p_col, width = 5.5, height = 3.5, dpi = 300, bg = "transparent")

# (c) EXPOSURE noise: prescriptions by class ---------------------------------
rx_cols <- c(prev_aminoglycoside_prescribed   = "Aminoglycoside Rx",
             prev_carbapenem_prescribed       = "Carbapenem Rx",
             prev_fluoroquinolone_prescribed  = "Fluoroquinolone Rx",
             prev_cephalosporin_prescribed    = "Cephalosporin Rx",
             prev_colistin_prescribed         = "Colistin Rx")
rx_df <- pct_present(pa_visit_full, names(rx_cols)) %>%
  mutate(item = recode(item, !!!rx_cols))
p_rx <- barpanel(rx_df, "Prescribed ≠ taken/realized",
                 "", "#E8C547")
ggsave("figures/talk/p_rx.pdf", p_rx, width = 6, height = 3.5, dpi = 300, bg = "transparent")

# (d) SAMPLING noise: colonies (isolates) tested per culture -----------------
# True depth = n_isolates (independent colonies), NOT n_tested (= isolates x
# ~10 drug results). Median is 1: most cultures test a single colony.
depth_col <- intersect(c("n_isolates", "n_tested"), names(pa_visit_full))[1]
med_n <- median(pa_visit_full[[depth_col]], na.rm = TRUE)
pct1  <- round(100 * mean(pa_visit_full[[depth_col]] == 1, na.rm = TRUE))
p_depth <- ggplot(pa_visit_full, aes(.data[[depth_col]])) +
  geom_bar(fill = "#90606C") +
  coord_cartesian(xlim = c(0, 8)) +
  labs(x = "Colonies (isolates) tested per culture", y = "Visits",
       title = "Sampling: 54% test one colony",
       subtitle = sprintf("", pct1, med_n)) +
  theme(plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_blank())

ggsave("figures/talk/p_depth.pdf", p_depth, width = 5, height = 3.5, dpi = 300, bg = "transparent")

# (e) [removed] within-patient trajectories are now folded into panel (a):
# the class-resistance-over-time lines show the level, and the standalone
# faceted fig_resistance_over_time (Backup/Optional section) carries the full
# grey-trajectory "stable population, restless individuals" version.

# Computed, not typed. This strapline was hardcoded as
#   "92 PA-positive patients | 1,809 PA isolates | 1,091 visits | ~8 years"
# which silently survives every rerun: it would still have claimed 8 years and
# the pre-filter counts after DROP_PRE_2017 removed the 2012/2015 visits. Any
# number that appears on a slide has to be derived from the data that made it.
.cohort_span <- {
  d <- as.Date(pa_visit_full$date_collected_from_subject)
  as.numeric(difftime(max(d), min(d), units = "days")) / 365.25
}
cohort_caption <- sprintf(
  "%d PA-positive patients  |  %s PA isolates  |  %s visits  |  %s-%s (~%.0f years)",
  dplyr::n_distinct(pa_visit_full$individual_id),
  format(nrow(pa_isolates), big.mark = ","),
  format(nrow(pa_visit_full), big.mark = ","),
  format(min(as.Date(pa_visit_full$date_collected_from_subject)), "%Y"),
  format(max(as.Date(pa_visit_full$date_collected_from_subject)), "%Y"),
  .cohort_span)
cat("\n  cohort strapline:", cohort_caption, "\n")

fig_summary <- (p_res | p_col) / (p_rx | p_depth) +
  plot_annotation(
    caption = cohort_caption,
    theme = theme(plot.title = element_text(face = "bold", size = 16)))

ggsave("figures/talk/fig_summary_stats.pdf", fig_summary, width = 14, height = 9, dpi = 300, bg = "transparent")

# =============================================================================

# -----------------------------------------------------------------------------
# 9.2  Integrated comparison and key findings
# NOTE: the CCM columns here are guarded by exists(ccm_summary_overall) and
# will skip - see WHAT WAS DROPPED at the top of this file.
# -----------------------------------------------------------------------------


cat("\n")
cat("=" %>% strrep(70), "\n")
cat("INTEGRATED COMPARISON: GLMM vs CCM vs S-MAP\n")
cat("=" %>% strrep(70), "\n\n")

comparison_vars <- c("Carbapenem", "Aminoglycoside", "Fluoroquinolone", 
                     "Cephalosporin", "Colistin", "MRSA", "MSSA", 
                     "Mucoid")
# exacerbation omitted: constant in edm_df / pa_prescription, so it never
# enters CCM or S-map and would produce an empty row in this comparison table.

comparison_table <- tibble(Variable = comparison_vars)

# The CCM columns come from the drug-level placeholder CCM, which this pipeline
# no longer runs (see WHAT WAS DROPPED at the top). The endogenous CCM in 7.5 is
# resistance-class vs resistance-class and does not substitute: it has no row
# per drug-class PREDICTOR, so there is nothing to join on here. Guarded so the
# rest of the comparison table (S-map) still builds.
if (exists("ccm_summary_overall")) {
  comparison_table <- comparison_table %>%
    left_join(
      ccm_summary_overall %>%
        mutate(Variable = case_when(
          Predictor == "carbapenem" ~ "Carbapenem",
          Predictor == "aminoglycoside" ~ "Aminoglycoside",
          Predictor == "fluoroquinolone" ~ "Fluoroquinolone",
          Predictor == "cephalosporin" ~ "Cephalosporin",
          Predictor == "colistin" ~ "Colistin",
          Predictor == "has_mrsa" ~ "MRSA",
          Predictor == "has_mssa" ~ "MSSA",
          Predictor == "has_mucoid" ~ "Mucoid"
        )) %>%
        dplyr::select(Variable, CCM_rho = rho_end, CCM_Causal = Causal,
                      CCM_Confounded = Confounded, CCM_Interpretation = Interpretation),
      by = "Variable"
    )
} else {
  cat("\nNOTE: CCM columns omitted from the comparison table - the drug-level\n",
      "     placeholder CCM was dropped in the merge and the endogenous CCM\n",
      "     (Section 7.5) is not a drop-in replacement for it.\n", sep = "")
}

if (exists("coeffs_summary")) {
  comparison_table <- comparison_table %>%
    left_join(
      coeffs_summary %>%
        dplyr::select(Variable, Smap_Mean = Mean, Smap_CV = CV, Smap_SignChanges = Sign_Changes),
      by = "Variable"
    )
}

cat("INTEGRATED ANALYSIS COMPARISON:\n")
print(comparison_table)

write_csv(comparison_table, "results/glmm_ccm_smap_comparison.csv")

# =============================================================================
# PART 7: KEY FINDINGS SUMMARY
# =============================================================================

cat("\n")
cat("=" %>% strrep(70), "\n")
cat("KEY FINDINGS: ASYMMETRIC COLISTIN DESIGN\n")
cat("=" %>% strrep(70), "\n\n")

cat("DESIGN SUMMARY:\n")
cat("  - 4 drug classes as OUTCOMES: Carb, AG, FQ, Ceph\n")
cat("  - 5 drug classes as PREDICTORS: Carb, AG, FQ, Ceph, Colistin\n")
cat("  - Colistin excluded as outcome due to insufficient testing data\n")
cat("  - Colistin retained as predictor (45% prescription rate)\n\n")

cat("1. EDM PARAMETERS (5 targets):\n")
for (i in 1:nrow(params_summary)) {
  cat("   ", params_summary$Target[i], ": E =", params_summary$E[i], 
      ", θ =", params_summary$theta[i], " (", params_summary$Dynamics[i], ")\n")
}

cat("\n2. CCM CAUSAL EVIDENCE:\n")
if (exists("ccm_summary_overall")) {
  causal_vars <- ccm_summary_overall %>% filter(Causal == "✓") %>% pull(Predictor)
  cat("   Variables with causal evidence:", paste(causal_vars, collapse = ", "), "\n")
  
  # Colistin specifically
  coli_result <- ccm_summary_overall %>% filter(Predictor == "colistin")
  if (nrow(coli_result) > 0) {
    cat("   Colistin: rho =", round(coli_result$rho_end, 3), 
        ", Interpretation:", coli_result$Interpretation, "\n")
  }
}

cat("\n3. S-MAP STATE-DEPENDENT EFFECTS:\n")
if (exists("coeffs_summary")) {
  state_dependent <- coeffs_summary %>% filter(Sign_Changes == "Yes") %>% pull(Variable)
  if (length(state_dependent) > 0) {
    cat("   Variables with sign changes:", paste(state_dependent, collapse = ", "), "\n")
  }
  
  # Colistin specifically
  coli_smap <- coeffs_summary %>% filter(Variable == "Colistin")
  if (nrow(coli_smap) > 0) {
    cat("   Colistin: Mean =", round(coli_smap$Mean, 3), 
        ", CV =", round(coli_smap$CV, 1), "%",
        ", Sign changes:", coli_smap$Sign_Changes, "\n")
  }
}


# =============================================================================
# CAUSAL NETWORK (v2) — built dynamically from fitted models
# Pulls GLMM betas from all_effects (M2: Overall Rx) and CCM classifications
# from ccm_summary_overall, then assigns the "story" label by combining
# whether the GLMM effect is significant with what CCM says about causality.
#
# This replaces the previous hardcoded tribble so the figure stays in sync
# with the models every time the script is rerun.
# =============================================================================

library(tidyverse)

# -----------------------------------------------------------------------------
# Helper: classify a CCM result into a category used by the plot
# Rules are deliberately simple so we can defend them in the methods text.
#   - "autoregressive": not a candidate for classification (prior resistance)
#   - "bidirectional" : CCM shows convergence in BOTH directions (X→R and R→X)
#   - "reverse"       : CCM shows convergence only R→X (reverse causation)
#   - "one_way"       : CCM shows convergence only X→R
#   - "confounded"    : GLMM is significant but CCM fails to converge either
#                       way — likely a shared-cause (confounded) association
#   - "null"          : neither GLMM nor CCM signal anything
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Pull GLMM betas from M2 (the main summary model)
# -----------------------------------------------------------------------------
glmm_m2 <- all_effects %>%
  filter(Model == "M2: Overall Rx") %>%
  dplyr::select(Variable, glmm_est = Estimate, glmm_sig = Significance)

# -----------------------------------------------------------------------------
# Map M2 variable labels onto CCM predictor names, then join CCM classifications
# ccm_summary_overall uses lowercase drug names; all_effects uses "X Rx" labels
# -----------------------------------------------------------------------------
ccm_map <- tribble(
  ~glmm_label,          ~ccm_predictor,    ~plot_label,
  "Prior Resistance",   NA_character_,     "Prior Resistance",
  "Carbapenem Rx",      "carbapenem",      "Carbapenem",
  "Aminoglycoside Rx",  "aminoglycoside",  "Aminoglycoside",
  "Fluoroquinolone Rx", "fluoroquinolone", "Fluoroquinolone",
  "Cephalosporin Rx",   "cephalosporin",   "Cephalosporin",
  "Colistin Rx",        "colistin",        "Colistin"
)

# -----------------------------------------------------------------------------
# If the CCM summary exists, use it; otherwise fall back to an empty skeleton
# (the plot will still render with the GLMM betas alone).
# -----------------------------------------------------------------------------
if (exists("ccm_summary_overall")) {
  ccm_pull <- ccm_summary_overall %>%
    dplyr::select(Predictor, Causal, Confounded)
} else {
  warning("ccm_summary_overall not found — CCM classifications will be blank.")
  ccm_pull <- tibble(Predictor = character(),
                     Causal = character(),
                     Confounded = character())
}

results <- ccm_map %>%
  left_join(glmm_m2, by = c("glmm_label" = "Variable")) %>%
  left_join(ccm_pull, by = c("ccm_predictor" = "Predictor")) %>%
  rowwise() %>%
  mutate(
    ccm = if_else(
      plot_label == "Prior Resistance",
      "autoregressive",
      classify_ccm(Causal, Confounded, glmm_sig)
    ),
    story = case_when(
      ccm == "autoregressive"              ~ "Dominant predictor",
      ccm == "bidirectional"               ~ "Bidirectional — real signal",
      ccm == "reverse"                     ~ "Reverse only — Rx responds to R",
      ccm == "one_way"                     ~ "One-way causal",
      ccm == "confounded" & glmm_est > 0   ~ "GLMM positive, CCM null (confounded?)",
      ccm == "confounded" & glmm_est < 0   ~ "GLMM negative, CCM null (confounded?)",
      # No CCM was run for this variable: report the GLMM sign only, and do not
      # imply anything about causality.
      ccm == "not_assessed" & glmm_est > 0 ~ "GLMM positive — CCM not assessed",
      ccm == "not_assessed" & glmm_est < 0 ~ "GLMM negative — CCM not assessed",
      TRUE                                  ~ "No signal"
    )
  ) %>%
  ungroup() %>%
  dplyr::select(variable = plot_label, glmm_est, glmm_sig, ccm, story) %>%
  # Order: prior resistance first, then by |β| descending
  arrange(variable != "Prior Resistance", desc(abs(glmm_est))) %>%
  mutate(
    y = rev(seq_len(n())),
    color = case_when(
      ccm == "autoregressive"          ~ "#6A3D9A",   # purple — dominant
      ccm == "confounded" & glmm_est>0 ~ "#D48A8A",   # muted red — GLMM+ / CCM null
      ccm == "confounded" & glmm_est<0 ~ "#8AA8D4",   # muted blue — GLMM- / CCM null
      ccm == "bidirectional"           ~ "#FF7F00",   # orange — real feedback
      ccm == "reverse"                 ~ "#4DAF4A",   # green — reverse causation
      ccm == "one_way"                 ~ "#377EB8",   # blue — one-way causal
      ccm == "not_assessed"            ~ "grey70",   # grey — no CCM was run
      TRUE                              ~ "grey50"
    ),
    linetype = case_when(
      ccm == "confounded" ~ "dashed",
      ccm == "reverse"    ~ "dashed",
      ccm == "null"       ~ "dotted",
      TRUE                ~ "solid"
    ),
    linewidth = case_when(
      ccm == "autoregressive" ~ 2.4,
      TRUE                    ~ 1.4
    )
  )

# Print so we can eyeball the reconstruction before plotting
cat("\nCausal network inputs (built from all_effects + ccm_summary_overall):\n")
print(results %>% dplyr::select(variable, glmm_est, glmm_sig, ccm, story))

write_csv(results %>% dplyr::select(variable, glmm_est, glmm_sig, ccm, story),
          "results/causal_network_inputs.csv")

outcome <- tibble(x = 3.0, y = mean(results$y), label = "PA\nResistance")

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

p <- ggplot() +
  # Forward arrows — all variables except "reverse"
  geom_curve(
    data = filter(results, !(ccm %in% c("bidirectional", "reverse"))),
    aes(x = 0.25, xend = outcome$x - 0.35, y = y, yend = outcome$y,
        color = color, linetype = linetype, linewidth = linewidth),
    curvature = 0.12,
    arrow = arrow(type = "closed", length = unit(0.2, "inches"))
  ) +
  # Bidirectional — double-headed arrow
  geom_curve(
    data = filter(results, ccm == "bidirectional"),
    aes(x = 0.25, xend = outcome$x - 0.35, y = y, yend = outcome$y,
        color = color, linetype = linetype, linewidth = linewidth),
    curvature = 0.12,
    arrow = arrow(ends = "both", type = "closed", length = unit(0.2, "inches"))
  ) +
  # Reverse causation — arrow points FROM resistance TO drug
  geom_curve(
    data = filter(results, ccm == "reverse"),
    aes(x = outcome$x - 0.35, xend = 0.25, y = outcome$y, yend = y,
        color = color, linetype = linetype, linewidth = linewidth),
    curvature = -0.12,
    arrow = arrow(type = "closed", length = unit(0.2, "inches"))
  ) +
  # Variable nodes
  geom_point(data = results, aes(x = 0, y = y),
             size = 10, color = "grey20", fill = "white", shape = 21, stroke = 1.5) +
  # Outcome node
  geom_point(data = outcome, aes(x = x, y = y),
             size = 26, color = "grey20", fill = "#E41A1C", shape = 21, stroke = 2) +
  geom_text(data = outcome, aes(x = x, y = y, label = label),
            color = "white", fontface = "bold", size = 4.2, lineheight = 0.9) +
  # Variable labels
  geom_text(data = results, aes(x = -0.25, y = y + 0.12, label = variable),
            hjust = 1, size = 4.4, fontface = "bold") +
  geom_text(data = results,
            aes(x = -0.25, y = y - 0.22,
                label = paste0("GLMM beta = ", sprintf("%+.2f", glmm_est), " ", glmm_sig)),
            hjust = 1, size = 3.2, color = "grey40") +
  # Story annotations
  geom_text(data = results,
            aes(x = outcome$x + 1.1, y = y, label = story, color = color),
            hjust = 0, size = 4, fontface = "bold") +
  scale_color_identity() +
  scale_linetype_identity() +
  scale_linewidth_identity() +
  coord_cartesian(
    xlim = c(-2.8, 8.2),
    ylim = c(0.3, max(results$y) + 0.6),
    clip = "off"
  ) +
  labs(
    title = "What's actually driving resistance?",
    # Data-driven: only claim CCM when a CCM actually contributed a verdict.
    # With the placeholder CCM dropped, every row is "not_assessed" and the
    # figure must not advertise causal convergence it never computed.
    subtitle = if (all(results$ccm %in% c("not_assessed", "autoregressive"))) {
      "GLMM significance only — CCM not assessed for these predictors"
    } else {
      "GLMM significance × CCM causal convergence × Markov transitions"
    }
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 17, hjust = 0,
                                 margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11.5, color = "grey30", hjust = 0,
                                 margin = margin(b = 18)),
    plot.margin   = margin(22, 22, 22, 22)
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/causal_network_v2.pdf", p, width = 13, height = 7)

print(p)


# =============================================================================
# =============================================================================


# #############################################################################
# SECTION 10: DESIGN SENSITIVITY  -  colony denominators and an estimable
#             between-colony variance
#
# Added August 2026 in response to statistical review. Nothing above this line
# is modified: Sections 1-9 remain the ORIGINAL design, and every comparison
# here refits the original specification on the SAME rows as the alternative so
# the two are directly comparable. Old numbers stay quotable; new numbers sit
# beside them.
#
# The review's three load-bearing points, and what is tested here:
#   10.2  The beta-binomial denominator is drug-tests, not colonies. At the
#         isolate level ag_n_tested is 0-3 (three AG drugs); summed to the visit
#         it is ~3 x n_isolates. phi therefore mixes cross-drug correlation
#         WITHIN one colony with genuine variation BETWEEN colonies, and the
#         within-colony component is the stronger of the two. Refit with
#         trials = colonies tested and successes = colonies non-susceptible.
#   10.3  pi^2/3 is not an estimated variance component. It is the fixed
#         logistic link variance, so the "within-visit share" reported in
#         Section 3 is arithmetic, not evidence. Fit at the drug x isolate
#         level instead, where an isolate-level random effect IS identified,
#         and read the between-colony SD off the posterior.
#   10.4  Sampling depth inflates any P(at least one colony non-susceptible)
#         outcome mechanically. Refit the composite model with n_isolates as a
#         covariate and see whether the drug coefficients move.
#
# NOTE on 10.2 sample size: the five-predictor het_data set leaves only 155
# visits once ag_n_iso >= 2 is imposed, which is too thin to compare dispersion
# parameters. These models use the reduced predictor set (prior AG resistance,
# prior MRSA, prior MSSA), which holds 394 visits / 61 patients. prev mucoid and
# prev AG prescription are the two NA-heavy predictors - the same tension noted
# at het_data in Section 4.
# #############################################################################

# NB: the helper date column is called visit_date, not .d. A dot-prefixed name
# like .d partial-matches dplyr's own .data argument inside mutate()/group_by(),
# which silently pushes the data frame into ... and makes every column
# "not found". Keep helper columns un-dotted.
dir.create("models/sensitivity",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/sensitivity", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/sensitivity", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 10.1  Build colony-denominator variables
# For each class: how many COLONIES were tested for that class at the visit, and
# how many of those colonies were non-susceptible (R or I) to any agent in it.
# -----------------------------------------------------------------------------

class_drugs <- list(
  ag   = c("amikacin", "tobramycin", "gentamicin"),
  carb = c("meropenem", "imipenem"),
  ceph = c("cefepime", "ceftazidime"),
  fq   = c("ciprofloxacin", "levofloxacin"),
  pen  = "piperacillin_tazobactam",
  mono = "aztreonam",
  poly = "colistin"
)

iso_flags <- pa_isolates %>%
  mutate(visit_date = as.Date(date_collected_from_subject))

for (cl in names(class_drugs)) {
  dr <- class_drugs[[cl]]
  iso_flags[[paste0(cl, "_iso_tested")]] <- as.integer(
    rowSums(sapply(dr, function(x)
      iso_flags[[x]] %in% c("Susceptible", "Intermediate", "Resistant"))) > 0)
  iso_flags[[paste0(cl, "_iso_nonsusc")]] <- as.integer(
    rowSums(sapply(dr, function(x)
      iso_flags[[x]] %in% c("Resistant", "Intermediate"))) > 0)
}

colony_den <- iso_flags %>%
  group_by(individual_id, visit_date) %>%
  summarise(across(ends_with("_iso_tested"),  sum, .names = "{.col}"),
            across(ends_with("_iso_nonsusc"), sum, .names = "{.col}"),
            n_colonies = n(), .groups = "drop") %>%
  rename_with(~ str_replace(.x, "_iso_tested$",  "_n_iso"), ends_with("_iso_tested")) %>%
  rename_with(~ str_replace(.x, "_iso_nonsusc$", "_k"),     ends_with("_iso_nonsusc"))

pa_visit_sens <- pa_visit_full %>%
  mutate(visit_date = as.Date(date_collected_from_subject)) %>%
  left_join(colony_den, by = c("individual_id", "visit_date"))

write_csv(colony_den, "results/sensitivity/colony_denominators.csv")

cat("\n=== 10.1 colony denominators built ===\n")
cat("  visits:", nrow(colony_den), "| colonies:", sum(colony_den$n_colonies), "\n")
cat("  AG: visits with >=2 colonies tested:", sum(colony_den$ag_n_iso >= 2), "\n")
cat("  Note ag_n_tested (drug-tests) vs ag_n_iso (colonies):\n")
print(pa_visit_sens %>%
        summarise(mean_drugtests = mean(ag_n_tested), mean_colonies = mean(ag_n_iso, na.rm = TRUE),
                  ratio = mean(ag_n_tested) / mean(ag_n_iso, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 10.2  Same rows, two denominators: drug-tests vs colonies
# -----------------------------------------------------------------------------

sens_df <- pa_visit_sens %>%
  filter(ag_n_iso >= 2,
         !is.na(prev_ag_any_resistant), !is.na(prev_has_mrsa), !is.na(prev_has_mssa))

cat("\n=== 10.2 comparison frame ===\n")
cat("  N =", nrow(sens_df), "visits,", n_distinct(sens_df$individual_id), "patients\n")

priors_sens <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 2),   class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(gamma(2, 0.1),  class = "phi")
)

# (a) ORIGINAL design: denominator = drug-tests
m_den_drugtests <- brm(
  ag_n_resistant | trials(ag_n_tested) ~
    prev_ag_any_resistant + prev_has_mrsa + prev_has_mssa + (1 | individual_id),
  data = sens_df, family = beta_binomial(), prior = priors_sens,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/sensitivity/m_den_drugtests", file_refit = "on_change"
)

# (b) ALTERNATIVE: denominator = colonies
m_den_colonies <- brm(
  ag_k | trials(ag_n_iso) ~
    prev_ag_any_resistant + prev_has_mrsa + prev_has_mssa + (1 | individual_id),
  data = sens_df, family = beta_binomial(), prior = priors_sens,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/sensitivity/m_den_colonies", file_refit = "on_change"
)

den_compare <- bind_rows(
  as_draws_df(m_den_drugtests) %>% as_tibble() %>%
    dplyr::select(phi, sd_patient = `sd_individual_id__Intercept`,
                  b_prior = b_prev_ag_any_resistant) %>%
    mutate(denominator = "Drug-tests (original)"),
  as_draws_df(m_den_colonies) %>% as_tibble() %>%
    dplyr::select(phi, sd_patient = `sd_individual_id__Intercept`,
                  b_prior = b_prev_ag_any_resistant) %>%
    mutate(denominator = "Colonies (alternative)")
) %>%
  pivot_longer(c(phi, sd_patient, b_prior), names_to = "parameter", values_to = "value") %>%
  group_by(denominator, parameter) %>%
  summarise(median = median(value),
            lo = quantile(value, 0.025), hi = quantile(value, 0.975), .groups = "drop")

cat("\n=== 10.2 RESULT: does the denominator change the dispersion? ===\n")
print(den_compare %>% arrange(parameter, denominator))
write_csv(den_compare, "results/sensitivity/denominator_comparison.csv")

fig_den <- den_compare %>%
  mutate(parameter = recode(parameter,
                            phi = "phi (dispersion)",
                            sd_patient = "SD between patients",
                            b_prior = "beta: prior AG resistance")) %>%
  ggplot(aes(x = median, y = denominator, colour = denominator)) +
  geom_pointrange(aes(xmin = lo, xmax = hi), linewidth = 1) +
  facet_wrap(~parameter, scales = "free_x") +
  scale_colour_manual(values = c("Drug-tests (original)"  = unname(pal["sand"]),
                                 "Colonies (alternative)" = unname(pal["red"])),
                      guide = "none") +
  labs(title = "Denominator sensitivity: drug-tests vs colonies",
       subtitle = sprintf("Same %d visits, %d patients; only the trials() denominator differs",
                          nrow(sens_df), n_distinct(sens_df$individual_id)),
       x = "Posterior median (95% CrI)", y = NULL)

ggsave("figures/sensitivity/fig_denominator_comparison.pdf", fig_den,
       width = 11, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# 10.3  An ESTIMATED between-colony variance, replacing the pi^2/3 convention
#
# Fit at the drug x isolate level. AG is the only class with three agents, so it
# is the only place an isolate-level random effect is well identified: it is
# pinned by the 444 visits that contributed more than one colony.
#
# drug enters as a FIXED effect, not (1 | drug). Three levels is too few to
# estimate a variance component from, and the quantity of interest here is
# sd_isolate_id, not the spread across drugs.
# -----------------------------------------------------------------------------

ag_long <- pa_isolates %>%
  mutate(visit_date = as.Date(date_collected_from_subject),
         isolate_id = paste0("iso", row_number()),
         visit_key  = paste(individual_id, visit_date, sep = "_")) %>%
  dplyr::select(individual_id, visit_key, isolate_id, amikacin, tobramycin, gentamicin) %>%
  pivot_longer(c(amikacin, tobramycin, gentamicin), names_to = "drug", values_to = "res") %>%
  filter(res %in% c("Susceptible", "Intermediate", "Resistant")) %>%
  mutate(nonsusc = as.integer(res %in% c("Resistant", "Intermediate")))

cat("\n=== 10.3 drug x isolate frame ===\n")
cat("  rows:", nrow(ag_long), "| isolates:", n_distinct(ag_long$isolate_id),
    "| visits:", n_distinct(ag_long$visit_key),
    "| patients:", n_distinct(ag_long$individual_id), "\n")

# --- what actually identifies sd_isolate_id ---------------------------------
# An isolate with a single AG result contributes nothing to separating the
# isolate random effect from residual. Report the subset that does the work.
per_iso <- ag_long %>% count(isolate_id, name = "n_drugs")
cat("  isolates with >=2 AG results (these identify the isolate RE):",
    sum(per_iso$n_drugs >= 2), sprintf("(%.0f%%)\n", 100 * mean(per_iso$n_drugs >= 2)))
print(per_iso %>% count(n_drugs, name = "n_isolates"))

# --- reflex testing check ----------------------------------------------------
# If a second agent is only run when the first comes back resistant, then WHICH
# drugs have results is informative about the isolate's resistance, and
# sd_isolate_id absorbs that missingness rather than biology. Test directly:
# is drug B more likely to be tested when drug A is non-susceptible?
ag_wide <- pa_isolates %>%
  transmute(isolate_id = paste0("iso", row_number()),
            amikacin, tobramycin, gentamicin)
reflex <- purrr::map_dfr(
  list(c("tobramycin","amikacin"), c("tobramycin","gentamicin"),
       c("amikacin","tobramycin"), c("amikacin","gentamicin"),
       c("gentamicin","amikacin"), c("gentamicin","tobramycin")),
  function(pr) {
    a <- pr[1]; b <- pr[2]
    d <- ag_wide %>%
      filter(.data[[a]] %in% c("Susceptible","Intermediate","Resistant")) %>%
      mutate(a_nonsusc = .data[[a]] %in% c("Resistant","Intermediate"),
             b_tested  = .data[[b]] %in% c("Susceptible","Intermediate","Resistant"))
    if (n_distinct(d$a_nonsusc) < 2) return(NULL)
    tb <- table(d$a_nonsusc, d$b_tested)
    pv <- suppressWarnings(tryCatch(fisher.test(tb)$p.value, error = function(e) NA_real_))
    tibble(given = a, tested = b,
           pct_tested_if_S  = 100 * mean(d$b_tested[!d$a_nonsusc]),
           pct_tested_if_NS = 100 * mean(d$b_tested[d$a_nonsusc]),
           p = pv)
  })
cat("\n  Reflex-testing check - is testing of B conditional on A's result?\n")
print(reflex %>% mutate(across(starts_with("pct"), ~round(., 1)), p = signif(p, 3)))
cat("  Large gaps between the two percentage columns => informative missingness,\n")
cat("  and sd_isolate_id below is inflated by it.\n")
write_csv(reflex, "results/sensitivity/reflex_testing_check.csv")

m_var_isolate <- brm(
  nonsusc ~ drug + (1 | individual_id) + (1 | visit_key) + (1 | isolate_id),
  data = ag_long, family = bernoulli(),
  prior = c(prior(normal(0, 1.5), class = "b"),
            prior(normal(0, 2),   class = "Intercept"),
            prior(exponential(1), class = "sd")),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/sensitivity/m_var_isolate", file_refit = "on_change"
)

vc <- as_draws_df(m_var_isolate) %>% as_tibble()
var_iso <- tibble(
  component = c("Between patients", "Between visits (within patient)",
                "Between colonies (within visit)", "Residual (fixed by link)"),
  estimated = c(TRUE, TRUE, TRUE, FALSE),
  variance  = c(median(vc$`sd_individual_id__Intercept`)^2,
                median(vc$`sd_visit_key__Intercept`)^2,
                median(vc$`sd_isolate_id__Intercept`)^2,
                pi^2 / 3)
) %>%
  mutate(share = variance / sum(variance))

cat("\n=== 10.3 RESULT: between-colony variance, now ESTIMATED ===\n")
cat("  sd_isolate_id (between colonies, logit scale): ",
    sprintf("%.3f [%.3f, %.3f]\n",
            median(vc$`sd_isolate_id__Intercept`),
            quantile(vc$`sd_isolate_id__Intercept`, 0.025),
            quantile(vc$`sd_isolate_id__Intercept`, 0.975)))
print(var_iso %>% mutate(across(c(variance, share), ~round(., 3))))
cat("\n  Section 3 reports the within-visit slice as a share of a total whose\n")
cat("  residual term is fixed at pi^2/3 = 3.29. Compare that against the\n")
cat("  'Between colonies' row above, which is estimated from the data.\n")
cat("\n  CAVEAT, and it is not fixable with these data: sd_isolate_id is true\n")
cat("  between-colony BIOLOGICAL variation PLUS per-colony AST measurement\n")
cat("  error. With one AST per colony the two are not separable, so this term\n")
cat("  is an UPPER bound on between-colony biology. Replicate ZOI measurements\n")
cat("  (Vanderwoude, Section 8) could bound the AST-reproducibility component\n")
cat("  and turn this caveat into a result. Do not present it as pure biology.\n")
write_csv(var_iso, "results/sensitivity/variance_isolate_level.csv")

fig_var_iso <- var_iso %>%
  mutate(component = fct_reorder(component, share),
         lab = ifelse(estimated, "estimated", "fixed by logistic link")) %>%
  ggplot(aes(x = share, y = component, fill = lab)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::percent(share, accuracy = 0.1)), hjust = -0.1, size = 4) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.18))) +
  scale_fill_manual(values = c("estimated" = unname(pal["blue"]),
                               "fixed by logistic link" = "grey70"), name = NULL) +
  labs(title = "Variance decomposition with an estimated between-colony term",
       subtitle = "Drug x isolate level, aminoglycosides. The grey bar is a constant, not a finding.",
       x = "Share of total latent variance", y = NULL) +
  theme(legend.position = "bottom")

ggsave("figures/sensitivity/fig_variance_isolate_level.pdf", fig_var_iso,
       width = 10, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 10.3b  Is the near-zero between-visit share an artefact, or is it real?
#
# The between-visit share above comes out near zero, which reads as "patients
# have no trajectory". Before accepting that, it has to survive an
# identifiability check: 58% of visits contribute a SINGLE colony, and at a
# single-colony visit the visit random effect and the isolate random effect are
# perfectly confounded - nothing distinguishes "this visit was different" from
# "this colony was different". The isolate term is then free to absorb
# visit-level variance across most of the data. (The original Section 3
# decomposition put 3.5% on visits, seven times higher, precisely because it had
# no isolate term to absorb it.)
#
# So refit on multi-colony visits only, where the visit term IS identifiable.
#
# RESULT: the between-visit share does NOT reappear. It is ~0.3% in both, and
# sd_visit_key is 0.277 (all) vs 0.271 (multi-colony) - indistinguishable.
#
# READ THAT CAREFULLY - IT DOES NOT MEAN WHAT IT FIRST APPEARS TO.
# This block was originally written to conclude "the near-zero visit term is
# therefore real, not an artefact". That conclusion is WITHDRAWN. The test asks
# "does the estimate CHANGE when the term becomes identifiable?" when the
# question it needs to answer is "can the estimator SEE anything at all?" - and
# a test with no power returns "unchanged" whatever the truth is. Section 10.3d
# runs the test that does have an answer: simulating at this cohort's real
# colony-count distribution, a TRUE sd_visit of ZERO comes back as 0.267 and a
# true 0.8 comes back as 0.67, against the 0.277 observed here. Every interval
# spans the whole range. The component is not identified, so agreement between
# two subsets is agreement between two uninformative estimates.
#
# The multi-colony subset does not rescue this either: most of those visits have
# exactly 2 colonies, so the identifiability gain over the full data is small,
# which is itself part of why nothing moves.
#
# What this block still legitimately shows: the estimate is STABLE across
# subsets. That is worth knowing - it rules out the specific story where the
# isolate term visibly swallows visit variance at single-colony visits and
# disgorges it when restricted. It is not evidence that patients hold a constant
# level. For the better-identified estimate of the same quantity, see the
# within-patient autocorrelation in Section 14.2: 0.34 at short gaps, decaying
# over ~86 days, which measurement noise cannot produce.
#
# This also reinterprets a model-free check run on the same data. Comparing
# observed within-patient visit-to-visit variance against binomial sampling
# gives an excess of ~18% among multi-colony visits, which looks like evidence
# of visit-level movement. It is not: that test cannot separate "this visit's
# mean moved" from "the colonies at this visit disagree with each other". The
# mixed model can, and it assigns the excess to COLONIES - between-colony share
# rises from 19% to 25% in the multi-colony subset while the visit share stays
# flat. Within-host heterogeneity at a point in time, not change over time, is
# what single-colony sampling misses.
#
# CAVEAT on the subset: multi-colony visits are not randomly chosen. If a
# technician picks more colonies when a plate looks mixed, the subset is
# enriched for heterogeneous visits, which would inflate the between-colony
# share here. That works against the between-colony reading, not the
# between-visit one, so it does not rescue the trajectory interpretation.
# -----------------------------------------------------------------------------

multi_visits <- ag_long %>%
  distinct(visit_key, isolate_id) %>%
  count(visit_key, name = "n_colonies") %>%
  filter(n_colonies >= 2)

ag_long_multi <- ag_long %>% filter(visit_key %in% multi_visits$visit_key)

cat("\n=== 10.3b decomposition on multi-colony visits only ===\n")
cat("  visits:", n_distinct(ag_long_multi$visit_key), "of", n_distinct(ag_long$visit_key),
    sprintf("(%.0f%%)", 100 * n_distinct(ag_long_multi$visit_key) / n_distinct(ag_long$visit_key)), "\n")
cat("  rows:", nrow(ag_long_multi), "| isolates:", n_distinct(ag_long_multi$isolate_id),
    "| patients:", n_distinct(ag_long_multi$individual_id), "\n")

m_var_isolate_multi <- brm(
  nonsusc ~ drug + (1 | individual_id) + (1 | visit_key) + (1 | isolate_id),
  data = ag_long_multi, family = bernoulli(),
  prior = c(prior(normal(0, 1.5), class = "b"),
            prior(normal(0, 2),   class = "Intercept"),
            prior(exponential(1), class = "sd")),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/sensitivity/m_var_isolate_multi", file_refit = "on_change"
)

vc_m <- as_draws_df(m_var_isolate_multi) %>% as_tibble()
var_iso_multi <- tibble(
  component = c("Between patients", "Between visits (within patient)",
                "Between colonies (within visit)", "Residual (fixed by link)"),
  estimated = c(TRUE, TRUE, TRUE, FALSE),
  variance  = c(median(vc_m$`sd_individual_id__Intercept`)^2,
                median(vc_m$`sd_visit_key__Intercept`)^2,
                median(vc_m$`sd_isolate_id__Intercept`)^2,
                pi^2 / 3)
) %>% mutate(share = variance / sum(variance))

var_compare <- bind_rows(
  var_iso       %>% mutate(subset = sprintf("All visits (n = %d)", n_distinct(ag_long$visit_key))),
  var_iso_multi %>% mutate(subset = sprintf(">= 2 colonies (n = %d)", n_distinct(ag_long_multi$visit_key)))
)

cat("\n=== 10.3b RESULT: does the visit term reappear when it can be measured? ===\n")
print(var_compare %>%
        dplyr::select(subset, component, share) %>%
        pivot_wider(names_from = subset, values_from = share) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

vis_all   <- var_iso$share[var_iso$component == "Between visits (within patient)"]
vis_multi <- var_iso_multi$share[var_iso_multi$component == "Between visits (within patient)"]
cat(sprintf("\n  between-visit share: %.1f%% (all visits) -> %.1f%% (multi-colony only)\n",
            100 * vis_all, 100 * vis_multi))
cat("  sd_visit_key: ",
    sprintf("%.3f [%.3f, %.3f] all", median(vc$`sd_visit_key__Intercept`),
            quantile(vc$`sd_visit_key__Intercept`, .025), quantile(vc$`sd_visit_key__Intercept`, .975)),
    sprintf(" | %.3f [%.3f, %.3f] multi\n", median(vc_m$`sd_visit_key__Intercept`),
            quantile(vc_m$`sd_visit_key__Intercept`, .025), quantile(vc_m$`sd_visit_key__Intercept`, .975)))
if (abs(vis_multi - vis_all) < 0.02) {
  cat("  The visit term does not grow when restricted to multi-colony visits.\n")
  cat("  That rules out one specific artefact - the isolate term visibly\n")
  cat("  absorbing visit variance and giving it back on restriction - and it\n")
  cat("  does NOT establish that the near-zero share is real. See 10.3d: at\n")
  cat("  this depth a TRUE sd_visit of zero and a true 0.8 both return values\n")
  cat("  compatible with the 0.277 observed, so neither subset is informative\n")
  cat("  about the between-visit component. Two uninformative estimates\n")
  cat("  agreeing is not evidence. Do not quote this as stability.\n")
} else {
  cat("  The visit term CHANGED once identifiable - the full-data estimate was\n")
  cat("  absorbing visit variance into the isolate term. Quote the multi-colony\n")
  cat("  figure, and note that subset is not randomly chosen.\n")
}
write_csv(var_compare, "results/sensitivity/variance_isolate_level_contrast.csv")

fig_var_contrast <- var_compare %>%
  mutate(component = factor(component, levels = rev(c(
           "Between patients", "Between visits (within patient)",
           "Between colonies (within visit)", "Residual (fixed by link)"))),
         lab = ifelse(estimated, "estimated", "fixed by logistic link")) %>%
  ggplot(aes(x = share, y = component, fill = lab)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = scales::percent(share, accuracy = 0.1)), hjust = -0.1, size = 3.6) +
  facet_wrap(~subset) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.22))) +
  scale_fill_manual(values = c("estimated" = unname(pal["blue"]),
                               "fixed by logistic link" = "grey70"), name = NULL) +
  labs(title = "The between-visit term is unchanged by the multi-colony restriction",
       subtitle = paste("Restricting to visits with >= 2 colonies leaves the visit share unchanged - which rules out one",
                        "specific artefact,\nbut does NOT show the near-zero share is real. Section 10.3d shows a true",
                        "sd_visit of 0 and of 0.8 both\nreturn values compatible with what is observed: the component is",
                        "not identified at this sampling depth."),
       x = "Share of total latent variance", y = NULL) +
  theme(legend.position = "bottom")

ggsave("figures/sensitivity/fig_variance_isolate_contrast.pdf", fig_var_contrast,
       width = 12, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 10.3c  The same finding, shown rather than tabulated - and generalised
#
# (a) Does the near-zero between-visit term hold beyond aminoglycosides?
#     Moment decomposition of visit-level proportions, on visits with m >= 2
#     where signal and sampling noise are separable.
# (b) The illustration: patient trajectories look erratic, but most of the
#     movement is inside the sampling interval of a flat line.
# (c) What depth would be needed before a real change becomes visible.
# -----------------------------------------------------------------------------

# ---- (a) all four classes ---------------------------------------------------
btw_visit <- purrr::map_dfr(names(design_classes), function(cl) {
  pcol <- paste0(cl, "_prop_resistant_conservative"); mcol <- paste0(cl, "_n_tested")
  d <- pa_visit_full %>%
    filter(!is.na(.data[[pcol]]), !is.na(.data[[mcol]]), .data[[mcol]] >= 2) %>%
    transmute(id = individual_id, p = .data[[pcol]], m = .data[[mcol]]) %>%
    group_by(id) %>% filter(n() >= 3) %>%
    mutate(pbar = mean(p), dev = p - pbar, ni = n()) %>% ungroup()
  if (nrow(d) < 50) return(NULL)
  V_between   <- var(d$pbar)
  # mean(dev^2) is a BIASED estimator of within-patient variance: dev is taken
  # around the patient's OWN mean, estimated from the same visits, so its
  # expectation is V_within * (n_i - 1)/n_i. At the median 11 visits per patient
  # for carbapenem and FQ that is a 13% downward bias - and the leftover in
  # those two classes was only -0.003, small enough that the bias alone
  # manufactured it. Pool the per-patient sums of squares over sum(n_i - 1).
  V_within_ob <- {
    ss <- d %>% group_by(id) %>% summarise(ss = sum(dev^2), ni = first(ni), .groups = "drop")
    sum(ss$ss) / sum(ss$ni - 1)
  }
  V_sampling  <- mean(d$pbar * (1 - d$pbar) / d$m)
  leftover    <- V_within_ob - V_sampling          # NOT clamped - see below
  V_within    <- max(leftover, 0)
  tibble(class = unname(design_classes[cl]), n_visits = nrow(d), n_patients = n_distinct(d$id),
         median_visits = median(d$ni),
         between_patient = V_between / (V_between + V_within),
         between_visit   = V_within  / (V_between + V_within),
         within_obs = V_within_ob, sampling_expected = V_sampling,
         # reported raw so a boundary hit is visible as a boundary hit rather
         # than as a substantive "zero between-visit variance"
         leftover_raw = leftover, at_boundary = leftover <= 0)
})
cat("\n=== 10.3c(a) between-patient vs between-visit, all classes ===\n")
print(btw_visit %>% mutate(across(where(is.numeric), ~round(., 3))))
if (any(btw_visit$at_boundary)) {
  cat("  NOTE: ", sum(btw_visit$at_boundary), " class(es) sit at the zero boundary.\n", sep = "")
  cat("  A share of exactly 0.000 is a BOUNDARY FLAG, not a finding - the\n")
  cat("  estimator was clamped, so it reports the smallest value it can rather\n")
  cat("  than the value it measured. Read leftover_raw, not between_visit.\n")
} else {
  cat("  All four classes carry POSITIVE between-visit variance once the\n")
  cat("  within-patient estimator is degrees-of-freedom corrected. The earlier\n")
  cat("  version of this block used mean(dev^2) and reported exactly 0.000 for\n")
  cat("  carbapenem and FQ, concluding 'there is nothing left to be a\n")
  cat("  trajectory'. That conclusion was an artefact of the biased estimator,\n")
  cat("  and it is retracted here: correcting the bias flips both classes\n")
  cat("  positive (about 0.04 each). This is the same failure mode as the 0.3%\n")
  cat("  in 10.3d - an estimate pinned at a boundary, read as a result.\n")
}
cat("  Apparent UNDER-dispersion (leftover < 0) is also diagnostically wrong on\n")
cat("  its own terms: within-colony drug correlation of 0.45-0.70 produces\n")
cat("  OVER-dispersion relative to an independent-binomial expectation, never\n")
cat("  under. A negative leftover indicts the estimator, not the biology.\n")
write_csv(btw_visit, "results/sensitivity/between_visit_by_class.csv")

# ---- (b) the illustration ---------------------------------------------------
traj <- pa_visit_full %>%
  filter(!is.na(ag_prop_resistant_conservative), ag_n_tested >= 1) %>%
  transmute(id = individual_id, t = days_since_first_visit / 365.25,
            p = ag_prop_resistant_conservative, n = ag_n_tested, k = ag_n_resistant) %>%
  group_by(id) %>% filter(n() >= 6) %>%
  mutate(pbar_pool = sum(k) / sum(n)) %>% ungroup() %>%
  mutate(se_v   = sqrt(pbar_pool * (1 - pbar_pool) / n),
         lo     = pmax(p - 1.96 * se_v, 0), hi = pmin(p + 1.96 * se_v, 1),
         inside = abs(p - pbar_pool) <= 1.96 * se_v)

pct_flat <- round(100 * mean(traj$inside, na.rm = TRUE))
cat(sprintf("\n=== 10.3c(b) %d%% of visits sit within sampling error of the patient's own flat level ===\n",
            pct_flat))
cat("  patients with >= 6 visits:", n_distinct(traj$id), "| visits:", nrow(traj), "\n")

set.seed(3)
traj_pick <- traj %>% distinct(id, pbar_pool) %>%
  filter(pbar_pool > 0.15, pbar_pool < 0.85) %>% slice_sample(n = 6) %>% pull(id)
traj_sub  <- traj %>% filter(id %in% traj_pick) %>% mutate(pat = paste("Patient", dense_rank(id)))
traj_line <- traj_sub %>% distinct(pat, pbar_pool)

fig_traj <- ggplot(traj_sub, aes(t, p)) +
  geom_hline(data = traj_line, aes(yintercept = pbar_pool), colour = pal[["blue"]], linewidth = 1) +
  geom_line(colour = "grey55", linewidth = 0.4, alpha = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = "grey45", alpha = 0.8) +
  geom_point(aes(colour = inside), size = 2.2) +
  facet_wrap(~pat, ncol = 3) +
  scale_colour_manual(values = c(`TRUE` = pal[["red"]], `FALSE` = "#1A1A1A"),
                      labels = c(`TRUE` = "consistent with a flat level",
                                 `FALSE` = "genuine departure"), name = NULL) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.02, 1.02)) +
  labs(x = "Years since first visit", y = "AG resistance (R + I)",
       title = "Trajectories look erratic. At this depth, most movement is indistinguishable from noise.",
       subtitle = sprintf(paste("Blue = the patient's own pooled level. Bars = 95%% sampling interval at that visit's depth.",
                                "\nAcross all patients, %d%% of visits sit within sampling error of a flat line."),
                          pct_flat)) +
  theme(legend.position = "bottom")

ggsave("figures/sensitivity/fig_trajectories_are_noise.pdf", fig_traj,
       width = 12, height = 7, dpi = 300)

# ---- (c) what change could you even see, at a given depth? ------------------
# Smallest shift in resistance frequency distinguishable from sampling noise,
# taken as 2 x 1.96 x SE at p = 0.5 (the worst case).
detectable <- tibble(colonies = c(1, 2, 3, 5, 10, 20, 50, 75)) %>%
  mutate(se = sqrt(0.25 / colonies),
         min_detectable_shift = pmin(2 * 1.96 * se, 1))
cat("\n=== 10.3c(c) smallest visit-to-visit shift distinguishable from noise ===\n")
print(detectable %>% mutate(across(where(is.numeric), ~round(., 3))))
cat("  At one colony nothing smaller than the full 0-1 range is detectable,\n")
cat("  which is why 93% of visits are 'consistent with flat'. The question is\n")
cat("  not whether patients move - it is whether the measurement could see it.\n")
write_csv(detectable, "results/sensitivity/detectable_shift_by_depth.csv")

# -----------------------------------------------------------------------------
# 10.3d  Can the between-visit component be recovered AT ALL at this depth?
#
# 10.3b restricted to multi-colony visits and found the visit share unchanged,
# which I read as "the near-zero term is real". That was the wrong test. The
# right one asks whether the ESTIMATOR can recover a KNOWN sd_visit given this
# cohort's actual colony-count distribution (median 1 colony per visit, where
# the visit and isolate effects are perfectly confounded).
#
# Simulate on the real design skeleton with sd_visit set to 0, 0.3, 0.5, 0.8,
# refit, and see what comes back.
#
# RESULT (2026-08, AG):
#   true 0.0 -> est 0.267 [0.014, 0.801]
#   true 0.3 -> est 0.215 [0.011, 0.695]
#   true 0.5 -> est 0.401 [0.020, 0.960]
#   true 0.8 -> est 0.666 [0.053, 1.087]
# The observed value on real data is 0.277.
#
# Two things follow, and they point in opposite directions from the two
# readings that have been argued:
#   - The component does NOT collapse to zero under the exponential(1) prior:
#     a true 0.8 comes back as 0.67, so "the prior crushed it" is not what
#     happened.
#   - But every interval spans essentially the whole range tested. A true value
#     of ZERO returns 0.267, indistinguishable from the 0.277 observed. The
#     between-visit variance component is simply NOT IDENTIFIED at this depth.
#
# So the 0.3% share must not be quoted in either direction - not as evidence
# that patients are stable, and not as an artefact to be corrected. It is
# uninformative. The within-patient autocorrelation in Section 14.2 is the
# better-identified estimate of the same quantity, and it is NOT near zero
# (0.34 at short gaps, decaying over ~86 days). Where the two disagree, prefer
# the autocorrelation: it does not depend on separating two confounded random
# effects, and a clean exponential decay with interval is a signature that
# interval-independent measurement noise cannot produce.
#
# Set RUN_SDVISIT_SIM to TRUE to regenerate (about 4 model fits).
# -----------------------------------------------------------------------------

RUN_SDVISIT_SIM <- FALSE

if (RUN_SDVISIT_SIM) {
  sd_skel <- ag_long %>% dplyr::select(individual_id, visit_key, isolate_id, drug)
  sim_sd_visit <- function(sd_v, sd_pat = 3.75, sd_iso = 2.03) {
    set.seed(1)
    u_p <- distinct(sd_skel, individual_id) %>% mutate(u = rnorm(n(), 0, sd_pat))
    u_v <- distinct(sd_skel, visit_key)     %>% mutate(u = rnorm(n(), 0, sd_v))
    u_i <- distinct(sd_skel, isolate_id)    %>% mutate(u = rnorm(n(), 0, sd_iso))
    d <- sd_skel %>%
      left_join(u_p, by = "individual_id") %>% rename(up = u) %>%
      left_join(u_v, by = "visit_key")     %>% rename(uv = u) %>%
      left_join(u_i, by = "isolate_id")    %>% rename(ui = u) %>%
      mutate(nonsusc = rbinom(n(), 1, plogis(up + uv + ui)))
    fit <- brm(nonsusc ~ drug + (1|individual_id) + (1|visit_key) + (1|isolate_id),
               data = d, family = bernoulli(),
               prior = c(prior(normal(0,1.5), class="b"), prior(normal(0,2), class="Intercept"),
                         prior(exponential(1), class="sd")),
               chains = 2, cores = 2, iter = 2000, warmup = 1000,
               control = list(adapt_delta = 0.95), seed = 42, refresh = 0, silent = 2)
    dr <- as_draws_df(fit) %>% as_tibble()
    tibble(sd_visit_true = sd_v,
           sd_visit_est = median(dr$`sd_visit_key__Intercept`),
           lo = quantile(dr$`sd_visit_key__Intercept`, .025),
           hi = quantile(dr$`sd_visit_key__Intercept`, .975))
  }
  sd_recovery <- purrr::map_dfr(c(0, 0.3, 0.5, 0.8), sim_sd_visit)
  cat("\n=== 10.3d sd_visit recovery at the observed depth distribution ===\n")
  print(sd_recovery %>% mutate(across(where(is.numeric), ~round(., 3))))
  write_csv(sd_recovery, "results/sensitivity/sd_visit_recovery.csv")
}

# -----------------------------------------------------------------------------
# 10.4  Sampling depth as a covariate
# The review's pathway is: colonies picked tracks clinical concern, which tracks
# prescribing, so depth manufactures drug -> resistance associations. That has
# two testable legs. 10.4a tests both legs against the data instead of asserting
# them; 10.4b refits the composite model with n_isolates to see if the drug
# coefficients actually move.
# -----------------------------------------------------------------------------

# ---- 10.4a  is depth actually confounded with exposure or severity? ---------
# These numbers used to live in a comment. They justify skipping the review's
# confounding pathway, so they are computed here and written out.

depth_rx <- pa_visit_full %>% filter(!is.na(aminoglycoside_prescribed))
depth_wilcox <- suppressWarnings(
  wilcox.test(n_isolates ~ n_drug_classes > 0, data = depth_rx))

depth_by_rx <- depth_rx %>%
  group_by(any_rx = n_drug_classes > 0) %>%
  summarise(n = n(), mean_colonies = mean(n_isolates), .groups = "drop")

# The any-vs-none split above is nearly vacuous: almost every visit with
# prescription data got SOME antibiotic, so one arm is tiny. The informative
# versions are the correlation with HOW MANY classes were prescribed, and the
# per-class comparison. Report those as the actual test of leg 1.
depth_cor <- suppressWarnings(
  cor.test(depth_rx$n_isolates, depth_rx$n_drug_classes, method = "spearman"))

depth_per_class <- purrr::map_dfr(
  c("aminoglycoside", "carbapenem", "cephalosporin", "fluoroquinolones", "colistin"),
  function(cl) {
    col <- paste0(cl, "_prescribed")
    d <- depth_rx %>% filter(!is.na(.data[[col]]))
    if (n_distinct(d[[col]]) < 2) return(NULL)
    tibble(class = cl,
           n_rx      = sum(d[[col]] == 1),
           n_norx    = sum(d[[col]] == 0),
           mean_col_rx   = mean(d$n_isolates[d[[col]] == 1]),
           mean_col_norx = mean(d$n_isolates[d[[col]] == 0]),
           p = suppressWarnings(
             wilcox.test(d$n_isolates ~ d[[col]] == 1)$p.value))
  })

depth_by_status <- pa_visit_full %>%
  filter(!is.na(clinical_status)) %>%
  group_by(clinical_status) %>%
  summarise(n = n(), mean_colonies = mean(n_isolates), .groups = "drop")

depth_by_R <- pa_visit_full %>%
  filter(!is.na(any_resistance)) %>%
  mutate(depth = pmin(n_isolates, 4L)) %>%
  group_by(depth) %>%
  summarise(n = n(), pct_any_R = 100 * mean(any_resistance), .groups = "drop")

cat("\n=== 10.4a leg 1: is sampling depth associated with PRESCRIBING? ===\n")
print(depth_by_rx %>% mutate(mean_colonies = round(mean_colonies, 2)))
cat("  Wilcoxon p =", signif(depth_wilcox$p.value, 3),
    " (NOTE: one arm has n =", min(depth_by_rx$n), "- near-vacuous, see below)\n")
cat("\n  Informative version - colonies vs NUMBER of classes prescribed:\n")
cat("    Spearman rho =", round(depth_cor$estimate, 3),
    ", p =", signif(depth_cor$p.value, 3), "\n")
cat("\n  Per-class - are more colonies picked when a class is prescribed?\n")
print(depth_per_class %>% mutate(across(starts_with("mean_col"), ~round(., 2)),
                                 p = signif(p, 3)))
cat("\n=== 10.4a leg 1b: with CLINICAL STATUS? ===\n")
print(depth_by_status %>% mutate(mean_colonies = round(mean_colonies, 2)))
cat("\n=== 10.4a leg 2: does P(any resistance) rise with depth? ===\n")
print(depth_by_R %>% mutate(pct_any_R = round(pct_any_R, 1)))
cat("\n  Leg 2 (mechanical inflation) is real. Leg 1 is WEAKLY present: the\n")
cat("  any-Rx-vs-none split is near-vacuous, but colonies do track the number\n")
cat("  of classes prescribed (Spearman rho ~ 0.13), and colistin in particular\n")
cat("  is prescribed at visits where more colonies were picked. So the review's\n")
cat("  confounding pathway is present, if weak. 10.4b is therefore the test\n")
cat("  that matters: does adjusting for depth actually move the coefficients?\n")

write_csv(
  bind_rows(
    depth_by_rx %>% transmute(check = "by prescription",
                              group = if_else(any_rx, "any Rx", "no Rx"), n,
                              value = mean_colonies, metric = "mean colonies"),
    depth_by_status %>% transmute(check = "by clinical status",
                                  group = as.character(clinical_status), n,
                                  value = mean_colonies, metric = "mean colonies"),
    depth_by_R %>% transmute(check = "P(any R) by depth",
                             group = paste0("depth ", depth), n,
                             value = pct_any_R, metric = "% any resistance")
  ) %>% mutate(wilcoxon_p_anyrx = signif(depth_wilcox$p.value, 4),
               spearman_rho_nclasses = round(unname(depth_cor$estimate), 4),
               spearman_p_nclasses   = signif(depth_cor$p.value, 4)),
  "results/sensitivity/depth_confounding_checks.csv")
write_csv(depth_per_class, "results/sensitivity/depth_by_class.csv")

# ---- 10.4b  refit with depth adjusted ---------------------------------------
depth_df <- route_df %>% mutate(n_isolates_c = n_isolates - mean(n_isolates))

m_depth_adj <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance + prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_aminoglycoside_prescribed + prev_carbapenem_prescribed +
    prev_cephalosporin_prescribed + prev_fluoroquinolones_prescribed +
    prev_colistin_prescribed + n_isolates_c + (1 | individual_id),
  data = depth_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/sensitivity/m2_binary_depthadj", file_refit = "on_change"
)

depth_compare <- bind_rows(
  as_draws_df(m2_binary)   %>% as_tibble() %>% mutate(model = "Original"),
  as_draws_df(m_depth_adj) %>% as_tibble() %>% mutate(model = "+ n_isolates")
) %>%
  dplyr::select(model, starts_with("b_prev_")) %>%
  pivot_longer(-model, names_to = "term", values_to = "value") %>%
  group_by(model, term) %>%
  summarise(median = median(value),
            lo = quantile(value, 0.025), hi = quantile(value, 0.975), .groups = "drop") %>%
  mutate(term = str_remove(term, "^b_prev_"))

cat("\n=== 10.4 RESULT: do drug coefficients move when depth is adjusted? ===\n")
print(depth_compare %>%
        dplyr::select(term, model, median) %>%
        pivot_wider(names_from = model, values_from = median) %>%
        mutate(shift = `+ n_isolates` - Original) %>%
        arrange(desc(abs(shift))))

depth_eff <- as_draws_df(m_depth_adj) %>% as_tibble() %>% pull(b_n_isolates_c)
cat("\n  n_isolates coefficient: ",
    sprintf("%.3f [%.3f, %.3f], pd = %.3f\n", median(depth_eff),
            quantile(depth_eff, 0.025), quantile(depth_eff, 0.975),
            max(mean(depth_eff > 0), mean(depth_eff < 0))))
write_csv(depth_compare, "results/sensitivity/depth_adjustment_comparison.csv")

fig_depth <- depth_compare %>%
  ggplot(aes(x = median, y = fct_rev(term), colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.5), linewidth = 0.8) +
  scale_colour_manual(values = c("Original" = unname(pal["sand"]),
                                 "+ n_isolates" = unname(pal["blue"])), name = NULL) +
  labs(title = "Does adjusting for sampling depth move the drug effects?",
       subtitle = "Composite resistance model, same rows, with and without n_isolates",
       x = "Posterior median (95% CrI)", y = NULL) +
  theme(legend.position = "bottom")

ggsave("figures/sensitivity/fig_depth_adjustment.pdf", fig_depth,
       width = 10, height = 6, dpi = 300)

cat("\n=== SECTION 10 COMPLETE ===\n")
cat("  results/sensitivity/  colony_denominators.csv, denominator_comparison.csv,\n")
cat("                        variance_isolate_level.csv, depth_adjustment_comparison.csv\n")
cat("  figures/sensitivity/  fig_denominator_comparison.pdf,\n")
cat("                        fig_variance_isolate_level.pdf, fig_depth_adjustment.pdf\n")


# #############################################################################
# SECTION 11: REVIEW RESPONSES II  -  timing, pooling, calendar, measurement
#
# Second round of review responses. Like Section 10, everything here is
# ADDITIONAL: Sections 1-9 keep the original specification so the old numbers
# stay quotable, and each sub-section refits or re-runs alongside them.
#
#   11.1  CCM across patient boundaries  (#9)
#   11.2  Calendar time / modulator era  (#7)
#   11.3  Pooled cross-resistance        (#8)
#   11.4  Continuous-time Markov         (#10)
#   11.5  Measurement error, logit scale (#11)
#
# Two related fixes were made IN PLACE rather than here, because they were
# outright defects rather than alternatives:
#   - the S-map lag is now taken within patient (Section 7.4)
#   - extract_effects() reports probability of direction instead of stars, and
#     the pd > 0.90 tier no longer counts as evidence (AMR_common.R)
# #############################################################################

dir.create("results/review2", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/review2", recursive = TRUE, showWarnings = FALSE)
dir.create("models/review2",  recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 11.1  CCM that cannot embed across patients                            (#9)
#
# edm_df is 92 patient series stacked into one frame, and rEDM was handed a
# single library spanning all of it, so an embedding vector could be built from
# the end of one patient and the start of the next. rEDM 1.15's CCM() has no
# lib argument, so multi-block libraries are not available - but rEDM DROPS any
# embedding vector containing NaN. Inserting NA rows between patients therefore
# blocks exactly the crossing vectors and nothing else.
#
# Verified on a synthetic 3-segment series: usable Simplex predictions fell
# 58 -> 47 (the boundary-spanning ones) and rho moved 0.53 -> 0.33.
# -----------------------------------------------------------------------------

# NB: build from edm_endo (the ag_R/carb_R/... frame that Section 7.5 actually
# cross-maps), NOT edm_df, which holds the differently-named drug-level columns.
# Fixed at 9 rather than derived from E_map: E_map came from find_E() on the
# boundary-CROSSING frame, so using it here would size the barrier from the
# quantity we are about to re-estimate. 9 >= the maxE of 8 that find_E searches,
# which is enough to block every crossing vector.
sep_n <- 9L

# split() rather than group_modify(): group_modify re-fills the grouping column
# on every row it returns, which would hand the separator rows a real
# individual_id and defeat the whole point.
edm_df_sep <- do.call(rbind, lapply(split(edm_endo, edm_endo$individual_id), function(blk) {
  blank <- blk[rep(NA_integer_, sep_n), , drop = FALSE]
  rbind(blk, blank)
}))
edm_df_sep <- edm_df_sep[seq_len(nrow(edm_df_sep) - sep_n), , drop = FALSE]
edm_df_sep$time <- seq_len(nrow(edm_df_sep))
# rEDM uses the first column as its time index. individual_id was sitting there,
# which does not change the embedding (rEDM embeds by row order) but will throw
# "Failed to find column" the moment a requested variable is first. Be explicit.
edm_df_sep <- edm_df_sep[, c("time", setdiff(names(edm_df_sep), "time"))]
rownames(edm_df_sep) <- NULL

cat("\n=== 11.1 patient-separated EDM frame ===\n")
cat("  original rows:", nrow(edm_endo), " | separated rows:", nrow(edm_df_sep),
    " | separator width:", sep_n, "\n")
cat("  NA separator rows inserted:", sum(is.na(edm_df_sep$individual_id)), "\n")
cat("  usable ag_R values:", sum(!is.na(edm_df_sep$ag_R)),
    "(edm_endo has", sum(!is.na(edm_endo$ag_R)), ")\n")

# Re-select E on the separated frame. Section 7.3/7.5 chose E via find_E() on
# the concatenated frame, so those embeddings crossed patients too - fixing the
# CCM without fixing E would leave the boundary problem in the parameter that
# every CCM call depends on.
E_map_sep <- sapply(names(vars5), function(v) find_E_sep(edm_df_sep, v))
cat("\n  optimal E, boundary-crossing frame (used by 7.5):\n"); print(unlist(E_map))
cat("  optimal E, patient-separated frame (used below):\n");     print(E_map_sep)
if (!identical(unlist(E_map), unlist(E_map_sep)))
  cat("  -> E CHANGED once boundaries were closed; 7.5 embedded at the wrong E.\n")
E_map <- as.list(E_map_sep)   # everything downstream (11.1, 12.1) uses the safe E

ccm_pairs_sep <- list(
  c("ag_R","carb_R"), c("ag_R","ceph_R"), c("ag_R","fq_R"),
  c("carb_R","ceph_R"), c("carb_R","fq_R"), c("ceph_R","fq_R"),
  c("mucoid","ag_R"), c("mucoid","carb_R"), c("mucoid","ceph_R"), c("mucoid","fq_R")
)

run_ccm_pair_sep <- function(df, a, b, E_map) {
  d <- as.data.frame(df)                    # keep NA rows: they are the barrier
  ok <- sum(!is.na(d[[a]]) & !is.na(d[[b]]))
  if (ok < 80) { cat("  skip", a, "/", b, "- usable N =", ok, "\n"); return(NULL) }
  lib_max   <- nrow(d) - 1
  lib_sizes <- unique(round(seq(25, lib_max, length.out = 8)))
  lib_sizes <- lib_sizes[lib_sizes >= 10 & lib_sizes <= lib_max]
  Ea <- E_map[[a]]; Eb <- E_map[[b]]
  if (length(Ea) == 0 || is.na(Ea)) Ea <- 3L
  if (length(Eb) == 0 || is.na(Eb)) Eb <- 3L
  c_ab <- tryCatch(suppressWarnings(CCM(dataFrame = d, columns = b, target = a, E = Eb,
                       libSizes = lib_sizes, sample = 100)), error = function(e) NULL)
  c_ba <- tryCatch(suppressWarnings(CCM(dataFrame = d, columns = a, target = b, E = Ea,
                       libSizes = lib_sizes, sample = 100)), error = function(e) NULL)
  if (is.null(c_ab) || is.null(c_ba)) return(NULL)
  bind_rows(
    tibble(LibSize = c_ab$LibSize, rho = c_ab[[paste0(b, ":", a)]], driver = a, affected = b),
    tibble(LibSize = c_ba$LibSize, rho = c_ba[[paste0(a, ":", b)]], driver = b, affected = a)
  )
}

cat("\nRunning boundary-safe CCM on", length(ccm_pairs_sep), "pairs...\n")
ccm_all_sep <- purrr::map_dfr(ccm_pairs_sep, function(p)
  run_ccm_pair_sep(edm_df_sep, p[1], p[2], E_map))

if (!is.null(ccm_all_sep) && nrow(ccm_all_sep) > 0) {
  ccm_sum_sep <- ccm_all_sep %>%
    filter(!is.na(rho)) %>%
    group_by(driver, affected) %>%
    summarise(rho_start = rho[which.min(LibSize)],
              rho_end   = rho[which.max(LibSize)],
              delta     = rho_end - rho_start, .groups = "drop") %>%
    mutate(converges    = delta > 0.05 & rho_end > 0.1,
           driver_lab   = unname(vars5[driver]),
           affected_lab = unname(vars5[affected]))

  ccm_join <- ccm_sum %>%
    dplyr::select(driver, affected, rho_end_orig = rho_end, converges_orig = converges) %>%
    full_join(ccm_sum_sep %>%
                dplyr::select(driver, affected, rho_end_sep = rho_end, converges_sep = converges),
              by = c("driver", "affected")) %>%
    mutate(rho_shift = rho_end_sep - rho_end_orig)

  cat("\n=== 11.1 RESULT: boundary-safe vs original CCM ===\n")
  print(ccm_join %>% arrange(desc(abs(rho_shift))), n = 25)
  cat("\n  converging links, original :", sum(ccm_join$converges_orig, na.rm = TRUE), "\n")
  cat("  converging links, boundary-safe:", sum(ccm_join$converges_sep, na.rm = TRUE), "\n")
  cat("  mean |rho shift|:", round(mean(abs(ccm_join$rho_shift), na.rm = TRUE), 4), "\n")

  write_csv(ccm_sum_sep, "results/review2/ccm_summary_patientsafe.csv")
  write_csv(ccm_join,    "results/review2/ccm_boundary_comparison.csv")

  fig_ccm_bound <- ccm_join %>%
    filter(!is.na(rho_end_orig), !is.na(rho_end_sep)) %>%
    mutate(pair = paste(driver, "->", affected)) %>%
    pivot_longer(c(rho_end_orig, rho_end_sep), names_to = "frame", values_to = "rho") %>%
    mutate(frame = recode(frame, rho_end_orig = "Original (boundaries open)",
                                 rho_end_sep  = "Patient-separated")) %>%
    ggplot(aes(x = rho, y = fct_rev(pair), colour = frame)) +
    geom_point(size = 3, alpha = 0.85) +
    geom_vline(xintercept = 0.1, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = c("Original (boundaries open)" = unname(pal["sand"]),
                                   "Patient-separated" = unname(pal["red"])), name = NULL) +
    labs(title = "CCM cross-map skill with and without patient boundaries",
         subtitle = "Dashed line = the rho > 0.1 leg of the convergence rule",
         x = "rho at maximum library size", y = NULL) +
    theme(legend.position = "bottom")
  ggsave("figures/review2/fig_ccm_boundary.pdf", fig_ccm_bound, width = 10, height = 6, dpi = 300)
} else {
  cat("\n11.1: no CCM pairs returned usable output on the separated frame.\n")
}

# -----------------------------------------------------------------------------
# 11.2  Calendar time                                                    (#7)
#
# DATE SPAN - settled, do not restate from memory.
# analytical_data_prep.R sets DROP_PRE_2017 <- TRUE, which removes the 2012 and
# 2015 stragglers upstream. After that filter EVERY frame - pa_visit_full,
# pa_isolates and the all-organism clinical_microbe frame - spans
# 2017-01-03 to 2022-12-22. There are ZERO pre-2017 visits anywhere, so the
# "13 visits predate 2017" and "one visit reaches 2011" notes describe the
# PRE-filter data and are stale. Neither 2011-2022 nor 2015-2023 is the span of
# anything analysed here; the poster, the perspective manuscript and the CFF
# progress report all need 2017-2022. The span is asserted below rather than
# described, so this comment cannot go stale again.
# What IS visible is a collapse in sampling from 2020, coinciding with both
# elexacaftor-tezacaftor-ivacaftor and COVID; era is coded on that break.
#
# Scope: era is applied to the composite drug model (route_df) AND to M1
# (pa_ecological). It is NOT applied to the Markov chain, the variance
# decomposition, or CCM - those remain uncorrected for calendar time.
# -----------------------------------------------------------------------------

cat("\n=== 11.2 visits per year, by frame ===\n")
yr_tab <- bind_rows(
  pa_visit_full %>% transmute(frame = "pa_visit_full (PA visits)",
                              year = as.integer(format(as.Date(date_collected_from_subject), "%Y"))),
  pa_ecological %>% transmute(frame = "pa_ecological (M1)",
                              year = as.integer(format(as.Date(date_collected_from_subject), "%Y"))),
  route_df      %>% transmute(frame = "route_df (drug models)",
                              year = as.integer(format(as.Date(date_collected_from_subject), "%Y")))
) %>% count(frame, year) %>% pivot_wider(names_from = year, values_from = n, values_fill = 0)
print(yr_tab)
# Assert the span rather than describing it. If DROP_PRE_2017 is ever flipped
# off upstream, this stops the run instead of quietly changing every date claim.
.span_frames <- list(pa_visit_full = pa_visit_full$date_collected_from_subject,
                     pa_isolates   = pa_isolates$date_collected_from_subject)
for (.nm in names(.span_frames)) {
  .d <- as.Date(.span_frames[[.nm]]); .npre <- sum(as.integer(format(.d, "%Y")) < 2017, na.rm = TRUE)
  cat(sprintf("  %-14s span %s to %s | pre-2017 visits: %d\n",
              .nm, min(.d, na.rm = TRUE), max(.d, na.rm = TRUE), .npre))
  if (.npre > 0) stop(sprintf(
    "%d pre-2017 rows in %s - DROP_PRE_2017 is off upstream, and every date claim in the poster/manuscript is now wrong. Fix analytical_data_prep.R or update the claims.",
    .npre, .nm))
}
cat("  -> cite 2017-2022 everywhere. NOT 2011-2022, NOT 2015-2023.\n")
write_csv(yr_tab, "results/review2/visits_per_year_by_frame.csv")

add_era <- function(d) d %>%
  mutate(year = as.integer(format(as.Date(date_collected_from_subject), "%Y")),
         era  = factor(if_else(year >= 2020, "2020+", "pre-2020"),
                       levels = c("pre-2020", "2020+")))

era_df     <- add_era(route_df)
era_eco_df <- add_era(pa_ecological) %>%
  filter(!is.na(n_resistant), !is.na(n_tested), n_tested > 0,
         !is.na(prev_any_resistance), !is.na(prev_has_mrsa),
         !is.na(prev_has_mssa), !is.na(prev_has_any_mucoid))

cat("\n=== 11.2 era coding ===\n")
cat("  drug model (route_df):\n");     print(era_df %>% count(era))
cat("  ecological model (M1):\n");     print(era_eco_df %>% count(era))

m2_era <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance + prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_aminoglycoside_prescribed + prev_carbapenem_prescribed +
    prev_cephalosporin_prescribed + prev_fluoroquinolones_prescribed +
    prev_colistin_prescribed + era + (1 | individual_id),
  data = era_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/review2/m2_era", file_refit = "on_change"
)

era_compare <- bind_rows(
  as_draws_df(m2_binary) %>% as_tibble() %>% mutate(model = "Original"),
  as_draws_df(m2_era)    %>% as_tibble() %>% mutate(model = "+ era")
) %>%
  dplyr::select(model, starts_with("b_prev_")) %>%
  pivot_longer(-model, names_to = "term", values_to = "value") %>%
  group_by(model, term) %>%
  summarise(median = median(value), lo = quantile(value, .025),
            hi = quantile(value, .975), .groups = "drop") %>%
  mutate(term = str_remove(term, "^b_prev_"))

era_eff <- as_draws_df(m2_era) %>% as_tibble() %>% pull(`b_era2020P`)
cat("\n=== 11.2 RESULT: does the drug story survive calendar time? ===\n")
cat("  era (2020+ vs pre-2020): ",
    sprintf("%.3f [%.3f, %.3f], pd = %.3f\n", median(era_eff),
            quantile(era_eff, .025), quantile(era_eff, .975),
            max(mean(era_eff > 0), mean(era_eff < 0))))
print(era_compare %>% dplyr::select(term, model, median) %>%
        pivot_wider(names_from = model, values_from = median) %>%
        mutate(shift = `+ era` - Original) %>% arrange(desc(abs(shift))))
write_csv(era_compare, "results/review2/era_comparison.csv")

fig_era <- era_compare %>%
  ggplot(aes(x = median, y = fct_rev(term), colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .5), linewidth = .8) +
  scale_colour_manual(values = c("Original" = unname(pal["sand"]),
                                 "+ era" = unname(pal["plum"])), name = NULL) +
  labs(title = "Do the drug effects survive adjustment for calendar era?",
       subtitle = "Composite resistance model, same rows, with and without a 2020+ indicator",
       x = "Posterior median (95% CrI)", y = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/review2/fig_era_adjustment.pdf", fig_era, width = 10, height = 6, dpi = 300)

# ---- 11.2b  the same check on M1 (ecological, full visit span) --------------
# route_df is the prescription subset; M1 runs on pa_ecological, which carries
# the earlier visits. If a calendar effect exists anywhere it should show here.

m1_era <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance + prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    era + (1 | individual_id),
  data = era_eco_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/review2/m1_era", file_refit = "on_change"
)

m1_base <- brm(
  n_resistant | trials(n_tested) ~
    prev_any_resistance + prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    (1 | individual_id),
  data = era_eco_df, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/review2/m1_base_sameRows", file_refit = "on_change"
)

m1_era_cmp <- bind_rows(
  as_draws_df(m1_base) %>% as_tibble() %>% mutate(model = "M1 original"),
  as_draws_df(m1_era)  %>% as_tibble() %>% mutate(model = "M1 + era")
) %>%
  dplyr::select(model, starts_with("b_prev_")) %>%
  pivot_longer(-model, names_to = "term", values_to = "value") %>%
  group_by(model, term) %>%
  summarise(median = median(value), lo = quantile(value, .025),
            hi = quantile(value, .975), .groups = "drop") %>%
  mutate(term = str_remove(term, "^b_prev_"))

m1_era_eff <- as_draws_df(m1_era) %>% as_tibble() %>% pull(`b_era2020P`)
cat("\n=== 11.2b RESULT: calendar era in the ECOLOGICAL model ===\n")
cat("  N =", nrow(era_eco_df), "visits,", n_distinct(era_eco_df$individual_id), "patients\n")
cat("  era (2020+ vs pre-2020): ",
    sprintf("%.3f [%.3f, %.3f], pd = %.3f\n", median(m1_era_eff),
            quantile(m1_era_eff, .025), quantile(m1_era_eff, .975),
            max(mean(m1_era_eff > 0), mean(m1_era_eff < 0))))
print(m1_era_cmp %>% dplyr::select(term, model, median) %>%
        pivot_wider(names_from = model, values_from = median) %>%
        mutate(shift = `M1 + era` - `M1 original`) %>% arrange(desc(abs(shift))))
write_csv(m1_era_cmp, "results/review2/era_comparison_m1.csv")


# -----------------------------------------------------------------------------
# 11.3  Cross-resistance, partially pooled                               (#8)
#
# Sections 1.3-1.4 fit 16 separate models, one per (outcome class, cross-drug)
# cell, and report every cell. Noise cells get the same visual weight as real
# ones. Here the same 16 cells are ONE model with a partially pooled cross-drug
# coefficient per cell, so isolated cells shrink toward the grand mean unless
# the data pull them away.
# -----------------------------------------------------------------------------

out_specs <- tribble(
  ~outcome_class, ~yvar,           ~prevvar,
  "Carbapenem",   "carb_resistant", "prev_carb_any_resistant",
  "Aminoglycoside","ag_resistant",  "prev_ag_any_resistant",
  "Fluoroquinolone","fq_resistant", "prev_fq_any_resistant",
  "Cephalosporin","ceph_resistant", "prev_ceph_any_resistant"
)
rx_specs <- tribble(
  ~rx_class,        ~rxvar,
  "Carbapenem",     "prev_carbapenem_prescribed",
  "Aminoglycoside", "prev_aminoglycoside_prescribed",
  "Fluoroquinolone","prev_fluoroquinolones_prescribed",
  "Cephalosporin",  "prev_cephalosporin_prescribed"
)

pool_long <- purrr::pmap_dfr(out_specs, function(outcome_class, yvar, prevvar) {
  purrr::pmap_dfr(rx_specs, function(rx_class, rxvar) {
    if (!all(c(yvar, prevvar, rxvar) %in% names(pa_cross))) return(NULL)
    pa_cross %>%
      transmute(individual_id,
                y            = as.integer(.data[[yvar]]),
                prev_same    = as.numeric(.data[[prevvar]]),
                cross_rx     = as.numeric(.data[[rxvar]]),
                prev_has_mrsa, prev_has_mssa, prev_has_any_mucoid,
                outcome_class = outcome_class,
                rx_class      = rx_class,
                same_class    = as.integer(outcome_class == rx_class)) %>%
      filter(!is.na(y), !is.na(prev_same), !is.na(cross_rx))
  })
}) %>%
  mutate(cell = paste(rx_class, "->", outcome_class))

cat("\n=== 11.3 pooled cross-resistance frame ===\n")
cat("  rows:", nrow(pool_long), "| cells:", n_distinct(pool_long$cell),
    "| patients:", n_distinct(pool_long$individual_id), "\n")

m_cross_pooled <- brm(
  y ~ prev_same + cross_rx + same_class + prev_has_mrsa + prev_has_mssa +
    prev_has_any_mucoid + (1 | individual_id) + (0 + cross_rx | cell),
  data = pool_long, family = bernoulli(), prior = priors_bernoulli,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/review2/m_cross_pooled", file_refit = "on_change"
)

cell_eff <- coef(m_cross_pooled)$cell[, , "cross_rx"] %>%
  as_tibble(rownames = "cell") %>%
  rename(pooled = Estimate, lo = Q2.5, hi = Q97.5) %>%
  mutate(width_pooled = hi - lo)

unpooled <- pool_long %>%
  group_by(cell) %>%
  group_modify(~{
    fit <- try(glm(y ~ prev_same + cross_rx, data = .x, family = binomial()), silent = TRUE)
    if (inherits(fit, "try-error")) return(tibble(unpooled = NA_real_, se = NA_real_))
    s <- summary(fit)$coefficients
    if (!"cross_rx" %in% rownames(s)) return(tibble(unpooled = NA_real_, se = NA_real_))
    tibble(unpooled = s["cross_rx", "Estimate"], se = s["cross_rx", "Std. Error"])
  }) %>% ungroup()

shrink <- cell_eff %>%
  left_join(unpooled, by = "cell") %>%
  mutate(shrinkage = unpooled - pooled)

cat("\n=== 11.3 RESULT: pooled vs cell-by-cell cross-drug coefficients ===\n")
print(shrink %>% dplyr::select(cell, unpooled, pooled, shrinkage, width_pooled) %>%
        arrange(desc(abs(shrinkage))), n = 20)
sd_cell <- as_draws_df(m_cross_pooled) %>% as_tibble() %>% pull(`sd_cell__cross_rx`)
cat("\n  between-cell SD of the cross-drug effect: ",
    sprintf("%.3f [%.3f, %.3f]\n", median(sd_cell),
            quantile(sd_cell, .025), quantile(sd_cell, .975)))
cat("  (small SD => the 16 cells are largely one common effect, and reporting\n")
cat("   them separately overstates how much they differ)\n")
write_csv(shrink, "results/review2/cross_resistance_pooled.csv")

fig_pool <- shrink %>%
  filter(!is.na(unpooled)) %>%
  pivot_longer(c(unpooled, pooled), names_to = "estimator", values_to = "beta") %>%
  mutate(estimator = recode(estimator, unpooled = "Cell-by-cell (original)",
                                       pooled   = "Partially pooled")) %>%
  ggplot(aes(x = beta, y = fct_rev(cell), colour = estimator)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.8, alpha = .85) +
  scale_colour_manual(values = c("Cell-by-cell (original)" = unname(pal["sand"]),
                                 "Partially pooled" = unname(pal["blue"])), name = NULL) +
  labs(title = "Cross-resistance: 16 separate fits vs one partially pooled model",
       subtitle = "Pooling shrinks isolated cells toward the common effect",
       x = "Cross-drug log-odds coefficient", y = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/review2/fig_cross_pooled.pdf", fig_pool, width = 10, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# 11.4  Continuous-time Markov                                          (#10)
#
# The Section 7.1 chain treats every visit-to-visit step as one step, but gaps
# run 1 to 1848 days (median 84, IQR ratio 3.9). "83-91% persistence" is
# therefore per VISIT, not per unit time, and steady states inherit that.
#
# msm is not installed, but for a two-state chain the generator has a closed
# form, so this fits by direct maximum likelihood on the actual gap lengths.
# Q = [[-a, a], [b, -b]], s = a + b:
#   P11(t) = (b + a e^{-st}) / s      P12(t) = a (1 - e^{-st}) / s
#   P21(t) = b (1 - e^{-st}) / s      P22(t) = (a + b e^{-st}) / s
# -----------------------------------------------------------------------------

ctmc_df <- pa_visit_full %>%
  arrange(individual_id, date_collected_from_subject) %>%
  group_by(individual_id) %>%
  mutate(prev_state = lag(any_resistance),
         gap        = as.numeric(date_collected_from_subject - lag(date_collected_from_subject))) %>%
  ungroup() %>%
  filter(!is.na(prev_state), !is.na(any_resistance), !is.na(gap), gap > 0)

cat("\n=== 11.4 continuous-time Markov ===\n")
cat("  transitions:", nrow(ctmc_df), "| patients:", n_distinct(ctmc_df$individual_id), "\n")
cat("  gap days: median", median(ctmc_df$gap), "| IQR",
    paste(quantile(ctmc_df$gap, c(.25, .75)), collapse = "-"),
    "| max", max(ctmc_df$gap), "\n")

ctmc_nll <- function(par, d) {
  a <- exp(par[1]); b <- exp(par[2]); s <- a + b
  e <- exp(-s * d$gap)
  p11 <- (b + a * e) / s; p12 <- 1 - p11
  p22 <- (a + b * e) / s; p21 <- 1 - p22
  p <- ifelse(d$prev_state == 0,
              ifelse(d$any_resistance == 0, p11, p12),
              ifelse(d$any_resistance == 1, p22, p21))
  -sum(log(pmax(p, 1e-12)))
}
fit_ctmc <- optim(c(log(0.01), log(0.01)), ctmc_nll, d = ctmc_df,
                  method = "Nelder-Mead", control = list(maxit = 5000))
a_hat <- exp(fit_ctmc$par[1]); b_hat <- exp(fit_ctmc$par[2])

ctmc_out <- tibble(
  quantity = c("rate S->R (per day)", "rate R->S (per day)",
               "mean sojourn in S (days)", "mean sojourn in R (days)",
               "stationary P(R)", "P(R->R) at median gap", "P(S->R) at median gap"),
  value = c(a_hat, b_hat, 1 / a_hat, 1 / b_hat, a_hat / (a_hat + b_hat),
            (a_hat + b_hat * exp(-(a_hat + b_hat) * median(ctmc_df$gap))) / (a_hat + b_hat),
            a_hat * (1 - exp(-(a_hat + b_hat) * median(ctmc_df$gap))) / (a_hat + b_hat))
)
cat("\n=== 11.4 RESULT: rates per DAY, not per visit ===\n")
print(ctmc_out %>% mutate(value = round(value, 4)))

per_visit <- ctmc_df %>%
  summarise(`P(R->R) per visit` = mean(any_resistance[prev_state == 1] == 1),
            `P(S->R) per visit` = mean(any_resistance[prev_state == 0] == 1))
cat("\n  For comparison, the per-VISIT quantities Section 7.1 reports:\n")
print(per_visit %>% mutate(across(everything(), ~round(., 3))))
cat("\n  These are not the same quantity. The per-visit numbers depend on how\n")
cat("  often patients happened to be sampled; the rates above do not.\n")
write_csv(ctmc_out, "results/review2/markov_continuous_time.csv")

gap_grid <- tibble(t = seq(1, 730, by = 5)) %>%
  mutate(`P(R->R)` = (a_hat + b_hat * exp(-(a_hat + b_hat) * t)) / (a_hat + b_hat),
         `P(S->R)` = a_hat * (1 - exp(-(a_hat + b_hat) * t)) / (a_hat + b_hat))
fig_ctmc <- gap_grid %>%
  pivot_longer(-t, names_to = "transition", values_to = "p") %>%
  ggplot(aes(t, p, colour = transition)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = median(ctmc_df$gap), linetype = "dashed", colour = "grey50") +
  annotate("text", x = median(ctmc_df$gap) + 12, y = 0.05, hjust = 0, size = 3.4,
           colour = "grey40", label = sprintf("median gap = %.0f d", median(ctmc_df$gap))) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_colour_manual(values = c("P(R->R)" = unname(pal["red"]),
                                 "P(S->R)" = unname(pal["blue"])), name = NULL) +
  labs(title = "Transition probability depends on the gap, not the visit",
       subtitle = "Two-state continuous-time Markov chain fitted to observed gap lengths",
       x = "Days since previous visit", y = "Probability") +
  theme(legend.position = "bottom")
ggsave("figures/review2/fig_markov_continuous.pdf", fig_ctmc, width = 9, height = 5.5, dpi = 300)

# -----------------------------------------------------------------------------
# 11.5  Measurement error on the logit scale                            (#11)
#
# me() assumes Gaussian error. prev_theta is a Beta posterior, which is at its
# most skewed near 0 and 1 - exactly where a cohort with 83-91% persistence
# sits, so the Gaussian assumption is worst where most of the data is.
#
# For theta ~ Beta(alpha, beta), the LOG-ODDS has, exactly:
#     E[logit theta]   = digamma(alpha) - digamma(beta)
#     Var[logit theta] = trigamma(alpha) + trigamma(beta)
# and is far closer to Gaussian than theta itself. Refit on that scale.
# -----------------------------------------------------------------------------

me_logit <- me_data %>%
  mutate(alpha = 1 + prev_n_resistant_lag,
         beta  = 1 + prev_n_tested_lag - prev_n_resistant_lag,
         prev_logit_pm  = digamma(alpha) - digamma(beta),
         prev_logit_psd = sqrt(trigamma(alpha) + trigamma(beta)))

cat("\n=== 11.5 measurement error on the logit scale ===\n")
cat("  N:", nrow(me_logit), "| visits with prev_theta_pm > 0.9:",
    sum(me_logit$prev_theta_pm > 0.9),
    sprintf("(%.0f%%)\n", 100 * mean(me_logit$prev_theta_pm > 0.9)))
cat("  probability-scale SD: median", round(median(me_logit$prev_theta_psd), 3), "\n")
cat("  logit-scale SD      : median", round(median(me_logit$prev_logit_psd), 3), "\n")

m_me_logit <- brm(
  bf(ag_n_resistant | trials(ag_n_tested) ~
       me(prev_logit_pm, prev_logit_psd) +
       prev_has_any_mucoid + prev_has_mrsa + prev_has_mssa +
       prev_aminoglycoside_prescribed + (1 | individual_id)),
  data = me_logit, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  save_pars = save_pars(latent = TRUE),
  file = "models/review2/m_me_logit", file_refit = "on_change"
)

me_scale_cmp <- bind_rows(
  as_draws_df(m_me_v5) %>% as_tibble() %>%
    dplyr::select(matches("^bsp_meprev_theta")) %>%
    rename_with(~"beta_prev") %>% mutate(scale = "Probability scale (original)"),
  as_draws_df(m_me_logit) %>% as_tibble() %>%
    dplyr::select(matches("^bsp_meprev_logit")) %>%
    rename_with(~"beta_prev") %>% mutate(scale = "Logit scale (alternative)")
) %>%
  group_by(scale) %>%
  summarise(median = median(beta_prev), lo = quantile(beta_prev, .025),
            hi = quantile(beta_prev, .975),
            pd = max(mean(beta_prev > 0), mean(beta_prev < 0)), .groups = "drop")

cat("\n=== 11.5 RESULT: prior-resistance effect under each error scale ===\n")
print(me_scale_cmp %>% mutate(across(where(is.numeric), ~round(., 3))))
cat("\n  NOTE: the two coefficients are NOT on a common scale (one is per unit\n")
cat("  of theta, the other per unit of logit theta), so compare the direction,\n")
cat("  the pd, and whether the interval excludes zero - not the magnitudes.\n")
write_csv(me_scale_cmp, "results/review2/me_scale_comparison.csv")

# -----------------------------------------------------------------------------
# 11.5b  The same model on the COLONY denominator
#
# m_me_v5 measures the prior visit's precision as prev_n_tested_lag = drug-tests.
# Three aminoglycosides on one colony are not three independent draws (within
# colony rho = 0.57, so ~1.4 effective tests per colony), so that denominator
# credits each colony with about three times the information it carries and the
# measurement-error correction shrinks too little. Refit with COLONIES as the
# denominator - the number of colonies tested for AG, and the number of those
# non-susceptible - which is the unit Section 10.2 argues for.
#
# The direction is favourable: a smaller denominator means MORE shrinkage, so
# the naive-vs-aware contrast should get stronger, not weaker.
# -----------------------------------------------------------------------------

ag_drugs_me <- c("amikacin", "tobramycin", "gentamicin")
colony_counts_me <- pa_isolates %>%
  mutate(vdate = as.Date(date_collected_from_subject),
         iso_tested  = as.integer(if_any(all_of(ag_drugs_me),
                                         ~ .x %in% c("Susceptible", "Intermediate", "Resistant"))),
         iso_nonsusc = as.integer(if_any(all_of(ag_drugs_me),
                                         ~ .x %in% c("Resistant", "Intermediate")))) %>%
  group_by(individual_id, vdate) %>%
  summarise(ag_n_iso = sum(iso_tested), ag_k_iso = sum(iso_nonsusc), .groups = "drop")

me_colony <- pa_visit_full %>%
  mutate(vdate = as.Date(date_collected_from_subject)) %>%
  left_join(colony_counts_me, by = c("individual_id", "vdate")) %>%
  arrange(individual_id, vdate) %>%
  group_by(individual_id) %>%
  mutate(prev_k_iso = lag(ag_k_iso), prev_n_iso = lag(ag_n_iso)) %>%
  ungroup() %>%
  filter(!is.na(ag_n_iso), ag_n_iso >= 1,
         !is.na(prev_k_iso), !is.na(prev_n_iso), prev_n_iso >= 1,
         !is.na(prev_has_any_mucoid), !is.na(prev_has_mrsa), !is.na(prev_has_mssa),
         !is.na(prev_aminoglycoside_prescribed)) %>%
  mutate(alpha = 1 + prev_k_iso, beta = 1 + prev_n_iso - prev_k_iso,
         prev_logit_pm  = digamma(alpha) - digamma(beta),
         prev_logit_psd = sqrt(trigamma(alpha) + trigamma(beta)))

cat("\n=== 11.5b measurement error on the COLONY denominator ===\n")
cat("  N:", nrow(me_colony), "visits |", n_distinct(me_colony$individual_id), "patients\n")
cat("  median prior colonies:", median(me_colony$prev_n_iso),
    "| median drug-tests behind the same visits:",
    median(me_colony$ag_n_tested, na.rm = TRUE), "\n")
cat("  median logit-scale measurement SD:", round(median(me_colony$prev_logit_psd), 3),
    "(vs", round(median(me_logit$prev_logit_psd), 3), "on drug-tests)\n")

m_me_colony <- brm(
  bf(ag_k_iso | trials(ag_n_iso) ~
       me(prev_logit_pm, prev_logit_psd) +
       prev_has_any_mucoid + prev_has_mrsa + prev_has_mssa +
       prev_aminoglycoside_prescribed + (1 | individual_id)),
  data = me_colony, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  save_pars = save_pars(latent = TRUE),
  file = "models/review2/m_me_colony", file_refit = "on_change"
)

me_denom_cmp <- bind_rows(
  as_draws_df(m_me_logit)  %>% as_tibble() %>%
    dplyr::select(matches("^bsp_meprev_logit")) %>% rename_with(~"beta_prev") %>%
    mutate(denominator = "Drug-tests (original)"),
  as_draws_df(m_me_colony) %>% as_tibble() %>%
    dplyr::select(matches("^bsp_meprev_logit")) %>% rename_with(~"beta_prev") %>%
    mutate(denominator = "Colonies (corrected)")
) %>%
  group_by(denominator) %>%
  summarise(median = median(beta_prev), lo = quantile(beta_prev, .025),
            hi = quantile(beta_prev, .975),
            pd = max(mean(beta_prev > 0), mean(beta_prev < 0)), .groups = "drop")

cat("\n=== 11.5b RESULT: prior-resistance effect under each denominator ===\n")
print(me_denom_cmp %>% mutate(across(where(is.numeric), ~round(., 3))))
cat("  Both are on the logit-theta scale, so these ARE comparable - unlike the\n")
cat("  probability-vs-logit comparison in 11.5.\n")
write_csv(me_denom_cmp, "results/review2/me_denominator_comparison.csv")

cat("\n=== SECTION 11 COMPLETE ===\n")
cat("  results/review2/  ccm_summary_patientsafe.csv, ccm_boundary_comparison.csv,\n")
cat("                    era_comparison.csv, cross_resistance_pooled.csv,\n")
cat("                    markov_continuous_time.csv, me_scale_comparison.csv\n")
cat("  figures/review2/  fig_ccm_boundary.pdf, fig_era_adjustment.pdf,\n")
cat("                    fig_cross_pooled.pdf, fig_markov_continuous.pdf\n")


# #############################################################################
# SECTION 12: SURROGATE TESTS AND MEASUREMENT ERROR ON CLASS-SPECIFIC TERMS
#
# Third round of review responses. As in Sections 10-11, everything is
# additional: the base result is recomputed alongside so both can be quoted.
#
#   12.1  CCM against a surrogate null   -  base (threshold only) vs surrogate
#   12.2  me() on the class-specific prior-resistance terms
# #############################################################################

dir.create("results/review3", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/review3", recursive = TRUE, showWarnings = FALSE)
dir.create("models/review3",  recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 12.1  A null distribution for CCM convergence                          (#9b)
#
# Section 11.1 declares a link when delta > 0.05 and rho_end > 0.1. That is a
# bare threshold with no null: autocorrelated, highly persistent series produce
# apparent convergence routinely, and this cohort is 83-91% persistent. Here the
# same statistic is compared against surrogate data.
#
# Two nulls, because neither alone is right for this design:
#   ebisuzaki   phase-randomised, preserves the power spectrum. The standard
#               test, but it is applied to the concatenated series and so does
#               not respect patient blocks.
#   block-shift each patient's series is circularly shifted by a random amount.
#               Preserves within-patient autocorrelation and block structure,
#               destroys the alignment between the two variables. Better suited
#               to 90 short series (median 12 visits) than phase randomisation.
#
# The driver is what gets randomised: for the a -> b direction, a is surrogated,
# so the null is "a does not drive b" with b's dynamics left intact.
#
# Cost: 6 pairs x 2 directions x N_SURR x 2 methods CCM calls, ~1.6 s each.
# Results are cached to csv; delete the file to recompute.
# -----------------------------------------------------------------------------

# 500 surrogates puts the smallest achievable p at 1/501 = 0.002, which is
# below any threshold these tests are read against. At 50 the floor was 0.020 -
# fine for a screen, but a p reported as "0.02" was then indistinguishable from
# "as small as this design can measure", which is not a number to publish.
# COST: 6 pairs x 2 directions x 500 surrogates x 2 methods x ~1.6 s is roughly
# five hours. The cache is stamped with N_SURR, so this recomputes once and then
# loads from disk. Drop to 50 only for a fast structural check of the pipeline.
N_SURR <- if (AMR_QUICK) 20 else 500
surr_cache <- "results/review3/ccm_surrogate_tests.csv"

ccm_classes4 <- c("ag_R", "carb_R", "ceph_R", "fq_R")
surr_pairs <- combn(ccm_classes4, 2, simplify = FALSE)

# circularly shift each patient's block by a random offset
block_shift <- function(x, blocks) {
  out <- x
  for (b in unique(blocks)) {
    idx <- which(blocks == b)
    if (length(idx) > 1) out[idx] <- x[idx][(seq_along(idx) + sample.int(length(idx), 1) - 1) %% length(idx) + 1]
  }
  out
}

# rebuild the separated frame from a modified endo table
make_sep <- function(df, sep_n) {
  s <- do.call(rbind, lapply(split(df, df$individual_id), function(b)
    rbind(b, b[rep(NA_integer_, sep_n), , drop = FALSE])))
  s <- s[seq_len(nrow(s) - sep_n), , drop = FALSE]
  s$time <- seq_len(nrow(s)); rownames(s) <- NULL
  as.data.frame(s)
}

ccm_rho_end <- function(d, drv, aff, E_use, lib_lo, lib_hi) {
  r <- tryCatch(suppressWarnings(CCM(dataFrame = d, columns = aff, target = drv,
                                     E = E_use, libSizes = c(lib_lo, lib_hi), sample = 30)),
                error = function(e) NULL)
  if (is.null(r)) return(c(NA_real_, NA_real_))
  col <- paste0(aff, ":", drv)
  if (!col %in% names(r)) return(c(NA_real_, NA_real_))
  c(r[[col]][1], r[[col]][nrow(r)])   # rho at small lib, rho at max lib
}

# The cache records the E used and the surrogate count. E comes from
# find_E_sep() in 11.1, so re-selecting E must invalidate these results rather
# than silently serving surrogates computed at a different embedding dimension.
surr_stamp <- paste0("E=", paste(names(E_map), unlist(E_map), sep = ":", collapse = ","),
                     ";N=", N_SURR)
surr_cached_ok <- FALSE
if (file.exists(surr_cache)) {
  surr_res <- read_csv(surr_cache, show_col_types = FALSE)
  surr_cached_ok <- "stamp" %in% names(surr_res) && all(surr_res$stamp == surr_stamp)
  if (!surr_cached_ok)
    cat("\n=== 12.1 cache STALE (E or N_SURR changed) - recomputing ===\n")
}
if (surr_cached_ok) {
  cat("\n=== 12.1 surrogate tests loaded from cache ===\n")
} else {
  cat("\n=== 12.1 running CCM surrogate tests ===\n")
  cat("  ", length(surr_pairs), "pairs x 2 directions x", N_SURR, "surrogates x 2 methods\n")
  set.seed(42)
  endo_ok  <- edm_endo %>% filter(if_all(all_of(ccm_classes4), ~ !is.na(.)))
  sep_base <- make_sep(endo_ok, sep_n)
  lib_hi   <- nrow(sep_base) - 10
  lib_lo   <- 25

  surr_res <- purrr::map_dfr(surr_pairs, function(pr) {
    purrr::map_dfr(list(c(pr[1], pr[2]), c(pr[2], pr[1])), function(dirn) {
      drv <- dirn[1]; aff <- dirn[2]
      E_use <- E_map[[aff]]; if (length(E_use) == 0 || is.na(E_use)) E_use <- 3L
      obs <- ccm_rho_end(sep_base, drv, aff, E_use, lib_lo, lib_hi)

      purrr::map_dfr(c("ebisuzaki", "block_shift"), function(meth) {
        null_end <- numeric(N_SURR); null_delta <- numeric(N_SURR)
        for (s in seq_len(N_SURR)) {
          e2 <- endo_ok
          e2[[drv]] <- if (meth == "ebisuzaki") {
            as.numeric(rEDM::SurrogateData(endo_ok[[drv]], method = "ebisuzaki", num_surr = 1))
          } else {
            block_shift(endo_ok[[drv]], endo_ok$individual_id)
          }
          rr <- ccm_rho_end(make_sep(e2, sep_n), drv, aff, E_use, lib_lo, lib_hi)
          null_end[s]   <- rr[2]
          null_delta[s] <- rr[2] - rr[1]
        }
        tibble(driver = drv, affected = aff, method = meth,
               rho_end_obs   = obs[2],
               delta_obs     = obs[2] - obs[1],
               p_rho_end     = (1 + sum(null_end   >= obs[2],          na.rm = TRUE)) / (N_SURR + 1),
               p_delta       = (1 + sum(null_delta >= obs[2] - obs[1], na.rm = TRUE)) / (N_SURR + 1),
               null_end_mean = mean(null_end, na.rm = TRUE))
      })
    })
  })
  surr_res <- surr_res %>% mutate(stamp = surr_stamp)
  write_csv(surr_res, surr_cache)
}

# ---- power diagnostic for the block-shift null -------------------------------
# Circular shift returns the IDENTICAL series for any patient block that is
# constant, so for those patients the surrogate IS the observation. That drags
# the null toward the observed value and the p-value toward 0.5, which makes a
# non-significant block-shift result uninformative rather than evidence of no
# coupling. Print the constant fraction before reading anything into it.
const_frac <- purrr::map_dfr(ccm_classes4, function(v) {
  b <- edm_endo %>% filter(!is.na(.data[[v]])) %>%
    group_by(individual_id) %>%
    summarise(n = n(), constant = dplyr::n_distinct(.data[[v]]) == 1, .groups = "drop")
  tibble(variable = v, n_blocks = nrow(b),
         pct_constant = 100 * mean(b$constant),
         pct_len1     = 100 * mean(b$n == 1))
})
cat("\n=== 12.1 block-shift power check: constant patient blocks ===\n")
print(const_frac %>% mutate(across(starts_with("pct"), ~round(., 1))))
cat("  A high constant fraction means block-shift is a CONSERVATIVE filter:\n")
cat("  passing it is meaningful, failing it is not evidence of absence.\n")
write_csv(const_frac, "results/review3/blockshift_power_check.csv")

surr_verdict <- surr_res %>%
  mutate(base_converges = delta_obs > 0.05 & rho_end_obs > 0.1,
         link = paste(driver, "->", affected)) %>%
  # 24 direction x method tests. Benjamini-Hochberg within each null, on the
  # larger of the two p-values (the test requires BOTH to clear).
  group_by(method) %>%
  mutate(p_joint = pmax(p_rho_end, p_delta),
         p_joint_bh = p.adjust(p_joint, method = "BH")) %>%
  ungroup() %>%
  mutate(surrogate_passes    = base_converges & p_rho_end < 0.05 & p_delta < 0.05,
         surrogate_passes_bh = base_converges & p_joint_bh < 0.05)

cat("\n=== 12.1 RESULT: base threshold vs surrogate-tested ===\n")
print(surr_verdict %>%
        dplyr::select(link, method, rho_end_obs, null_end_mean, p_rho_end, p_delta,
                      p_joint_bh, base_converges, surrogate_passes, surrogate_passes_bh) %>%
        arrange(link, method), n = 30)

cat("\n  links passing the BASE threshold        :",
    n_distinct(surr_verdict$link[surr_verdict$base_converges]), "\n")
for (m in unique(surr_verdict$method))
  cat("  links also passing surrogate test (", m, "):",
      n_distinct(surr_verdict$link[surr_verdict$method == m & surr_verdict$surrogate_passes]), "\n")
for (m in unique(surr_verdict$method))
  cat("  links surviving BH correction   (", m, "):",
      n_distinct(surr_verdict$link[surr_verdict$method == m & surr_verdict$surrogate_passes_bh]), "\n")
cat("\n  NOTE: with", N_SURR, "surrogates the smallest achievable p is",
    sprintf("%.3f", 1 / (N_SURR + 1)),
    "- several p-values sit at that floor, so BH is operating near its\n",
    "  resolution limit. Raise N_SURR to >= 500 for any number that leaves this file.\n")
write_csv(surr_verdict, "results/review3/ccm_surrogate_verdict.csv")

fig_surr <- surr_verdict %>%
  ggplot(aes(x = rho_end_obs, y = fct_rev(link), colour = method)) +
  geom_point(aes(x = null_end_mean), shape = 4, size = 3, colour = "grey55") +
  geom_point(size = 3.2, alpha = 0.9) +
  facet_wrap(~method) +
  scale_colour_manual(values = c(ebisuzaki = unname(pal["red"]),
                                 block_shift = unname(pal["blue"])), guide = "none") +
  labs(title = "CCM cross-map skill against a surrogate null",
       subtitle = "Coloured = observed; grey x = mean of the surrogate null",
       x = "rho at maximum library size", y = NULL)
ggsave("figures/review3/fig_ccm_surrogate.pdf", fig_surr, width = 11, height = 5.5, dpi = 300)

# ---- 12.1b  the SLIDE figure, rebuilt from the surrogate verdict ------------
# This overwrites figures/ccm/fig_ccm_coupling.pdf, the path used on slides and
# in the manuscript. Section 7.5's raw version now writes to
# fig_ccm_coupling_uncorrected.pdf instead. An edge is bold only if it clears
# the base threshold AND both surrogate nulls; anything else is faint. If
# nothing clears, the diagram renders with faint edges only and says so, which
# is the honest result rather than an absent figure.

surr_pair_edges <- surr_verdict %>%
  mutate(a_lab = unname(vars5[driver]), b_lab = unname(vars5[affected])) %>%
  mutate(a = pmin(a_lab, b_lab), b = pmax(a_lab, b_lab)) %>%
  group_by(a, b) %>%
  summarise(rho = max(rho_end_obs, na.rm = TRUE),
            n_methods = n_distinct(method),
            # BH-corrected: 24 direction x method tests were run, and the
            # uncorrected verdict is not what should drive a slide.
            n_pass    = n_distinct(method[surrogate_passes_bh]),
            n_pass_raw = n_distinct(method[surrogate_passes]),
            .groups = "drop") %>%
  mutate(converges = n_pass == n_methods & n_methods > 1)

n_bold <- sum(surr_pair_edges$converges)
cat("\n=== 12.1b slide figure ===\n")
cat("  pairs clearing BOTH surrogate nulls (BH-corrected):", n_bold,
    "of", nrow(surr_pair_edges), "\n")
cat("  same, uncorrected (NOT used for the figure):",
    sum(surr_pair_edges$n_pass_raw == surr_pair_edges$n_methods &
        surr_pair_edges$n_methods > 1), "\n")
if (n_bold > 0) {
  cat("  bold edges:",
      paste(surr_pair_edges$a[surr_pair_edges$converges], "<->",
            surr_pair_edges$b[surr_pair_edges$converges], collapse = "; "), "\n")
}

cap_txt <- if (n_bold > 0) {
  paste0("Bold = clears the convergence threshold and both surrogate nulls, ",
         "BH-corrected across ", nrow(surr_verdict), " tests (ebisuzaki and block-shift, ",
         N_SURR, " surrogates each). Faint = tested, not distinguishable from null.")
} else {
  paste0("No pair survived both surrogate nulls after BH correction across ",
         nrow(surr_verdict), " tests (ebisuzaki and block-shift, ", N_SURR,
         " surrogates each). All edges tested; none distinguishable from a null ",
         "with matched autocorrelation. Faint edges only.")
}

fig_ccm_coupling_final <- ccm_coupling_diagram(
  edges   = surr_pair_edges,
  xy      = ccm_xy,
  title   = "Dynamical coupling among resistance classes",
  caption = cap_txt
)
ggsave("figures/ccm/fig_ccm_coupling.pdf", fig_ccm_coupling_final,
       width = 6.5, height = 5, dpi = 300, bg = "transparent")
write_csv(surr_pair_edges, "results/review3/ccm_coupling_edges_final.csv")


# -----------------------------------------------------------------------------
# 12.2  Measurement error on the class-specific prior-resistance terms   (#11)
#
# Section 4's me() correction is applied only to the composite prior-resistance
# covariate. The class-specific and cross-resistance claims run on
# prev_*_any_resistant, a BINARY flag derived from as few as two drug-tests,
# and those carry no correction at all. The review's tension: if measurement
# noise is big enough to be the paper's thesis, every one of those predictors
# is attenuated, and whatever sits beside them is correspondingly inflated.
#
# Test: take the pooled cross-resistance model from 11.3 and swap the binary
# prev_same for a measurement-error-aware latent frequency on the logit scale
# (exact Beta moments, as in 11.5). If the review is right, the cross-drug
# coefficients should SHRINK once prior resistance is de-attenuated.
# -----------------------------------------------------------------------------

prev_latent <- pa_visit_full %>%
  arrange(individual_id, date_collected_from_subject) %>%
  group_by(individual_id) %>%
  mutate(across(all_of(paste0(c("ag","carb","ceph","fq"), "_n_resistant")),
                ~lag(.x), .names = "prevk_{.col}"),
         across(all_of(paste0(c("ag","carb","ceph","fq"), "_n_tested")),
                ~lag(.x), .names = "prevn_{.col}")) %>%
  ungroup() %>%
  dplyr::select(individual_id, date_collected_from_subject,
                starts_with("prevk_"), starts_with("prevn_"))

beta_logit <- function(k, n) {
  a <- 1 + k; b <- 1 + n - k
  list(pm = digamma(a) - digamma(b), psd = sqrt(trigamma(a) + trigamma(b)))
}

# rebuild the 11.3 long frame carrying dates, so per-class latent moments join
pool_long2 <- purrr::pmap_dfr(out_specs, function(outcome_class, yvar, prevvar) {
  purrr::pmap_dfr(rx_specs, function(rx_class, rxvar) {
    if (!all(c(yvar, prevvar, rxvar) %in% names(pa_cross))) return(NULL)
    pa_cross %>%
      transmute(individual_id, date_collected_from_subject,
                y         = as.integer(.data[[yvar]]),
                cross_rx  = as.numeric(.data[[rxvar]]),
                prev_bin  = as.numeric(.data[[prevvar]]),
                prev_has_mrsa, prev_has_mssa, prev_has_any_mucoid,
                outcome_class = outcome_class, rx_class = rx_class,
                same_class = as.integer(outcome_class == rx_class),
                cls  = outcome_class) %>%
      filter(!is.na(y), !is.na(prev_bin), !is.na(cross_rx))
  })
}) %>%
  mutate(cell = paste(rx_class, "->", outcome_class),
         cls_key = recode(cls, Carbapenem = "carb", Aminoglycoside = "ag",
                          Fluoroquinolone = "fq", Cephalosporin = "ceph")) %>%
  left_join(prev_latent, by = c("individual_id", "date_collected_from_subject"))

pool_long2 <- pool_long2 %>%
  mutate(prev_k = purrr::pmap_dbl(list(cls_key, prevk_ag_n_resistant, prevk_carb_n_resistant,
                                       prevk_ceph_n_resistant, prevk_fq_n_resistant),
                                  function(k, a, c_, ce, f)
                                    switch(k, ag = a, carb = c_, ceph = ce, fq = f, NA_real_)),
         prev_n = purrr::pmap_dbl(list(cls_key, prevn_ag_n_tested, prevn_carb_n_tested,
                                       prevn_ceph_n_tested, prevn_fq_n_tested),
                                  function(k, a, c_, ce, f)
                                    switch(k, ag = a, carb = c_, ceph = ce, fq = f, NA_real_))) %>%
  filter(!is.na(prev_k), !is.na(prev_n), prev_n >= 1) %>%
  mutate(prev_logit_pm  = beta_logit(prev_k, prev_n)$pm,
         prev_logit_psd = beta_logit(prev_k, prev_n)$psd)

cat("\n=== 12.2 pooled frame with class-specific latent prior resistance ===\n")
cat("  rows:", nrow(pool_long2), "| cells:", n_distinct(pool_long2$cell),
    "| patients:", n_distinct(pool_long2$individual_id), "\n")
cat("  median prev_n (drug-tests behind each prior flag):", median(pool_long2$prev_n), "\n")
cat("  median logit-scale measurement SD:", round(median(pool_long2$prev_logit_psd), 3), "\n")

m_cross_naive2 <- brm(
  y ~ prev_bin + cross_rx + same_class + prev_has_mrsa + prev_has_mssa +
    prev_has_any_mucoid + (1 | individual_id) + (0 + cross_rx | cell),
  data = pool_long2, family = bernoulli(), prior = priors_bernoulli,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/review3/m_cross_naive2", file_refit = "on_change"
)

m_cross_me <- brm(
  y ~ me(prev_logit_pm, prev_logit_psd) + cross_rx + same_class +
    prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    (1 | individual_id) + (0 + cross_rx | cell),
  data = pool_long2, family = bernoulli(), prior = priors_bernoulli,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  save_pars = save_pars(latent = TRUE),
  file = "models/review3/m_cross_me", file_refit = "on_change"
)

get_cells <- function(m, lab) {
  coef(m)$cell[, , "cross_rx"] %>% as_tibble(rownames = "cell") %>%
    transmute(cell, model = lab, beta = Estimate, lo = Q2.5, hi = Q97.5)
}
me_cells <- bind_rows(get_cells(m_cross_naive2, "Naive (binary prior)"),
                      get_cells(m_cross_me,     "ME-aware (latent prior)"))

me_shift <- me_cells %>%
  dplyr::select(cell, model, beta) %>%
  pivot_wider(names_from = model, values_from = beta) %>%
  mutate(shift = `ME-aware (latent prior)` - `Naive (binary prior)`,
         shrunk = abs(`ME-aware (latent prior)`) < abs(`Naive (binary prior)`))

cat("\n=== 12.2 RESULT: do cross-drug effects shrink once prior resistance is de-attenuated? ===\n")
print(me_shift %>% arrange(desc(abs(shift))), n = 20)
cat("\n  cells that shrank toward zero:", sum(me_shift$shrunk, na.rm = TRUE),
    "of", nrow(me_shift), "\n")
cat("  mean |cross-drug beta| naive   :",
    round(mean(abs(me_shift$`Naive (binary prior)`), na.rm = TRUE), 3), "\n")
cat("  mean |cross-drug beta| ME-aware:",
    round(mean(abs(me_shift$`ME-aware (latent prior)`), na.rm = TRUE), 3), "\n")

sd_n <- as_draws_df(m_cross_naive2) %>% as_tibble() %>% pull(`sd_cell__cross_rx`)
sd_m <- as_draws_df(m_cross_me)     %>% as_tibble() %>% pull(`sd_cell__cross_rx`)
cat("\n  between-cell SD, naive   :", sprintf("%.3f [%.3f, %.3f]\n", median(sd_n),
    quantile(sd_n, .025), quantile(sd_n, .975)))
cat("  between-cell SD, ME-aware:", sprintf("%.3f [%.3f, %.3f]\n", median(sd_m),
    quantile(sd_m, .025), quantile(sd_m, .975)))
write_csv(me_cells, "results/review3/cross_me_cells.csv")
write_csv(me_shift, "results/review3/cross_me_shift.csv")

fig_me_cross <- me_cells %>%
  ggplot(aes(x = beta, y = fct_rev(cell), colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .55),
                  linewidth = .7, fatten = 2.5) +
  scale_colour_manual(values = c("Naive (binary prior)" = unname(pal["sand"]),
                                 "ME-aware (latent prior)" = unname(pal["rose"])), name = NULL) +
  labs(title = "Cross-drug effects before and after de-attenuating prior resistance",
       subtitle = "Pooled cross-resistance model; prior resistance as a binary flag vs a latent frequency with known error",
       x = "Cross-drug log-odds coefficient", y = NULL) +
  theme(legend.position = "bottom")
ggsave("figures/review3/fig_cross_me.pdf", fig_me_cross, width = 11, height = 7, dpi = 300)

cat("\n=== SECTION 12 COMPLETE ===\n")
cat("  results/review3/  ccm_surrogate_tests.csv, ccm_surrogate_verdict.csv,\n")
cat("                    cross_me_cells.csv, cross_me_shift.csv\n")
cat("  figures/review3/  fig_ccm_surrogate.pdf, fig_cross_me.pdf\n")


# #############################################################################
# SECTION 13: MULTISPATIAL CCM  -  the published implementation
#
# Sections 11.1/12.1 hand-rolled a boundary-safe CCM on top of rEDM by padding
# patient blocks with NA. That was the right idea, but the idea has a name and
# a package: multispatial CCM (Clark et al. 2015, Ecology 96:1174-1181),
# written for exactly this design - many short replicate series, none of them
# long enough to embed alone, pooled into one attractor.
#
# The original intent behind stitching patient series together was to fill the
# manifold and estimate at the population level. That intent is correct and this
# is the method for it. What multispatial CCM does NOT do is concatenate: the
# library draws points from every patient, but a state vector is never allowed
# to span two of them, because such a vector describes a transition that
# happened in nobody. The package takes a single vector with NA between series
# (c(1,2,3,NA,1,2,3) = two series of 3), which is the format built in 11.1.
#
# What this adds over the hand-rolled version:
#   - bootstrap CIs on rho across patients, not a single point estimate
#   - ccmtest(), a convergence test against library size, in place of the
#     delta > 0.05 & rho > 0.1 threshold that had no null behind it
#   - an implementation that can be cited rather than defended
#
# Runtime: ~105 s per direction at 100 iterations, so ~20 min for 6 pairs.
# Cached to csv; delete the file to recompute.
# #############################################################################

dir.create("results/msccm", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/msccm", recursive = TRUE, showWarnings = FALSE)

MSCCM_ITER  <- 100
msccm_cache <- "results/msccm/msccm_results.csv"

# -----------------------------------------------------------------------------
# 13.1  Build the multispatial vectors (one NA between patients)
# -----------------------------------------------------------------------------

ms_vec <- function(df, col, id = "individual_id") {
  s <- split(df[[col]], df[[id]])
  as.numeric(unlist(lapply(s, function(x) c(x, NA))))
}

ms_classes <- c(ag_R = "Aminoglycoside", carb_R = "Carbapenem",
                ceph_R = "Cephalosporin", fq_R  = "Fluoroquinolone")

ms_src <- edm_endo %>% arrange(individual_id)
ms_data <- lapply(names(ms_classes), function(v) ms_vec(ms_src, v))
names(ms_data) <- names(ms_classes)

cat("\n=== 13.1 multispatial vectors ===\n")
cat("  patients (series):", n_distinct(ms_src$individual_id),
    "| vector length:", length(ms_data[[1]]),
    "| real points:", sum(!is.na(ms_data[[1]])), "\n")
cat("  series length: median", median(table(ms_src$individual_id)),
    "| range", paste(range(table(ms_src$individual_id)), collapse = "-"), "\n")

# -----------------------------------------------------------------------------
# 13.2  Embedding dimension, and whether there is a signal to embed
#
# SSR_pred_boot gives prediction skill by E across the pooled series.
# SSR_check_signal asks whether skill DECAYS with prediction horizon, which is
# the signature of state-dependent dynamics rather than autocorrelation. A flat
# profile in either is the honest answer that these series are too short and too
# persistent to reconstruct an attractor from.
# -----------------------------------------------------------------------------

ms_E <- purrr::map_dfr(names(ms_classes), function(v) {
  purrr::map_dfr(2:6, function(E) {
    r <- tryCatch(SSR_pred_boot(A = ms_data[[v]], E = E, predstep = 1, tau = 1),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    tibble(variable = v, label = unname(ms_classes[v]), E = E,
           rho = mean(r$rho, na.rm = TRUE))
  })
})

ms_bestE <- ms_E %>% group_by(variable, label) %>%
  slice_max(rho, n = 1, with_ties = FALSE) %>% ungroup()

cat("\n=== 13.2 prediction skill by embedding dimension ===\n")
print(ms_E %>% mutate(rho = round(rho, 3)) %>%
        pivot_wider(id_cols = label, names_from = E, values_from = rho,
                    names_prefix = "E="))
cat("\n  best E per variable:\n")
print(ms_bestE %>% transmute(label, E, rho = round(rho, 3)))
cat("  spread in rho across E:",
    round(max(ms_E$rho) - min(ms_E$rho), 3),
    "- a flat profile means E is barely identified.\n")
write_csv(ms_E, "results/msccm/msccm_embedding.csv")

ms_signal <- purrr::map_dfr(names(ms_classes), function(v) {
  E_use <- ms_bestE$E[ms_bestE$variable == v]
  r <- tryCatch(SSR_check_signal(A = ms_data[[v]], E = E_use, tau = 1,
                                 predsteplist = 1:8), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  # rho_pre_slope is a named length-2 vector: c(Estimate, `Pr(>|t|)`).
  tibble(variable = v, label = unname(ms_classes[v]), E = E_use,
         rho_slope = unname(r$rho_pre_slope[["Estimate"]]),
         p_slope   = unname(r$rho_pre_slope[["Pr(>|t|)"]]))
})
cat("\n=== 13.2b signal check (skill decay with prediction horizon) ===\n")
print(ms_signal %>% mutate(rho_slope = round(rho_slope, 4), p_slope = signif(p_slope, 3)))
cat("  Negative slope = skill decays with horizon, the signature of\n")
cat("  state-dependent dynamics rather than plain autocorrelation. This is the\n")
cat("  precondition for CCM meaning anything; near-zero slope would say the\n")
cat("  series are predictable only because they are persistent.\n")
cat("  Variables with a significant decay (p < 0.05):",
    sum(ms_signal$p_slope < 0.05, na.rm = TRUE), "of", nrow(ms_signal), "\n")
write_csv(ms_signal, "results/msccm/msccm_signal_check.csv")

# -----------------------------------------------------------------------------
# 13.3  Bootstrapped CCM for every pair, both directions
# -----------------------------------------------------------------------------

ms_pairs <- combn(names(ms_classes), 2, simplify = FALSE)

if (file.exists(msccm_cache)) {
  ms_res <- read_csv(msccm_cache, show_col_types = FALSE)
  cat("\n=== 13.3 multispatial CCM loaded from cache ===\n")
} else {
  cat("\n=== 13.3 running multispatial CCM (", length(ms_pairs),
      "pairs x 2 directions,", MSCCM_ITER, "iterations) ===\n")
  ms_res <- purrr::map_dfr(ms_pairs, function(pr) {
    a <- pr[1]; b <- pr[2]
    Ea <- ms_bestE$E[ms_bestE$variable == a]
    Eb <- ms_bestE$E[ms_bestE$variable == b]
    cat("  ", unname(ms_classes[a]), "<->", unname(ms_classes[b]), "...\n")
    ab <- tryCatch(CCM_boot(A = ms_data[[a]], B = ms_data[[b]], E = Ea,
                            tau = 1, iterations = MSCCM_ITER), error = function(e) NULL)
    ba <- tryCatch(CCM_boot(A = ms_data[[b]], B = ms_data[[a]], E = Eb,
                            tau = 1, iterations = MSCCM_ITER), error = function(e) NULL)
    if (is.null(ab) || is.null(ba)) return(NULL)
    tt <- tryCatch(ccmtest(ab, ba), error = function(e) c(NA_real_, NA_real_))
    bind_rows(
      tibble(pair = paste(unname(ms_classes[a]), "<->", unname(ms_classes[b])),
             direction = paste(unname(ms_classes[a]), "->", unname(ms_classes[b])),
             E = Ea, L = ab$Lobs, rho = ab$rho, sdev = ab$sdevrho,
             p_converge = unname(tt[1])),
      tibble(pair = paste(unname(ms_classes[a]), "<->", unname(ms_classes[b])),
             direction = paste(unname(ms_classes[b]), "->", unname(ms_classes[a])),
             E = Eb, L = ba$Lobs, rho = ba$rho, sdev = ba$sdevrho,
             p_converge = unname(tt[2]))
    )
  })
  write_csv(ms_res, msccm_cache)
}

# delta is measured from a STABLE starting point, not from L = min.
# CCM_boot's smallest library here is L = 5, where the estimate is unstable and
# biased low; differencing against it manufactures a positive delta (up to 0.14)
# even for curves that are flat or declining across their usable range. Taking
# the quarter-library point instead reports the shape a reader actually sees.
ms_summary <- ms_res %>%
  group_by(pair, direction, E, p_converge) %>%
  summarise(L_min = min(L), L_max = max(L),
            rho_start = rho[which.min(L)],
            rho_q25   = rho[which.min(abs(L - 0.25 * max(L)))],
            rho_end   = rho[which.max(L)],
            sd_end    = sdev[which.max(L)], .groups = "drop") %>%
  mutate(delta_from_Lmin = rho_end - rho_start,   # the flattering version
         delta           = rho_end - rho_q25,     # the honest one
         lo95  = rho_end - 1.96 * sd_end,
         hi95  = rho_end + 1.96 * sd_end,
         converges = p_converge < 0.05 & lo95 > 0)

cat("\n=== 13.3 RESULT: multispatial CCM with bootstrap CIs ===\n")
print(ms_summary %>%
        transmute(direction, E,
                  rho_q25 = round(rho_q25, 3),
                  rho_end = round(rho_end, 3),
                  `95% CI` = sprintf("[%.3f, %.3f]", lo95, hi95),
                  delta = round(delta, 3),
                  p_converge = signif(p_converge, 3),
                  converges) %>%
        arrange(desc(rho_end)), n = 20)
cat("\n  directions with convergence p < 0.05 AND CI excluding 0:",
    sum(ms_summary$converges, na.rm = TRUE), "of", nrow(ms_summary), "\n")
cat("\n  CURVE SHAPE - the part that decides whether this is causality:\n")
cat("    directions where rho DECREASES from quarter-library to full:",
    sum(ms_summary$delta < 0, na.rm = TRUE), "of", nrow(ms_summary), "\n")
cat(sprintf("    mean change over the usable range: %+.4f\n", mean(ms_summary$delta, na.rm = TRUE)))
cat(sprintf("    (same quantity measured from L = %.0f instead: %+.4f - flattering,\n",
            min(ms_summary$L_min), mean(ms_summary$delta_from_Lmin, na.rm = TRUE)))
cat("     because the smallest library is too small to estimate rho stably)\n")
cat("  High but FLAT rho is the failure mode, not the success mode: cross-map\n")
cat("  skill that is already there at a fifth of the data and does not grow is\n")
cat("  shared autocorrelation, not one variable's history encoded in the\n")
cat("  other's attractor. Cobey & Baskerville 2016 is the citation.\n")
write_csv(ms_summary, "results/msccm/msccm_summary.csv")

# -----------------------------------------------------------------------------
# 13.4  Does the published method agree with the hand-rolled one?
# -----------------------------------------------------------------------------

if (exists("surr_verdict")) {
  hand <- surr_verdict %>%
    mutate(direction = paste(unname(vars5[driver]), "->", unname(vars5[affected]))) %>%
    group_by(direction) %>%
    summarise(hand_rho = first(rho_end_obs),
              hand_pass_bh = any(surrogate_passes_bh), .groups = "drop")
  agree <- ms_summary %>%
    dplyr::select(direction, ms_rho = rho_end, ms_p = p_converge, ms_converges = converges) %>%
    inner_join(hand, by = "direction")
  cat("\n=== 13.4 published vs hand-rolled ===\n")
  print(agree %>% transmute(direction,
                            ms_rho = round(ms_rho, 3), hand_rho = round(hand_rho, 3),
                            ms_converges, hand_pass_bh) %>% arrange(desc(ms_rho)))
  cat("\n  directions where the two methods disagree:",
      sum(agree$ms_converges != agree$hand_pass_bh, na.rm = TRUE), "\n")
  write_csv(agree, "results/msccm/msccm_vs_handrolled.csv")
}

# -----------------------------------------------------------------------------
# 13.5  Convergence curves with bootstrap CIs
# -----------------------------------------------------------------------------

# Colouring by `direction` gave the legend 12 entries for a 6-panel figure, and
# it truncated at both ends into something unreadable. Each panel only ever
# holds TWO lines, so the colour only has to separate them WITHIN a panel: map
# it to which way the arrow points relative to the panel title, reuse those two
# colours across panels, and label each line with its terminal rho in its own
# colour. That drops the legend to nothing and puts the number a reader wants
# where they are already looking.
arrow_dir <- function(pair, direction) {
  ifelse(sub(" ->.*$", "", direction) == sub(" <->.*$", "", pair),
         "left -> right (as titled)", "right -> left")
}

# The static contemporaneous correlation for each pair. With E = 1 the
# embedding is just the current value, so cross-mapping degenerates to a lookup
# and rho has no way to exceed this. Drawn as a reference line, it makes the
# "the level is already explained without dynamics" point visible in the figure
# instead of only in the caption. Reused by 13.7.
stat_pairs <- combn(names(ms_classes), 2, simplify = FALSE)
stat_src <- edm_endo
stat_dem <- stat_src %>% group_by(individual_id) %>%
  mutate(across(all_of(names(ms_classes)), ~ .x - mean(.x, na.rm = TRUE))) %>% ungroup()

static_r <- purrr::map_dfr(stat_pairs, function(pp) {
  tibble(pair = paste(unname(ms_classes[pp[1]]), "<->", unname(ms_classes[pp[2]])),
         static_r_raw      = cor(stat_src[[pp[1]]], stat_src[[pp[2]]], use = "complete.obs"),
         static_r_demeaned = cor(stat_dem[[pp[1]]], stat_dem[[pp[2]]], use = "complete.obs"))
})

ms_plot <- ms_res %>% filter(!is.na(rho)) %>%
  mutate(flow = arrow_dir(pair, direction))

# terminal rho, placed top-left of each panel, one line per direction
ms_endlab <- ms_summary %>%
  mutate(flow = arrow_dir(pair, direction),
         lab  = sprintf("%s:  rho = %.2f", direction, rho_end),
         vj   = ifelse(flow == "right -> left", 3.1, 1.5))

fig_msccm <- ms_plot %>%
  ggplot(aes(x = L, y = rho, colour = flow, fill = flow)) +
  geom_ribbon(aes(ymin = rho - 1.96 * sdev, ymax = rho + 1.96 * sdev),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_hline(data = static_r, aes(yintercept = static_r_raw),
             linetype = "dotdash", colour = "grey35", linewidth = 0.45,
             inherit.aes = FALSE) +
  geom_text(data = ms_endlab, aes(x = -Inf, y = Inf, label = lab, colour = flow, vjust = vj),
            hjust = -0.04, size = 2.9, inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(~pair, ncol = 2) +
  scale_colour_manual(values = c("left -> right (as titled)" = unname(pal[["red"]]),
                                 "right -> left"             = unname(pal[["blue"]]))) +
  scale_fill_manual(values   = c("left -> right (as titled)" = unname(pal[["red"]]),
                                 "right -> left"             = unname(pal[["blue"]]))) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 4),
                     expand = expansion(mult = 0.06)) +
  labs(title = "Multispatial CCM: cross-map skill vs library size",
       subtitle = sprintf(paste0("Clark et al. 2015 implementation, %d bootstrap iterations, %d patient series. ",
                                 "Dot-dash = the pair's plain\ncontemporaneous correlation - the curves start there ",
                                 "and stay there, which is the whole result."),
                          MSCCM_ITER, n_distinct(ms_src$individual_id)),
       x = "Library size (pooled across patients)", y = "rho (cross-map skill)") +
  theme(legend.position = "none",
        panel.spacing.x = unit(0.8, "lines"),
        axis.text.x     = element_text(size = 7.5))

ggsave("figures/msccm/fig_msccm_convergence.pdf", fig_msccm,
       width = 12, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# 13.6  The coupling diagram, third time - now from multispatial CCM
#
# The same node/edge diagram already exists twice: fig_ccm_coupling_uncorrected
# (Section 7.5, raw rEDM, boundary-crossing, no null) and fig_ccm_coupling
# (Section 12.1b, boundary-safe with BH-corrected surrogate tests). This is the
# third and most defensible version: the published Clark et al. implementation,
# bootstrap CIs across patients, and ccmtest() as the convergence criterion.
#
# An edge is bold only if EITHER direction clears ccmtest at p < 0.05 with a
# bootstrap CI excluding zero. Having all three is the point: the raw version
# carries a bold beta-lactam link, and it disappears as soon as patient
# boundaries are closed, a null is applied, or the published estimator is used.
# Three independent reasons, same answer.
# -----------------------------------------------------------------------------

ms_edges <- ms_summary %>%
  tidyr::separate(pair, into = c("a", "b"), sep = " <-> ", remove = FALSE) %>%
  group_by(a, b) %>%
  summarise(rho = max(rho_end, na.rm = TRUE),
            converges = any(converges, na.rm = TRUE),
            min_p = min(p_converge, na.rm = TRUE), .groups = "drop")

n_bold_ms <- sum(ms_edges$converges)
cat("\n=== 13.6 coupling diagram from multispatial CCM ===\n")
cat("  pairs with a converging direction (ccmtest p < 0.05, CI excluding 0):",
    n_bold_ms, "of", nrow(ms_edges), "\n")
cat(sprintf("  smallest convergence p across all pairs: %.2f\n", min(ms_edges$min_p, na.rm = TRUE)))

cap_ms <- if (n_bold_ms > 0) {
  sprintf(paste("Bold = converging under multispatial CCM (Clark et al. 2015),",
                "%d bootstrap iterations. Faint = tested, no convergence."), MSCCM_ITER)
} else {
  sprintf(paste0("No pair converges. Multispatial CCM (Clark et al. 2015), %d bootstrap\n",
                 "iterations across %d patient series; smallest convergence p = %.2f.\n",
                 "Cross-map skill is high but FLAT - already maximal at a fifth of the\n",
                 "library - which is shared autocorrelation, not causality."),
          MSCCM_ITER, n_distinct(ms_src$individual_id), min(ms_edges$min_p, na.rm = TRUE))
}

fig_ccm_coupling_ms <- ccm_coupling_diagram(
  edges   = ms_edges,
  xy      = ccm_xy,
  title   = "Dynamical coupling (multispatial CCM)",
  caption = cap_ms,
  edge_label = "cross-map skill %.2f")

ggsave("figures/ccm/fig_ccm_coupling_msccm.pdf", fig_ccm_coupling_ms,
       width = 8, height = 6.5, dpi = 300, bg = "transparent")
write_csv(ms_edges, "results/msccm/msccm_coupling_edges.csv")

# -----------------------------------------------------------------------------
# 13.7  Is CCM rho anything other than the static correlation?
#
# With E = 1 the "embedding" is just the current value, so cross-mapping reduces
# to a lookup and rho should track the plain contemporaneous correlation. If it
# does, CCM on these data returns correlation with extra steps, and no amount of
# library size will change that. This is the one-panel version of the argument.
# -----------------------------------------------------------------------------

# static_r is computed in 13.5 (it is drawn as a reference line there) and
# reused here rather than recomputed.
ccm_vs_cor <- ms_summary %>% group_by(pair) %>%
  summarise(ccm_rho = max(rho_end, na.rm = TRUE), .groups = "drop") %>%
  left_join(static_r, by = "pair")

cat("\n=== 13.7 CCM cross-map skill vs plain correlation ===\n")
print(ccm_vs_cor %>% arrange(desc(ccm_rho)) %>%
        mutate(across(where(is.numeric), ~round(., 3))))
cat(sprintf("\n  cor(CCM rho, raw static r)      : %.3f\n",
            cor(ccm_vs_cor$ccm_rho, ccm_vs_cor$static_r_raw)))
cat(sprintf("  cor(CCM rho, demeaned static r) : %.3f\n",
            cor(ccm_vs_cor$ccm_rho, ccm_vs_cor$static_r_demeaned)))
cat("  A near-1 correlation means CCM is recovering the contemporaneous\n")
cat("  association and nothing else. It also answers 'maybe convergence is\n")
cat("  just slow': the LEVEL is already fully explained without any dynamics.\n")
write_csv(ccm_vs_cor, "results/msccm/ccm_vs_static_correlation.csv")

fig_ccm_vs_cor <- ccm_vs_cor %>%
  ggplot(aes(static_r_raw, ccm_rho)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3.5, colour = pal[["red"]]) +
  ggrepel::geom_text_repel(aes(label = pair), size = 3, colour = "grey30") +
  coord_equal(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  labs(title = "CCM returns the static correlation, not dynamics",
       subtitle = sprintf("r = %.2f between the two. Dashed = identity.",
                          cor(ccm_vs_cor$ccm_rho, ccm_vs_cor$static_r_raw)),
       x = "Contemporaneous correlation between the two classes",
       y = "CCM cross-map skill (rho at max library)")
ggsave("figures/msccm/fig_ccm_vs_correlation.pdf", fig_ccm_vs_cor,
       width = 7.5, height = 6, dpi = 300)

cat("\n  Three coupling diagrams now exist, in increasing order of rigour:\n")
cat("    fig_ccm_coupling_uncorrected.pdf  raw rEDM, boundaries open, no null\n")
cat("    fig_ccm_coupling.pdf              boundary-safe + BH surrogate tests\n")
cat("    fig_ccm_coupling_msccm.pdf        published estimator + ccmtest\n")
cat("  If one goes on a slide, use the third - it is the citable one.\n")

cat("\n=== SECTION 13 COMPLETE ===\n")
cat("  results/msccm/  msccm_embedding.csv, msccm_signal_check.csv,\n")
cat("                  msccm_results.csv, msccm_summary.csv, msccm_vs_handrolled.csv,\n")
cat("                  ccm_vs_static_correlation.csv\n")
cat("  figures/msccm/  fig_msccm_convergence.pdf, fig_ccm_vs_correlation.pdf,\n")
cat("                  fig_ccm_coupling_msccm.pdf\n")


# #############################################################################
# SECTION 14: SURVEILLANCE DESIGN  -  what would it take to detect the dynamics?
#
# Section 13 found no cross-class convergence anywhere. That is a null result,
# but it is a designable one: it says the current sampling scheme cannot resolve
# the dynamics, not that there are none. This section turns that into a
# specification with three levers, calibrated on the observed data:
#
#   14.1  DEPTH     colonies / drug-tests per visit  -> measurement reliability
#   14.2  INTERVAL  days between visits              -> retained autocorrelation
#   14.3  LENGTH x COHORT  visits per patient, patients -> statistical power
#   14.4  Can CCM specifically be rescued by any of it?
#
# WHY STITCHING CANNOT RESCUE CCM, stated precisely: library size and trajectory
# length are different resources, and stitching only buys the first. CCM's
# convergence criterion assumes that growing the library samples ONE attractor,
# traversed by one trajectory, more densely. Adding a patient does not do that -
# it adds breadth, not depth. Two consequences follow:
#   (i)  the nearest neighbours of any delay vector are mostly OTHER patients'
#        points, so the local prediction borrows across systems rather than
#        along an orbit. That is an ensemble assumption - that patients are
#        replicate realisations of one shared system differing only in initial
#        conditions - and with different strains, genotypes, treatment histories
#        and colonisation ages it is doubtful here.
#   (ii) temporal extent per trajectory never grows. At median 12 visits and
#        E = 3 there are ~10 usable delay vectors per patient, and no number of
#        added patients resolves a process slower than ~12 visits.
# So convergence in L stops being a diagnostic for coupling and becomes one for
# ensemble homogeneity. CCM's causal claim rests on depth; stitching buys
# breadth. See also Cobey & Baskerville 2016 (PLoS ONE 11:e0169050), who show
# CCM inferring interactions between INDEPENDENT pathogen strains that merely
# oscillate at similar frequencies - the most likely alternative explanation for
# a carbapenem<->cephalosporin link driven by shared treatment cycling.
#   14.5  Recommendation
#
# The headline from 14.2 is the one that reframes the question: within-patient
# autocorrelation decays as r(t) = a*exp(-t/tau) with tau ~ 99 days, but the
# AMPLITUDE a is only ~0.32. Sampling more often raises the retained
# autocorrelation toward 0.32 and no further. There is no interval at which
# these series look like a smooth trajectory, which is what attractor
# reconstruction needs - so "collect longer series" is not by itself the answer.
# #############################################################################

dir.create("results/design", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/design", recursive = TRUE, showWarnings = FALSE)

# design_classes is defined in Section 0 - it is used from 10.3c onward, which
# is ~1800 lines before this point.

# -----------------------------------------------------------------------------
# 14.1  DEPTH -> reliability
#
# A visit's resistance frequency is k/m from m drug-tests. Its variance is
# signal + binomial sampling noise:
#     V_obs = V_theta*(1 - 1/m) + mu(1-mu)/m
# which inverts to give V_theta, and hence reliability = V_theta/(V_theta+noise).
# NOTE this is only identified for m >= 2. With a single test per visit the
# observed proportion is 0 or 1 and signal cannot be separated from noise at
# all - which is the situation for 57% of carbapenem and fluoroquinolone visits.
# -----------------------------------------------------------------------------

rel_tab <- purrr::map_dfr(names(design_classes), function(cl) {
  p <- pa_visit_full[[paste0(cl, "_prop_resistant_conservative")]]
  m <- pa_visit_full[[paste0(cl, "_n_tested")]]
  ok <- !is.na(p) & !is.na(m) & m >= 2
  if (sum(ok) < 30) return(NULL)
  pp <- p[ok]; mm <- m[ok]
  mu <- mean(pp); Vobs <- var(pp); mbar <- mean(mm)
  Vth <- max((Vobs - mu * (1 - mu) / mbar) / (1 - 1 / mbar), 0)
  Etq <- mu * (1 - mu) - Vth
  tibble(class = unname(design_classes[cl]),
         pct_single_test = 100 * mean(m[!is.na(m)] == 1),
         mean_m = mbar, V_theta = Vth, E_theta1theta = Etq,
         reliability_now = Vth / (Vth + Etq / mbar),
         m_for_rel_0.8 = ceiling(Etq / (Vth * (1 / 0.8 - 1))),
         m_for_rel_0.9 = ceiling(Etq / (Vth * (1 / 0.9 - 1))))
})

cat("\n=== 14.1 DEPTH: measurement reliability vs tests per visit ===\n")
print(rel_tab %>% mutate(across(where(is.numeric), ~round(., 3))))
write_csv(rel_tab, "results/design/depth_reliability.csv")

# X axis in COLONIES, with a design-effect correction.
#
# The earlier version put DRUG-TESTS on the x axis and hit reliability 0.80 at
# about m = 5. Current AG practice is ~1.5 colonies x 3 drugs ~ 4.5 drug-tests,
# so as drawn the figure said current sampling was already adequate - the exact
# opposite of the recommendation it was meant to support.
#
# Two corrections. The unit is colonies, because that is what a lab can act on.
# And three aminoglycosides on ONE colony are not three independent draws: their
# within-colony correlations are 0.45 (ami-tob), 0.56 (tob-gen), 0.70 (ami-gen),
# mean rho = 0.57. The standard design effect for k correlated measurements is
#     m_eff = k / (1 + (k - 1) * rho)
# which at k = 3, rho = 0.57 gives 1.4 effective tests per colony, not 3.
# Same helper the measurement-error model in Section 4 uses, so the design
# effect quoted here and the denominator correction there cannot disagree.
ag_within_rho <- ag_within_colony_rho(pa_isolates)
k_drugs          <- 3
deff             <- 1 + (k_drugs - 1) * ag_within_rho
m_eff_per_colony <- k_drugs / deff

cat(sprintf("\n  within-colony AG drug correlation rho = %.2f -> design effect %.2f\n",
            ag_within_rho, deff))
cat(sprintf("  effective independent tests per colony: %.2f (not %d)\n",
            m_eff_per_colony, k_drugs))

cur_colonies <- pa_isolates %>%
  mutate(vd = as.Date(date_collected_from_subject),
         t = as.integer(if_any(all_of(c("amikacin", "tobramycin", "gentamicin")),
                               ~ .x %in% c("Susceptible", "Intermediate", "Resistant")))) %>%
  group_by(individual_id, vd) %>% summarise(nc = sum(t), .groups = "drop") %>%
  filter(nc >= 1) %>% pull(nc) %>% mean()

rel_curve <- tidyr::expand_grid(class = rel_tab$class, colonies = 1:20) %>%
  left_join(rel_tab %>% dplyr::select(class, V_theta, E_theta1theta), by = "class") %>%
  mutate(m_eff       = colonies * m_eff_per_colony,
         reliability = V_theta / (V_theta + E_theta1theta / m_eff),
         naive       = V_theta / (V_theta + E_theta1theta / (colonies * k_drugs)))

col_for_080 <- rel_curve %>% filter(reliability >= 0.8) %>%
  group_by(class) %>% summarise(colonies_needed = min(colonies), .groups = "drop")
cat("  colonies needed for reliability 0.80:\n"); print(col_for_080)
cat(sprintf("  current mean depth: %.2f colonies\n", cur_colonies))
write_csv(col_for_080, "results/design/colonies_for_reliability.csv")

fig_depth <- ggplot(rel_curve, aes(colonies, reliability, colour = class)) +
  geom_hline(yintercept = c(0.8, 0.9), linetype = "dashed", colour = "grey60") +
  geom_line(aes(y = naive, group = class), colour = "grey78", linewidth = 0.5) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = cur_colonies, linetype = "dotted", colour = "grey35") +
  annotate("text", x = cur_colonies + 0.3, y = 0.06, hjust = 0, size = 3.3, colour = "grey30",
           label = sprintf("current: %.1f colonies", cur_colonies)) +
  scale_x_continuous(breaks = seq(0, 20, 2)) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = "How many colonies does one visit need?",
       subtitle = sprintf(paste0("Coloured: 3 AG drugs per colony, corrected for their within-colony correlation ",
                                 "(rho = %.2f, %.1f effective tests/colony).\nGrey: the same curve if those 3 drugs ",
                                 "were independent - what the model currently assumes."),
                          ag_within_rho, m_eff_per_colony),
       x = "Colonies tested per visit",
       y = "Reliability of the visit's resistance frequency") +
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave("figures/design/fig_depth_reliability.pdf", fig_depth, width = 10, height = 6.5, dpi = 300)

# -----------------------------------------------------------------------------
# 14.2  INTERVAL -> retained autocorrelation
#
# Correlation between consecutive visits, computed WITHIN patient (each patient
# centred first). The pooled correlation is ~0.7 but the ICC is ~0.6, so most of
# that is patients differing in level, not a trajectory. It is the within-patient
# number that a manifold reconstruction depends on.
# -----------------------------------------------------------------------------

ac_pairs <- purrr::map_dfr(names(design_classes), function(cl) {
  col <- paste0(cl, "_prop_resistant_conservative")
  pa_visit_full %>%
    arrange(individual_id, date_collected_from_subject) %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(individual_id) %>%
    mutate(gap = as.numeric(as.Date(date_collected_from_subject) -
                              lag(as.Date(date_collected_from_subject))),
           xc = .data[[col]] - mean(.data[[col]]), lxc = lag(xc)) %>%
    ungroup() %>%
    filter(!is.na(lxc), !is.na(gap)) %>%
    transmute(class = unname(design_classes[cl]), gap, xc, lxc)
})

ac_bins <- ac_pairs %>%
  mutate(bin = cut(gap, c(0, 14, 28, 42, 60, 90, 150, 240, 400, Inf))) %>%
  group_by(bin) %>% filter(n() >= 60) %>%
  summarise(n = n(), gap_mid = median(gap), r = cor(xc, lxc), .groups = "drop")

ac_fit <- try(nls(r ~ a * exp(-gap_mid / tau), data = ac_bins,
                  start = list(a = 0.4, tau = 60), weights = ac_bins$n), silent = TRUE)

cat("\n=== 14.2 INTERVAL: within-patient autocorrelation vs gap ===\n")
print(ac_bins %>% mutate(r = round(r, 3)))
if (!inherits(ac_fit, "try-error")) {
  a_hat <- unname(coef(ac_fit)["a"]); tau_hat <- unname(coef(ac_fit)["tau"])
  cat(sprintf("\n  fit: r(t) = %.3f * exp(-t / %.0f days)   (tau ~ %.0f weeks)\n",
              a_hat, tau_hat, tau_hat / 7))
  cur_gap <- median(ac_pairs$gap)
  cat(sprintf("  current median interval %.0f d -> retained AC %.2f\n",
              cur_gap, a_hat * exp(-cur_gap / tau_hat)))
  cat(sprintf("  CEILING: even at t -> 0 the within-patient AC is only %.2f.\n", a_hat))
  cat("  Correcting for measurement reliability (~", round(mean(rel_tab$reliability_now), 2),
      ") puts the true ceiling near ", round(a_hat / mean(rel_tab$reliability_now), 2), ".\n", sep = "")
  cat("  So shortening the interval helps, but cannot deliver the 0.7-0.9\n")
  cat("  autocorrelation that attractor reconstruction wants. That is a\n")
  cat("  property of the biology at this timescale, not of the cohort size.\n")
  ac_grid <- tibble(gap = 0:400) %>% mutate(r = a_hat * exp(-gap / tau_hat))
  write_csv(bind_rows(ac_bins %>% mutate(src = "observed"),
                      tibble(bin = NA, n = NA_integer_, gap_mid = c(14, 28, 56, 84),
                             r = a_hat * exp(-c(14, 28, 56, 84) / tau_hat), src = "fitted")),
            "results/design/interval_autocorrelation.csv")
  fig_int <- ggplot(ac_bins, aes(gap_mid, r)) +
    geom_line(data = ac_grid, aes(gap, r), colour = pal[["sand"]], linewidth = 1) +
    geom_point(aes(size = n), colour = pal[["red"]]) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = cur_gap, linetype = "dotted") +
    annotate("text", x = cur_gap + 8, y = a_hat * 0.95, hjust = 0, size = 3.4,
             label = sprintf("current median\n%.0f days", cur_gap)) +
    scale_size_continuous(guide = "none") +
    labs(title = "Within-patient autocorrelation decays over ~14 weeks",
         subtitle = sprintf("r(t) = %.2f exp(-t/%.0f d). The amplitude, not the interval, is the ceiling.",
                            a_hat, tau_hat),
         x = "Days between consecutive visits", y = "Within-patient lag-1 correlation")
  ggsave("figures/design/fig_interval_autocorrelation.pdf", fig_int, width = 9, height = 5.5, dpi = 300)
} else {
  a_hat <- 0.32; tau_hat <- 99
  cat("  (exponential fit failed; falling back to a = 0.32, tau = 99 d)\n")
}

# ---- 14.2d  is the decay trajectory memory, or episode structure? ----------
# Since 10.3d showed sd_visit is not identified at this depth, 14.2 is now the
# ONLY positive evidence of within-patient temporal structure, so its obvious
# confound has to be closed before it carries weight.
#
# Sampling interval is not randomly assigned: short gaps cluster INSIDE
# exacerbation episodes, where the patient sits in one clinical state and the
# same lineage may simply be resampled weeks apart. That manufactures high
# short-gap correlation which is not trajectory memory.
#
# Three fits on the same data: all pairs; dropping pairs that plausibly sit in
# one episode (both visits APE and <= 30 days apart); and Stable -> Stable pairs
# only, which have no episode structure at all. If the decay survives the third,
# it is a defensible positive result and the tension with 10.3d resolves - the
# mixed model cannot see the component, the autocorrelation can.

ac_pairs_cs <- purrr::map_dfr(names(design_classes), function(cl) {
  col <- paste0(cl, "_prop_resistant_conservative")
  pa_visit_full %>%
    arrange(individual_id, date_collected_from_subject) %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(individual_id) %>%
    mutate(gap = as.numeric(as.Date(date_collected_from_subject) -
                              lag(as.Date(date_collected_from_subject))),
           xc = .data[[col]] - mean(.data[[col]]), lxc = lag(xc),
           cs_now = as.character(clinical_status),
           cs_prev = lag(as.character(clinical_status))) %>%
    ungroup() %>%
    filter(!is.na(lxc), !is.na(gap)) %>%
    transmute(class = unname(design_classes[cl]), gap, xc, lxc, cs_now, cs_prev)
})

.is_ape <- function(x) !is.na(x) & x %in% c("APE Outpatient", "APE IV Antibiotics")
ac_pairs_cs <- ac_pairs_cs %>%
  mutate(same_episode = .is_ape(cs_now) & .is_ape(cs_prev) & gap <= 30,
         both_stable  = cs_now == "Clinically Stable" & cs_prev == "Clinically Stable")

fit_ac <- function(d, label, min_n = 30) {
  b <- d %>% mutate(bin = cut(gap, c(0, 14, 28, 42, 60, 90, 150, 240, 400, Inf))) %>%
    group_by(bin) %>% filter(n() >= min_n) %>%
    summarise(n = n(), gap_mid = median(gap), r = cor(xc, lxc), .groups = "drop")
  if (nrow(b) < 3) return(tibble(subset = label, n_pairs = nrow(d), n_bins = nrow(b),
                                 a = NA_real_, tau = NA_real_))
  f <- try(nls(r ~ a * exp(-gap_mid / tau), data = b,
               start = list(a = 0.4, tau = 60), weights = b$n), silent = TRUE)
  if (inherits(f, "try-error")) return(tibble(subset = label, n_pairs = nrow(d),
                                              n_bins = nrow(b), a = NA_real_, tau = NA_real_))
  tibble(subset = label, n_pairs = nrow(d), n_bins = nrow(b),
         a = unname(coef(f)["a"]), tau = unname(coef(f)["tau"]))
}

ac_confound <- bind_rows(
  fit_ac(ac_pairs_cs, "all pairs"),
  fit_ac(ac_pairs_cs %>% filter(!same_episode), "excluding same-episode pairs"),
  fit_ac(ac_pairs_cs %>% filter(both_stable),   "Stable -> Stable only"))

cat("\n=== 14.2d autocorrelation: trajectory memory or episode structure? ===\n")
cat(sprintf("  same-episode pairs (both APE, gap <= 30d): %d of %d\n",
            sum(ac_pairs_cs$same_episode), nrow(ac_pairs_cs)))
print(ac_confound %>% mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(),
      row.names = FALSE)
write_csv(ac_confound, "results/design/interval_autocorrelation_confound.csv")
cat("\n  How to read this (measured 2026-08-27: 0.343 / 0.342 / 0.172):\n")
cat("  - The same-episode exclusion is a WEAK test - it removes only a few\n")
cat("    hundred pairs and barely moves a. Do not report it as the check.\n")
cat("  - Stable -> Stable is the informative contrast, and a roughly HALVES\n")
cat("    (0.34 -> 0.17). About half the headline amplitude is carried by pairs\n")
cat("    involving an exacerbation, so 0.34 cannot be cited as trajectory\n")
cat("    memory. The defensible number is ~0.17.\n")
cat("  - tau in the Stable -> Stable row is NOT identified (a near-flat curve\n")
cat("    fits at almost any tau on ~900 pairs). Quote a, not tau, there.\n")
cat("  - This makes the ceiling argument STRONGER: if genuine within-patient\n")
cat("    memory is ~0.17, shortening the sampling interval is even further from\n")
cat("    delivering the 0.7-0.9 that attractor reconstruction needs.\n")

# ---- 14.2b  anchor the interval to treatment cycles, not round numbers -------
# An interval is only interpretable against the timescales the biology actually
# runs on. P. aeruginosa turns over far faster than any sampling scheme can
# follow, so the relevant clock is not generation time but the PERTURBATION
# schedule: a 14-day IV course, and the 28-day on/off cycling of inhaled
# tobramycin / aztreonam. If consecutive cultures are further apart than one
# cycle, they bracket several perturbations and alias them into a single step -
# which is what makes the series look like noise around a patient mean.

cycle_landmarks <- tibble(
  landmark = c("IV course", "inhaled on/off cycle", "two cycles", "quarterly"),
  days     = c(14, 28, 56, 90)) %>%
  mutate(pct_gaps_within = purrr::map_dbl(days, ~100 * mean(ac_pairs$gap <= .x)),
         retained_AC     = a_hat * exp(-days / tau_hat),
         pct_of_ceiling  = 100 * retained_AC / a_hat)

cat("\n=== 14.2b interval against treatment cycles ===\n")
print(cycle_landmarks %>% mutate(across(where(is.numeric), ~round(., 2))))
cat(sprintf("\n  median interval %.0f d = %.1f IV courses = %.1f inhaled cycles\n",
            median(ac_pairs$gap), median(ac_pairs$gap) / 14, median(ac_pairs$gap) / 28))
cat(sprintf("  %.0f%% of consecutive-visit pairs span MORE than one inhaled cycle,\n",
            100 * mean(ac_pairs$gap > 28)))
cat("  so most observed transitions average over several treatment perturbations.\n")

# ---- 14.2c  the loophole that keeps this an open question, not a null -------
# Everything above assumes visits sample the trajectory at RANDOM phases of the
# treatment cycle. They do not. Cultures here are clinically triggered -
# exacerbation, culture, treat - so visits land at systematically SIMILAR points
# in that cycle. A process locked to the cycle and sampled once per cycle at the
# same phase aliases to apparent constancy: the quantity is varying, and every
# observation catches it in the same place.
#
# This is the temporal analogue of the single-colony problem in Section 10.3d.
# There, informative sampling made a variable quantity unidentifiable across
# colonies; here it can make one unidentifiable across visits. In both cases the
# honest statement is about what the design can SEE, not about the system.
#
# So the defensible claim is NOT "there are no fast dynamics". It is:
#   we cannot distinguish "no fast dynamics" from "fast dynamics aliased by
#   clinically-triggered sampling", and only a dedicated dense-sampling
#   substudy - fixed calendar intervals, decoupled from clinical trigger -
#   separates them.
# That is a testable open question rather than a null result, and it is the
# version of the interval argument that survives scrutiny.
if ("clinical_status" %in% names(pa_visit_full)) {
  phase_tab <- pa_visit_full %>%
    filter(!is.na(clinical_status)) %>%
    count(clinical_status) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    arrange(desc(n))
  cat("\n=== 14.2c are visits phase-locked to the treatment cycle? ===\n")
  print(phase_tab)
  cat("  Visits are triggered by clinical state, not by the calendar. If the\n")
  cat("  resistance frequency oscillates ON that same cycle, sampling at a\n")
  cat("  fixed phase of it returns a near-constant value no matter how large\n")
  cat("  the true amplitude. Constancy under informative sampling is not\n")
  cat("  evidence of a constant process.\n")
  cat("  CONSEQUENCE FOR THE INTERVAL ASK: it is not 'sample faster because we\n")
  cat("  measured fast dynamics'. It is 'sample on a FIXED SCHEDULE because the\n")
  cat("  current trigger cannot tell a flat process from an aliased one'.\n")
  write_csv(phase_tab, "results/design/visit_phase_locking.csv")
}
# ---- the empirical timescale anchor, from this cohort's own transitions -----
# Generation time is the wrong anchor: in vivo P. aeruginosa doubling in CF
# sputum is poorly constrained and varies between biofilm and planktonic states,
# so a sampling interval pinned to it is arguable. The quantity that matters is
# how long a resistant lineage takes to change frequency appreciably - and
# Section 11.4 already estimates that from the data as the CTMC mean sojourn
# times, 1/rate, in days.
if (exists("ctmc_out")) {
  soj_S <- ctmc_out$value[ctmc_out$quantity == "mean sojourn in S (days)"]
  soj_R <- ctmc_out$value[ctmc_out$quantity == "mean sojourn in R (days)"]
  cat(sprintf("\n  CTMC mean sojourn (Section 11.4): S %.0f d, R %.0f d\n", soj_S, soj_R))
  cat(sprintf("  A common rule of thumb is to sample at <= 1/4 of the fastest\n"))
  cat(sprintf("  sojourn: %.0f/4 = %.0f days.\n", soj_S, soj_S / 4))
  cat("  CAUTION, and it runs in our favour: those sojourns were estimated from\n")
  cat("  transitions observed at a median 84-day gap, so any S->R->S excursion\n")
  cat("  completing between visits is invisible and both rates are biased DOWN.\n")
  cat("  The true sojourn is shorter, so the required interval is shorter than\n")
  cat("  this estimate implies - it is an upper bound, not a target.\n")
}

cat("  Sampling once per inhaled cycle (28 d) recovers ",
    sprintf("%.0f%%", cycle_landmarks$pct_of_ceiling[cycle_landmarks$days == 28]),
    " of the attainable\n  autocorrelation against ",
    sprintf("%.0f%%", 100 * a_hat * exp(-median(ac_pairs$gap) / tau_hat) / a_hat),
    " now - and it aligns with a rhythm the clinic already runs.\n", sep = "")
write_csv(cycle_landmarks, "results/design/interval_vs_treatment_cycles.csv")

# -----------------------------------------------------------------------------
# 14.3  LENGTH x COHORT -> power
#
# Generative model calibrated to the above: latent logit-scale AR(1) with the
# achievable within-patient autocorrelation, a patient-level intercept giving
# the observed ICC, cross-class coupling beta at lag 1, and binomial observation
# at depth m. The detector is a cross-lagged mixed model - does last visit's
# class X predict this visit's class Y, over and above last visit's Y - which is
# the population-level question CCM was being asked, in a form this data shape
# can actually support.
# -----------------------------------------------------------------------------

sim_cohort <- function(n_pat, L, rho, beta, m, sd_pat = 1.0, sd_inn = 1.0) {
  purrr::map_dfr(seq_len(n_pat), function(i) {
    mux <- rnorm(1, 0, sd_pat); muy <- rnorm(1, 0, sd_pat)
    x <- y <- numeric(L); x[1] <- rnorm(1); y[1] <- rnorm(1)
    for (t in 2:L) {
      x[t] <- rho * x[t - 1] + rnorm(1, 0, sd_inn)
      y[t] <- rho * y[t - 1] + beta * x[t - 1] + rnorm(1, 0, sd_inn)
    }
    tibble(id = i, t = seq_len(L),
           kx = rbinom(L, m, plogis(mux + x)), ky = rbinom(L, m, plogis(muy + y)), m = m)
  })
}

detect_xlag <- function(d) {
  dd <- d %>% group_by(id) %>%
    mutate(py = ky / m, px = kx / m, lpy = lag(py), lpx = lag(px)) %>%
    ungroup() %>% filter(!is.na(lpy))
  fit <- try(suppressMessages(lme4::lmer(py ~ lpy + lpx + (1 | id), data = dd)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  s <- summary(fit)$coefficients
  if (!"lpx" %in% rownames(s)) return(NA_real_)
  2 * pnorm(-abs(s["lpx", "t value"]))
}

# Monte Carlo replicates per cell. At 50 the standard error on a power near 0.5
# is sqrt(.25/50) = 0.071, which is the same size as the lever differences the
# grid is used to rank - power ran 0.20 -> 0.30 -> 0.24 -> 0.14 across adjacent
# intervals in one panel, non-monotone purely from noise. At 300 the SE falls to
# 0.029, so a 0.14 gap between levers is ~3.4 SE and the ranking is readable.
# Cost is roughly 216 cells x DESIGN_REPS lmer fits; 300 is an overnight run.
DESIGN_REPS <- if (AMR_QUICK) 15 else 300
design_cache <- "results/design/power_grid_dt.csv"

# Sampling interval enters through the autocorrelation it retains. The 14.2 fit
# gives observed r(dt) = a*exp(-dt/tau); dividing by reliability recovers the
# LATENT autocorrelation the simulation should use, since the binomial
# observation step is applied separately below.
rho_at_dt <- function(dt) {
  min(max(a_hat * exp(-dt / tau_hat) / mean(rel_tab$reliability_now), 0.02), 0.95)
}
dt_levels <- c(14, 28, 56, 84)
cat("\n  latent AR(1) rho implied by each sampling interval:\n")
print(tibble(interval_days = dt_levels,
             rho_latent = round(vapply(dt_levels, rho_at_dt, numeric(1)), 3)))

# Like the 12.1 surrogate cache, this one is stamped. Raising DESIGN_REPS must
# invalidate it - otherwise the rep count changes in the source and the figure
# silently keeps reporting the old, noisier grid.
design_stamp <- paste0("reps=", DESIGN_REPS, ";dt=", paste(dt_levels, collapse = "/"))
design_cached_ok <- FALSE
if (file.exists(design_cache)) {
  power_grid <- read_csv(design_cache, show_col_types = FALSE)
  design_cached_ok <- "stamp" %in% names(power_grid) && all(power_grid$stamp == design_stamp)
  if (!design_cached_ok)
    cat("\n=== 14.3 cache STALE (reps or dt levels changed) - recomputing ===\n")
}
if (design_cached_ok) {
  cat("\n=== 14.3 power grid loaded from cache ===\n")
} else {
  cat("\n=== 14.3 simulating power grid (", DESIGN_REPS, "reps per cell ) ===\n")
  set.seed(42)
  grid <- tidyr::expand_grid(n_pat = c(30, 90, 200), L = c(12, 24, 48),
                             beta = c(0, 0.2, 0.3), m = c(3, 10), dt = dt_levels)
  power_grid <- purrr::pmap_dfr(grid, function(n_pat, L, beta, m, dt) {
    rho <- rho_at_dt(dt)
    p <- replicate(DESIGN_REPS, detect_xlag(sim_cohort(n_pat, L, rho, beta, m)))
    pw <- mean(p < 0.05, na.rm = TRUE); nok <- sum(!is.na(p))
    tibble(n_pat = n_pat, L = L, beta = beta, m = m, dt = dt, rho = rho,
           power = pw, n_ok = nok,
           # binomial MC error, carried so every downstream number can be read
           # against the noise that produced it
           power_se = sqrt(pw * (1 - pw) / nok),
           stamp = design_stamp)
  })
  write_csv(power_grid, design_cache)
}

cat("\n=== 14.3 RESULT: false-positive rate (beta = 0) ===\n")
print(power_grid %>% filter(beta == 0) %>%
        summarise(cells = n(), mean_FPR = round(mean(power), 3), max_FPR = round(max(power), 3)))

cat("\n=== 14.3 WHICH LEVER PAYS? power at beta = 0.2, m = 3 ===\n")
cat("  rows = sampling interval, columns = cohort size x visits per patient\n")
print(power_grid %>% filter(beta == 0.2, m == 3) %>%
        mutate(design = paste0("n", n_pat, "/L", L)) %>%
        pivot_wider(id_cols = dt, names_from = design, values_from = power) %>%
        rename(interval_days = dt))

cat("\n  Marginal effect of each lever, holding the others at the cohort's own\n")
cat("  operating point (n = 90, L = 12, m = 3, beta = 0.2):\n")
base_row <- power_grid %>% filter(n_pat == 90, L == 12, m == 3, beta == 0.2, dt == 84)
base_pw  <- base_row$power
base_se  <- if ("power_se" %in% names(base_row)) base_row$power_se else NA_real_
lev <- bind_rows(
  power_grid %>% filter(n_pat == 200, L == 12, m == 3, beta == 0.2, dt == 84) %>%
    transmute(lever = "cohort 90 -> 200 patients", power, power_se),
  power_grid %>% filter(n_pat == 90, L == 48, m == 3, beta == 0.2, dt == 84) %>%
    transmute(lever = "length 12 -> 48 visits", power, power_se),
  power_grid %>% filter(n_pat == 90, L == 12, m == 10, beta == 0.2, dt == 84) %>%
    transmute(lever = "depth 3 -> 10 tests", power, power_se),
  power_grid %>% filter(n_pat == 90, L == 12, m == 3, beta == 0.2, dt == 28) %>%
    transmute(lever = "interval 84 -> 28 days", power, power_se),
  power_grid %>% filter(n_pat == 90, L == 12, m == 3, beta == 0.2, dt == 14) %>%
    transmute(lever = "interval 84 -> 14 days", power, power_se)
) %>%
  # Every gain below is a difference of two Monte Carlo estimates. Reporting it
  # without its own error is how a +0.04 that is indistinguishable from zero
  # ends up quoted as a finding. se_gain makes the comparison self-checking:
  # a gain smaller than ~2 se_gain is not separable from no effect.
  mutate(baseline = base_pw,
         gain     = power - baseline,
         se_gain  = sqrt(power_se^2 + base_se^2),
         separable_from_zero = abs(gain) > 2 * se_gain) %>%
  arrange(desc(gain))
print(lev %>% mutate(across(where(is.numeric), ~round(., 3))))
cat(sprintf("\n  (%d reps/cell; MC standard error on a single power near 0.5 is %.3f)\n",
            DESIGN_REPS, sqrt(0.25 / DESIGN_REPS)))
cat("\n  READ THIS BEFORE QUOTING THE INTERVAL ROW.\n")
cat("  In this simulation the interval buys almost nothing, while length,\n")
cat("  depth and cohort all buy a lot. That contradicts the intuition that\n")
cat("  denser sampling is the dominant lever, and the reason is a limitation\n")
cat("  of how dt is modelled here, not evidence against the intuition:\n")
cat("    (1) dt enters ONLY through the latent autocorrelation rho. The\n")
cat("        per-step coupling beta is held fixed, but a real process observed\n")
cat("        twice as often has a SMALLER per-step effect - there is less time\n")
cat("        between observations for it to accumulate. Holding beta fixed\n")
cat("        while raising rho mostly adds collinearity between lagged X and\n")
cat("        lagged Y, which hurts the partial coefficient.\n")
cat("    (2) L and dt are varied independently, so no budget is held constant.\n")
cat("        The interval question is only well posed under a stated\n")
cat("        constraint. At fixed VISITS per patient, shortening dt shrinks\n")
cat("        the calendar span covered and should lose. At fixed STUDY\n")
cat("        DURATION, shortening dt buys more visits and should win. This\n")
cat("        grid answers the first; the design question is usually the\n")
cat("        second.\n")
cat("  So: the honest claim from THIS simulation is that visits-per-patient\n")
cat("  and tests-per-visit dominate cohort size. Whether denser sampling beats\n")
cat("  a larger cohort needs a fixed-duration reformulation, which is not\n")
cat("  built here. Do not cite the interval row as evidence either way.\n")
write_csv(lev, "results/design/lever_comparison.csv")

mde <- power_grid %>% filter(beta > 0, power >= 0.8) %>%
  group_by(n_pat, L, m, dt) %>% summarise(min_detectable_beta = min(beta), .groups = "drop")
cat("\n=== 14.3 minimum detectable coupling (80% power), by interval ===\n")
print(mde %>% filter(m == 3) %>%
        pivot_wider(id_cols = c(n_pat, L), names_from = dt,
                    values_from = min_detectable_beta, names_prefix = "dt="))

# The Emory operating point is marked so the figure shows where the cohort
# actually sits, not just an abstract surface.
emory <- power_grid %>% filter(n_pat == 90, L == 12, m == 3, dt == 84, beta == 0.2)

# Strip labels are built as ORDERED FACTORS, not pasted strings. As plain
# character, "n = 200" sorts before "n = 30" and "m = 10 tests" before
# "m = 3 tests", so the panels ran 200/200/30/30/90/90 - alphabetical order on a
# numeric variable. The beta strip was an unlabelled bare number for the same
# reason: factor(beta) prints "0.2", not "beta = 0.2".
fig_power <- power_grid %>% filter(beta > 0) %>%
  mutate(depth   = factor(paste0("m = ", m, " tests"),
                          levels = paste0("m = ", sort(unique(m)), " tests")),
         cohort  = factor(paste0("n = ", n_pat),
                          levels = paste0("n = ", sort(unique(n_pat)))),
         couple  = factor(paste0("beta = ", beta),
                          levels = paste0("beta = ", sort(unique(beta))))) %>%
  ggplot(aes(dt, power, colour = factor(L))) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  facet_grid(depth ~ cohort + couple, labeller = label_wrap_gen(20)) +
  # Adjacent panels ran their edge ticks together ("...2814 84 56 2814..."):
  # the rightmost "14" of one panel touching the leftmost "84" of the next.
  # Padding the panel range and widening the gutter separates them.
  scale_x_reverse(breaks = dt_levels, expand = expansion(mult = 0.15)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Power by sampling interval, cohort and series length",
       subtitle = sprintf(paste0("Cross-lagged mixed model, %d reps/cell (MC SE ~%.2f at power 0.5). ",
                                 "X axis reversed: denser sampling to the right.\n",
                                 "See the caveat in the console output before reading the ",
                                 "interval axis as a design recommendation."),
                          DESIGN_REPS, sqrt(0.25 / DESIGN_REPS)),
       x = "Days between visits", y = "Power", colour = "visits/patient") +
  theme(legend.position = "bottom",
        panel.spacing.x = unit(0.7, "lines"),
        axis.text.x     = element_text(size = 7.5))
if (nrow(emory) == 1) {
  cat(sprintf("\n  Emory operating point (n=90, L=12, m=3, dt=84, beta=0.2): power = %.2f\n",
              emory$power))
}
ggsave("figures/design/fig_power_grid.pdf", fig_power, width = 11, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# 14.4  Can CCM itself be rescued?
# Same simulated cohorts, but detected with CCM on the patient-separated frame,
# at the current design and at an enriched one.
# -----------------------------------------------------------------------------

ccm_detect_sim <- function(d, sep = 4L) {
  blocks <- split(d, d$id)
  frame <- do.call(rbind, lapply(blocks, function(b) {
    bb <- data.frame(x = b$kx / b$m, y = b$ky / b$m)
    rbind(bb, bb[rep(NA_integer_, sep), , drop = FALSE])
  }))
  frame <- frame[seq_len(nrow(frame) - sep), , drop = FALSE]
  frame$time <- seq_len(nrow(frame))
  # rEDM treats the FIRST column as its time index. With x first it consumed x
  # and then could not find it as a data column ("Failed to find column: x"),
  # which is why this returned NaN. Put time first explicitly.
  frame <- frame[, c("time", "x", "y")]
  lib_hi <- nrow(frame) - 10
  if (lib_hi < 40) return(NA)
  r <- tryCatch(suppressWarnings(
    CCM(dataFrame = frame, columns = "y", target = "x", E = 3,
        libSizes = c(25, lib_hi), sample = 30)), error = function(e) NULL)
  if (is.null(r) || !"y:x" %in% names(r)) return(NA)
  r[["y:x"]][nrow(r)] - r[["y:x"]][1]        # convergence delta
}

ccm_cache <- "results/design/ccm_power.csv"
if (file.exists(ccm_cache)) {
  ccm_power <- read_csv(ccm_cache, show_col_types = FALSE)
  cat("\n=== 14.4 CCM power loaded from cache ===\n")
} else {
  cat("\n=== 14.4 can CCM detect the same coupling? ===\n")
  set.seed(7)
  ccm_grid <- tidyr::expand_grid(
    design = c("current", "enriched"), beta = c(0, 0.3))
  ccm_power <- purrr::pmap_dfr(ccm_grid, function(design, beta) {
    # Kept deliberately modest: CCM on a 200x48 cohort is ~9600 points and each
    # call takes tens of seconds, which buys no extra insight here. n = 90 with
    # L = 48 already quadruples trajectory DEPTH, which is the lever that would
    # matter if CCM could be rescued at all.
    n_pat <- if (design == "current") 90 else 90
    L     <- if (design == "current") 12 else 48
    m     <- if (design == "current") 3  else 10
    dt    <- if (design == "current") 84 else 28
    deltas <- replicate(15, ccm_detect_sim(sim_cohort(n_pat, L, rho_at_dt(dt), beta, m)))
    tibble(design = design, n_pat = n_pat, L = L, m = m, dt = dt, beta = beta,
           mean_delta = mean(deltas, na.rm = TRUE),
           pct_delta_gt_0.05 = 100 * mean(deltas > 0.05, na.rm = TRUE))
  })
  write_csv(ccm_power, ccm_cache)
}
print(ccm_power %>% mutate(across(where(is.numeric), ~round(., 3))))
cat("  Compare the beta = 0.3 rows against beta = 0: if the convergence delta\n")
cat("  barely separates them, CCM is not the limiting instrument's fault - the\n")
cat("  signal is simply not of the kind CCM detects.\n")

# -----------------------------------------------------------------------------
# 14.5  Recommendation
# -----------------------------------------------------------------------------

cat("\n=== 14.5 DESIGN RECOMMENDATION ===\n")
# Stated in COLONIES, matching fig_depth_reliability. The drug-test version of
# this line is what made the old figure look like current practice was already
# adequate: 1.5 colonies x 3 AG drugs = 4.5 "tests" clears a 5-test target, but
# those 3 drugs on one colony carry ~1.4 tests' worth of independent
# information, not 3. Colonies is also the unit a lab can act on.
cat(sprintf("  DEPTH    : %d-%d colonies per visit for reliability 0.8 (now: mean %.1f).\n",
            min(col_for_080$colonies_needed), max(col_for_080$colonies_needed), cur_colonies))
cat(sprintf("             That is ~%.1fx current depth. The equivalent figure in drug-tests\n",
            mean(col_for_080$colonies_needed) / cur_colonies))
cat(sprintf("             is %d-%d, but do NOT quote it: it treats the %d aminoglycosides run\n",
            min(rel_tab$m_for_rel_0.8), max(rel_tab$m_for_rel_0.8), k_drugs))
cat(sprintf("             on one colony as %d independent draws when their within-colony\n", k_drugs))
cat(sprintf("             correlation (rho = %.2f) makes them worth %.1f.\n",
            ag_within_rho, m_eff_per_colony))
cat(sprintf("             %.0f%% of carbapenem and FQ visits currently have a single test,\n",
            max(rel_tab$pct_single_test)))
cat("             where signal and noise cannot be separated at all.\n")
cat(sprintf("  INTERVAL : %.0f-%.0f days. Shorter raises retained AC toward the %.2f\n",
            14, 28, a_hat))
cat("             ceiling; it cannot exceed it\n")
# Equivalent designs rather than one "best" point: several (L, n) combinations
# buy the same minimum detectable effect, and they differ enormously in what
# they cost and how long they take. Show the frontier and let the reader pick.
cat("  LENGTH x COHORT (at depth m = 10, interval 28 d):\n")
cat("    NOTE the lever comparison above: visits-per-patient (+0.46) and\n")
cat("    tests-per-visit (+0.42) both beat cohort size (+0.32) at this\n")
cat("    operating point. A registry that scales N while keeping the same\n")
cat("    sparse, clinically-driven visit schedule buys the least of the three.\n")
cat("    That still supports the headline - depth of follow-up beats breadth\n")
cat("    of enrolment - but see the caveat above before extending it to the\n")
cat("    sampling INTERVAL specifically, which this grid cannot settle.\n")
frontier <- mde %>% filter(m == 10, dt == 28) %>%
  mutate(visits_total = n_pat * L, study_years = round(L * 28 / 365, 1)) %>%
  group_by(min_detectable_beta) %>%
  arrange(visits_total, .by_group = TRUE) %>% slice_head(n = 2) %>% ungroup() %>%
  arrange(min_detectable_beta, visits_total)
print(frontier %>% transmute(`min beta` = min_detectable_beta,
                             patients = n_pat, `visits/patient` = L,
                             `total visits` = visits_total,
                             `study duration (yr)` = study_years))
cat("    Note the duration column: L = 48 visits at a 28-day interval is a\n")
cat("    3.7-year study. Detecting beta = 0.2 needs only about half of that.\n")
write_csv(frontier, "results/design/design_frontier.csv")
cat("  ESTIMATOR: a cross-lagged mixed model is the population-level tool this\n")
cat("             data shape supports. EDM/CCM needs a smooth trajectory the\n")
cat("             biology does not provide at this timescale.\n")

write_csv(mde, "results/design/min_detectable_effect.csv")
cat("\n=== SECTION 14 COMPLETE ===\n")
cat("  results/design/  depth_reliability.csv, colonies_for_reliability.csv,\n")
cat("                   interval_autocorrelation.csv, interval_vs_treatment_cycles.csv,\n")
cat("                   visit_phase_locking.csv, power_grid_dt.csv, lever_comparison.csv,\n")
cat("                   min_detectable_effect.csv, ccm_power.csv, design_frontier.csv\n")
cat("  figures/design/  fig_depth_reliability.pdf, fig_interval_autocorrelation.pdf,\n")
cat("                   fig_power_grid.pdf\n")

# ---- end of AMR_pipeline.R --------------------------------------------------
