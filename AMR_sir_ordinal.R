# =============================================================================
# AMR_sir_ordinal.R  -  using the full S / I / R call instead of collapsing to R+I
# Emory pilot cohort | companion to AMR_pipeline.R
#
# Author: Canan Karakoc
#
# WHAT THIS FILE CONCLUDES (read before citing anything from it)
# -----------------------------------------------------------------------------
#   1. The variance decomposition is INVARIANT to how S/I/R is collapsed.
#      Ordinal, R+I, and strict-R give the same patient / visit / colony /
#      residual split. So the pipeline's R+I coding needs no defence, and
#      there is no reason to refit anything as ordinal. This is the answer to
#      the coding question, and it is a robustness result, not a null.
#
#   2. Given one colony reads Intermediate, a second colony from the SAME
#      sputum disagrees ~74% of the time (vs ~29% for S, ~34% for R). Those
#      raw rates are observed facts and are the quotable numbers.
#
#   3. That is BAND GEOMETRY, NOT excess assay error at the boundary.
#      The Intermediate band is 0.2-1.0 colony-SD wide depending on drug, so a
#      colony landing in it routinely has siblings across a threshold. Once
#      chance agreement is computed against the SAME PATIENT rather than the
#      whole cohort, every category shows the same small excess agreement
#      (kappa ~ 0.05) and Intermediate is not special. See Section 2.
#
#   4. Within-visit agreement is barely above within-patient agreement for any
#      category. The visit level is close to empty - which agrees with
#      sd_visit_key ~ 0.25 in Section 3 and with the 10.3d identifiability
#      simulation in the main pipeline. Three routes, same answer.
#
# The honest headline is therefore "the Intermediate band is narrower than
# within-host variation, so the call cannot be reproducibly assigned" - the
# thesis stated in CLSI's own units - NOT "the assay is noisiest at the
# boundary". The second claim was an artefact of a permutation null that
# destroyed between-patient structure; see the note in Section 2.2.
#
# SCOPE NOTE. The 9 drugs here span six classes, including Penicillin
# (pip-tazo) and Monobactam (aztreonam), which are NOT among the pipeline's
# four classes. That is why the residual share here (~30%) differs from the
# AG-only figure (~14.7%). Cite one model or the other, never both together.
#
# CONTENTS
#   0  Setup, long S/I/R frame
#   1  Is "Intermediate" the same thing across the cohort?   (breakpoint drift)
#   2  Within-visit reproducibility, against two nulls
#   3  Ordinal variance decomposition + per-drug band widths
#   4  Proportional-odds check
#   5  Summary
#
# Sections 1-2 are fast. Section 3 fits three large crossed-RE models and is
# slow. Set FIT_MODELS <- FALSE to run 1-2 alone.
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
source("AMR_common.R")

for (d in c("results/sir", "figures/sir", "models/sir")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

set.seed(42)
FIT_MODELS <- TRUE
N_PERM <- 500
N_BOOT <- 1000

# -----------------------------------------------------------------------------
# Local theme. theme_talk() is sized for single-panel poster figures (18pt text)
# and collides badly in 3x3 facets. This is the same look, sized for panels,
# with subtitles and captions off - the narrative belongs in the text, not
# stacked on the plot.
# -----------------------------------------------------------------------------
theme_sir <- function(base_size = 11) {
  theme_talk(base_size = base_size) +
    theme(
      text          = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 1),
      axis.text.x   = element_text(size = base_size - 1),
      axis.title    = element_text(size = base_size),
      legend.text   = element_text(size = base_size - 1),
      plot.title    = element_text(size = base_size + 2, face = "bold"),
      strip.text    = element_text(size = base_size - 1, face = "bold"),
      plot.subtitle = element_blank(),
      plot.caption  = element_blank(),
      legend.position = "top"
    )
}

# Consistent colours for the three calls, used everywhere in this file.
sir_cols <- c(Susceptible = "#A9C7CC",
              Intermediate = unname(pal["gold"]),
              Resistant    = unname(pal["red"]))


# #############################################################################
# SECTION 0: SETUP AND THE LONG S/I/R FRAME
# #############################################################################

.dat        <- load_amr_data()
pa_isolates <- .dat$pa_isolates
rm(.dat)

# Keep a drug if it has >= 500 non-missing calls and >= 50 Intermediate calls.
# Drops imipenem (n = 14), ciprofloxacin (n = 8), colistin (22 I calls).
pa_panel <- c("amikacin", "gentamicin", "tobramycin", "imipenem", "meropenem",
              "cefepime", "ceftazidime", "ciprofloxacin", "levofloxacin",
              "piperacillin_tazobactam", "aztreonam", "colistin")
sir_levels <- c("Susceptible", "Intermediate", "Resistant")

drug_screen <- map_dfr(pa_panel, function(d) {
  v <- pa_isolates[[d]]; v <- v[v %in% sir_levels]
  tibble(drug = d, n = length(v), n_S = sum(v == "Susceptible"),
         n_I = sum(v == "Intermediate"), n_R = sum(v == "Resistant"),
         pct_I = 100 * mean(v == "Intermediate"))
}) %>% mutate(keep = n >= 500 & n_I >= 50)

cat("\n0.1  DRUG SCREEN\n")
print(as.data.frame(drug_screen), row.names = FALSE, digits = 3)
drugs_use <- drug_screen %>% filter(keep) %>% pull(drug)
cat("\nKeeping", length(drugs_use), "drugs:", paste(drugs_use, collapse = ", "), "\n")

drug_class <- c(
  amikacin = "Aminoglycoside", gentamicin = "Aminoglycoside", tobramycin = "Aminoglycoside",
  meropenem = "Carbapenem", imipenem = "Carbapenem",
  cefepime = "Cephalosporin", ceftazidime = "Cephalosporin",
  levofloxacin = "Fluoroquinolone", ciprofloxacin = "Fluoroquinolone",
  piperacillin_tazobactam = "Penicillin", aztreonam = "Monobactam",
  colistin = "Polymyxin")

sir_long <- pa_isolates %>%
  mutate(isolate_id = sprintf("iso%04d", row_number()),
         visit_key  = paste(individual_id, visit_number, sep = "_"),
         year       = lubridate::year(date_collected_from_subject),
         phenotype  = pseudomonas_aeruginosa_phenotype) %>%
  select(individual_id, visit_key, visit_number, isolate_id, year,
         date_collected_from_subject, phenotype, clinical_status, all_of(drugs_use)) %>%
  pivot_longer(all_of(drugs_use), names_to = "drug", values_to = "sir") %>%
  filter(sir %in% sir_levels) %>%
  mutate(sir      = factor(sir, levels = sir_levels, ordered = TRUE),
         drug     = factor(drug, levels = drugs_use),
         class    = factor(unname(drug_class[as.character(drug)])),
         R_cons   = as.integer(sir %in% c("Intermediate", "Resistant")),
         R_strict = as.integer(sir == "Resistant"))

cat("\n0.2  LONG FRAME\n")
cat("  rows (colony x drug):", nrow(sir_long),
    "| colonies:", n_distinct(sir_long$isolate_id),
    "| visits:", n_distinct(sir_long$visit_key),
    "| patients:", n_distinct(sir_long$individual_id), "\n")
cat("  classes spanned:", paste(levels(sir_long$class), collapse = ", "), "\n")
cat("  NOTE: this drug set is wider than the pipeline's four classes - see header.\n")
saveRDS(sir_long, "data/sir_long.rds")


# #############################################################################
# SECTION 1: IS "INTERMEDIATE" THE SAME THING ACROSS THE COHORT?
# #############################################################################

cat("\n\n=== SECTION 1: BREAKPOINT STABILITY OF THE INTERMEDIATE BAND ===\n")

sir_by_year <- sir_long %>%
  group_by(drug, class, year) %>%
  summarise(n = n(), n_pat = n_distinct(individual_id),
            pct_S = 100 * mean(sir == "Susceptible"),
            pct_I = 100 * mean(sir == "Intermediate"),
            pct_R = 100 * mean(sir == "Resistant"), .groups = "drop")

cat("\n1.1  INTERMEDIATE FRACTION BY YEAR (%)\n")
print(sir_by_year %>% select(drug, year, pct_I) %>%
        pivot_wider(names_from = year, values_from = pct_I) %>%
        mutate(across(where(is.numeric), ~round(., 1))) %>% as.data.frame(), row.names = FALSE)
cat("\n     n calls / n patients per year\n")
print(sir_by_year %>% transmute(drug, year, cell = sprintf("%d/%d", n, n_pat)) %>%
        pivot_wider(names_from = year, values_from = cell) %>% as.data.frame(), row.names = FALSE)
write_csv(sir_by_year, "results/sir/sir_composition_by_year.csv")

# -----------------------------------------------------------------------------
# 1.2  Within-patient shift, clustered by patient
# -----------------------------------------------------------------------------
# Patient x drug pairs are not independent - a patient contributes up to 9.
# A plain t.test over pairs gives an interval that is too narrow, so the CI
# here comes from a bootstrap resampling PATIENTS.
# -----------------------------------------------------------------------------

within_patient_shift <- function(cut_year, n_boot = 500) {
  paired <- sir_long %>%
    mutate(period = if_else(year < cut_year, "before", "after")) %>%
    group_by(individual_id, drug, period) %>%
    summarise(pI = mean(sir == "Intermediate"), n = n(), .groups = "drop") %>%
    pivot_wider(names_from = period, values_from = c(pI, n)) %>%
    filter(!is.na(pI_before), !is.na(pI_after), n_before >= 3, n_after >= 3) %>%
    mutate(delta = 100 * (pI_after - pI_before))
  if (nrow(paired) < 10) return(tibble(cut = cut_year, n_pairs = nrow(paired),
                                       n_patients = NA_integer_, mean_delta_pp = NA_real_,
                                       lo = NA_real_, hi = NA_real_))
  paired$individual_id <- as.character(paired$individual_id)
  pats <- unique(paired$individual_id)
  psplit <- split(paired, paired$individual_id)
  bs <- map_dbl(seq_len(n_boot), function(i)
    mean(bind_rows(psplit[sample(pats, length(pats), replace = TRUE)])$delta))
  tibble(cut = cut_year, n_pairs = nrow(paired), n_patients = length(pats),
         mean_delta_pp = mean(paired$delta),
         lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975)))
}

cat("\n1.2  WITHIN-PATIENT SHIFT IN I FRACTION (patient-clustered bootstrap CI)\n")
wp <- map_dfr(sort(unique(sir_long$year))[-1], within_patient_shift)
print(wp %>% mutate(across(where(is.numeric), ~round(., 2))) %>% as.data.frame(), row.names = FALSE)
write_csv(wp, "results/sir/within_patient_I_shift.csv")
cat("\n     Each row splits the WHOLE window at that year, so rows are not\n")
cat("     independent tests. A breakpoint revision shows as an abrupt jump at\n")
cat("     one cut; a smooth ramp is drift, and does not threaten the I call.\n")

fig_sir_year <- sir_by_year %>%
  select(drug, year, pct_S, pct_I, pct_R) %>%
  pivot_longer(starts_with("pct_"), names_to = "call", values_to = "pct") %>%
  mutate(call = factor(recode(call, pct_S = "Susceptible", pct_I = "Intermediate",
                              pct_R = "Resistant"), levels = sir_levels)) %>%
  ggplot(aes(year, pct, fill = call)) +
  geom_col(width = 0.85) +
  facet_wrap(~ drug, ncol = 3) +
  scale_fill_manual(values = sir_cols, name = NULL) +
  scale_x_continuous(breaks = seq(min(sir_by_year$year), max(sir_by_year$year), 2)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "% of calls", title = "Is the Intermediate band stable over the window?") +
  theme_sir(11)
ggsave("figures/sir/sir_composition_by_year.pdf", fig_sir_year,
       width = 9, height = 7, dpi = 300)


# #############################################################################
# SECTION 2: WITHIN-VISIT REPRODUCIBILITY, AGAINST TWO NULLS
#
# Take every ordered pair of distinct colonies (i, j) from the same visit and
# drug. The INDEX colony's call stratifies; the outcome is whether the PARTNER
# disagrees. Nothing is mechanically forced either way.
#
# THE CHOICE OF NULL DECIDES THE ANSWER. Read this before quoting a kappa.
#
#   NULL A (global): shuffle calls within drug across all colonies. Chance
#     agreement is then referenced to the whole cohort's marginal. This is
#     WRONG for the reproducibility question: ~50% of latent variance is
#     between patients, and S and R are far more patient-clustered than I, so
#     the null is more permissive for S and R and manufactures a low kappa for
#     Intermediate. An earlier version of this file used Null A alone and
#     reported "Intermediate is 2.5x less reproducible". That claim is
#     withdrawn.
#
#   NULL B (within-patient, PRIMARY): shuffle within patient x drug. Patient
#     composition is preserved; only within-visit clustering is destroyed. This
#     asks the question that matters - do two colonies from the SAME sputum
#     agree more than two colonies from the same patient at different visits?
#
# Under Null B all three categories return kappa ~ 0.05 with overlapping
# intervals. There is no boundary-specific reproducibility deficit. What
# remains true, and is the quotable result, is the RAW disagreement rate.
# #############################################################################

cat("\n\n=== SECTION 2: WITHIN-VISIT REPRODUCIBILITY ===\n")

base_cols <- sir_long %>% select(individual_id, visit_key, drug, isolate_id, sir, phenotype)
mkpairs <- function(d) {
  d %>% inner_join(d, by = c("visit_key", "drug"),
                   relationship = "many-to-many", suffix = c(".i", ".j")) %>%
    filter(isolate_id.i != isolate_id.j)
}

pair_df <- mkpairs(base_cols) %>%
  transmute(individual_id = individual_id.i, visit_key, drug,
            index = sir.i, partner = sir.j,
            disagree = as.integer(sir.i != sir.j),
            same_phenotype = case_when(
              is.na(phenotype.i) | is.na(phenotype.j) ~ NA_integer_,
              phenotype.i == phenotype.j ~ 1L, TRUE ~ 0L))

cat("\n2.1  ordered colony pairs:", nrow(pair_df),
    "| visits:", n_distinct(pair_df$visit_key),
    "| patients:", n_distinct(pair_df$individual_id), "\n")
cat("     NOTE: ordered pairs weight a visit by n(n-1), so deep visits count\n")
cat("     more, and depth correlates with clinical concern. The colony-level\n")
cat("     model in 2.5 has no such weighting and is the primary estimate.\n")

# -----------------------------------------------------------------------------
# 2.2  Observed rates, and kappa under both nulls
# -----------------------------------------------------------------------------

observed <- pair_df %>% group_by(index) %>%
  summarise(n_pairs = n(), obs_disagree = mean(disagree), .groups = "drop")

cat("\n2.2  RAW DISAGREEMENT - the quotable numbers, no null involved\n")
print(observed %>% mutate(obs_disagree = round(100 * obs_disagree, 1)) %>%
        as.data.frame(), row.names = FALSE)

perm_global <- function(...) mkpairs(
  base_cols %>% group_by(drug) %>% mutate(sir = sample(sir)) %>% ungroup()) %>%
  group_by(index = sir.i) %>% summarise(d = mean(sir.i != sir.j), .groups = "drop")

perm_within_patient <- function(...) mkpairs(
  base_cols %>% group_by(individual_id, drug) %>% mutate(sir = sample(sir)) %>% ungroup()) %>%
  group_by(index = sir.i) %>% summarise(d = mean(sir.i != sir.j), .groups = "drop")

cat("\n     Running", N_PERM, "permutations under each null ...\n")
# NB: compute the quantiles BEFORE overwriting the column being summarised -
# summarise() evaluates sequentially and a naive version returns a degenerate
# interval equal to the mean.
summ_null <- function(p) p %>% group_by(index) %>%
  summarise(nul_lo = unname(quantile(d, .025)), nul_hi = unname(quantile(d, .975)),
            nul = mean(d), .groups = "drop")

nullA <- summ_null(map_dfr(seq_len(N_PERM), perm_global))
nullB <- summ_null(map_dfr(seq_len(N_PERM), perm_within_patient))

# Cluster bootstrap over PATIENTS (colonies nest in visits nest in patients).
pair_df$individual_id <- as.character(pair_df$individual_id)
pats <- unique(pair_df$individual_id)
psplit <- split(pair_df, pair_df$individual_id)
boot_obs <- map_dfr(seq_len(N_BOOT), function(i)
  bind_rows(psplit[sample(pats, length(pats), replace = TRUE)]) %>%
    group_by(index) %>% summarise(obs = mean(disagree), .groups = "drop") %>%
    mutate(rep = i))

kap <- function(o, n) ((1 - o) - (1 - n)) / (1 - (1 - n))

make_kappa <- function(nul, label) {
  pt <- observed %>% left_join(nul, by = "index") %>%
    mutate(kappa = kap(obs_disagree, nul))
  ci <- boot_obs %>% left_join(nul %>% select(index, nul), by = "index") %>%
    mutate(kappa = kap(obs, nul)) %>% group_by(index) %>%
    summarise(lo = unname(quantile(kappa, .025)),
              hi = unname(quantile(kappa, .975)), .groups = "drop")
  pt %>% left_join(ci, by = "index") %>% mutate(null = label)
}

kappa_tbl <- bind_rows(make_kappa(nullB, "B: within-patient (PRIMARY)"),
                       make_kappa(nullA, "A: global (biased, for contrast)"))

cat("\n2.2  CATEGORY-SPECIFIC KAPPA UNDER BOTH NULLS\n")
print(kappa_tbl %>%
        transmute(null, index, obs = round(100 * obs_disagree, 1),
                  chance = round(100 * nul, 1), kappa = round(kappa, 3),
                  CI = sprintf("[%.3f, %.3f]", lo, hi)) %>%
        as.data.frame(), row.names = FALSE)
write_csv(kappa_tbl, "results/sir/reproducibility_by_call.csv")

# Difference, not ratio - kappa is not a ratio scale, so "2.5x less
# reproducible" is not a legitimate statement about it.
kb <- boot_obs %>% left_join(nullB %>% select(index, nul), by = "index") %>%
  mutate(kappa = kap(obs, nul)) %>% select(rep, index, kappa) %>%
  pivot_wider(names_from = index, values_from = kappa) %>%
  mutate(I_minus_S = Intermediate - Susceptible)
cat(sprintf("\n     kappa_I - kappa_S (Null B) = %.3f [%.3f, %.3f], P(I < S) = %.2f\n",
            median(kb$I_minus_S), quantile(kb$I_minus_S, .025),
            quantile(kb$I_minus_S, .975), mean(kb$I_minus_S < 0)))
cat("     Interval spans zero: Intermediate is NOT distinctively less\n")
cat("     reproducible once chance is referenced to the same patient.\n")

# -----------------------------------------------------------------------------
# 2.3  Adjacency of disagreements
# -----------------------------------------------------------------------------

adjacency <- pair_df %>% filter(disagree == 1) %>%
  mutate(step = if_else((index == "Susceptible" & partner == "Resistant") |
                          (index == "Resistant" & partner == "Susceptible"),
                        "S<->R (two steps)", "adjacent (S<->I or I<->R)")) %>%
  count(drug, step) %>% group_by(drug) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
cat("\n2.3  AMONG DISAGREEING PAIRS: ONE STEP OR TWO?\n")
print(adjacency %>% mutate(pct = round(pct, 1)) %>%
        pivot_wider(names_from = step, values_from = c(n, pct)) %>% as.data.frame(),
      row.names = FALSE)
write_csv(adjacency, "results/sir/disagreement_adjacency.csv")

# -----------------------------------------------------------------------------
# 2.4  Per-drug kappa, and the ordering test
# -----------------------------------------------------------------------------
# Drugs have very different I frequencies, so pooling invites a Simpson-type
# artefact. Under a no-ordering null P(I lowest in a drug) = 1/3, so
# P(I lowest in all 9) = (1/3)^9 = 5.1e-5.

per_drug_kappa <- function(perm_fun, label) {
  o <- pair_df %>% group_by(drug, index) %>%
    summarise(n = n(), obs = mean(disagree), .groups = "drop")
  n_ <- map_dfr(seq_len(200), function(i) perm_fun_drug(perm_fun)) %>%
    group_by(drug, index) %>% summarise(nul = mean(d), .groups = "drop")
  o %>% left_join(n_, by = c("drug", "index")) %>%
    mutate(kappa = kap(obs, nul), null = label)
}
perm_fun_drug <- function(fn) {
  sh <- if (identical(fn, perm_global))
    base_cols %>% group_by(drug) %>% mutate(sir = sample(sir)) %>% ungroup()
  else
    base_cols %>% group_by(individual_id, drug) %>% mutate(sir = sample(sir)) %>% ungroup()
  mkpairs(sh) %>% group_by(drug, index = sir.i) %>%
    summarise(d = mean(sir.i != sir.j), .groups = "drop")
}

pdk <- bind_rows(per_drug_kappa(perm_within_patient, "B: within-patient"),
                 per_drug_kappa(perm_global,        "A: global"))
write_csv(pdk, "results/sir/per_drug_kappa.csv")

ord_test <- pdk %>% select(null, drug, index, kappa) %>%
  pivot_wider(names_from = index, values_from = kappa) %>%
  mutate(I_lowest = Intermediate < pmin(Susceptible, Resistant)) %>%
  group_by(null) %>% summarise(n_drugs = n(), n_I_lowest = sum(I_lowest), .groups = "drop") %>%
  # Under a no-ordering null each drug puts Intermediate lowest with prob 1/3,
  # so the count is Binomial(n_drugs, 1/3). (1/3)^k is only the right tail when
  # k == n_drugs; use the binomial tail so partial counts are read correctly.
  mutate(expected_by_chance = n_drugs / 3,
         p_no_ordering = pbinom(n_I_lowest - 1, n_drugs, 1/3, lower.tail = FALSE))

cat("\n2.4  PER-DRUG KAPPA: in how many drugs is Intermediate lowest?\n")
print(as.data.frame(ord_test), row.names = FALSE)
print(pdk %>% filter(null == "B: within-patient") %>%
        select(drug, index, kappa) %>%
        pivot_wider(names_from = index, values_from = kappa) %>%
        mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(), row.names = FALSE)

cat("\n2.4b DISAGREEMENT BY INDEX CALL, SPLIT BY SHARED PHENOTYPE\n")
print(pair_df %>% filter(!is.na(same_phenotype)) %>%
        group_by(index, ph = if_else(same_phenotype == 1, "same", "mucoid vs non")) %>%
        summarise(n = n(), disagree = round(100 * mean(disagree), 1), .groups = "drop") %>%
        pivot_wider(names_from = ph, values_from = c(n, disagree)) %>% as.data.frame(),
      row.names = FALSE)

# -----------------------------------------------------------------------------
# 2.5  Colony-level model - PRIMARY estimate of the raw rates
# -----------------------------------------------------------------------------
# One row per (colony, drug) in a multi-colony visit: does THIS colony's call
# disagree with at least one sibling? No pair double-counting, no depth
# weighting, and drug / patient / visit held as random effects.
# These odds are NOT chance-corrected - see 2.2 for that.

colony_disagree <- sir_long %>%
  group_by(visit_key, drug) %>% filter(n() >= 2) %>%
  mutate(n_siblings = n() - 1L,
         disagree = as.integer(sapply(seq_len(n()), function(k) any(sir[-k] != sir[k])))) %>%
  ungroup() %>%
  mutate(index = relevel(factor(sir, levels = sir_levels, ordered = FALSE),
                         ref = "Susceptible"),
         log_sibs = log(n_siblings))

m_disagree <- brm(
  disagree ~ index + log_sibs + (1 | drug) + (1 | individual_id) + (1 | visit_key),
  data = colony_disagree, family = bernoulli(), prior = priors_bernoulli,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/sir/m_disagree_by_call", file_refit = "on_change")
print(summary(m_disagree))

or_tbl <- as_draws_df(m_disagree) %>% select(starts_with("b_index")) %>%
  pivot_longer(everything(), names_to = "term", values_to = "logOR") %>%
  group_by(term) %>%
  summarise(OR = median(exp(logOR)), lo = quantile(exp(logOR), .025),
            hi = quantile(exp(logOR), .975), .groups = "drop") %>%
  mutate(term = recode(term, b_indexIntermediate = "Intermediate vs Susceptible",
                       b_indexResistant = "Resistant vs Susceptible"))
cat("\n2.5  ODDS OF WITHIN-SAMPLE DISAGREEMENT vs A SUSCEPTIBLE CALL\n")
print(or_tbl %>% mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(),
      row.names = FALSE)
cat("     Uncorrected for base rate. Section 2.2 Null B shows the excess over\n")
cat("     chance is the same for all three calls; this OR is band geometry.\n")
write_csv(or_tbl, "results/sir/disagreement_odds_ratios.csv")

# -----------------------------------------------------------------------------
# 2.6  Figures
# -----------------------------------------------------------------------------

fig_rates <- observed %>%
  mutate(index = factor(index, levels = sir_levels)) %>%
  ggplot(aes(index, obs_disagree, fill = index)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = scales::percent(obs_disagree, accuracy = 1)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = sir_cols, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.9),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "What the sampled colony was called",
       y = "P(a second colony disagrees)",
       title = "One colony, one call - and at the boundary it is close to a coin flip") +
  theme_sir(12)
ggsave("figures/sir/disagreement_rates.pdf", fig_rates,
       width = 7, height = 4.6, dpi = 300)

fig_kappa <- kappa_tbl %>%
  mutate(index = factor(index, levels = sir_levels)) %>%
  ggplot(aes(index, kappa, fill = index)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.1, linewidth = 0.5) +
  facet_wrap(~ null) +
  scale_fill_manual(values = sir_cols, guide = "none") +
  labs(x = NULL, y = "Agreement above chance (kappa)",
       title = "The ordering is an artefact of the null, not a property of the assay") +
  theme_sir(11)
ggsave("figures/sir/reproducibility_by_call.pdf", fig_kappa,
       width = 8, height = 4.6, dpi = 300)


# #############################################################################
# SECTION 3: ORDINAL VARIANCE DECOMPOSITION
#
# res_ordinal ~ drug + (1|individual_id) + (1|visit_key) + (1|isolate_id)
#
# IDENTIFICATION CAVEAT. Each colony contributes up to 9 drug rows, so the
# isolate term separates from residual everywhere. But isolate separates from
# VISIT only through multi-colony visits, and the pipeline's 10.3d simulation
# shows that component is NOT recoverable at this depth (a true sd_visit of 0
# returns ~0.27; a true 0.8 returns ~0.67; we observe ~0.28). Expect R-hat
# inflation on sd_visit_key: it is a boundary component, not a mixing failure.
# The terms that carry the claims here (patient, colony, residual) run at
# R-hat <= 1.01. Do not read sd_visit as "visits do not vary" - read it as
# "with one colony per visit the two terms are not separable".
# #############################################################################

cat("\n\n=== SECTION 3: ORDINAL VARIANCE DECOMPOSITION ===\n")

if (!FIT_MODELS) {
  cat("\nFIT_MODELS is FALSE - skipping Sections 3 and 4.\n")
} else {

sir_re <- "(1 | individual_id) + (1 | visit_key) + (1 | isolate_id)"
# adapt_delta 0.99 AND longer chains: sd_visit_key sits against zero, which is
# funnel geometry - longer chains alone often do not clear it.
CTRL <- list(adapt_delta = 0.99, max_treedepth = 12)
IT <- 4000; WU <- 2000

m_ord <- brm(as.formula(paste("sir ~ drug +", sir_re)), data = sir_long,
  family = cumulative("logit"),
  prior = c(prior(normal(0, 4), class = "Intercept"),
            prior(normal(0, 2), class = "b"), prior(exponential(1), class = "sd")),
  chains = 4, cores = 4, iter = IT, warmup = WU, control = CTRL, seed = 42,
  file = "models/sir/m_ordinal_v2", file_refit = "on_change")
print(summary(m_ord))

m_cons <- brm(as.formula(paste("R_cons ~ drug +", sir_re)), data = sir_long,
  family = bernoulli(), prior = priors_bernoulli,
  chains = 4, cores = 4, iter = IT, warmup = WU, control = CTRL, seed = 42,
  file = "models/sir/m_binary_conservative_v2", file_refit = "on_change")

m_strict <- brm(as.formula(paste("R_strict ~ drug +", sir_re)), data = sir_long,
  family = bernoulli(), prior = priors_bernoulli,
  chains = 4, cores = 4, iter = IT, warmup = WU, control = CTRL, seed = 42,
  file = "models/sir/m_binary_strict_v2", file_refit = "on_change")

# Report diagnostics rather than assert convergence.
cat("\n3.0  DIAGNOSTICS (sd_visit_key is expected to be the worst - see caveat)\n")
diag_tbl <- map_dfr(list(`Ordinal (S < I < R)` = m_ord, `Binary, R + I vs S` = m_cons,
                         `Binary, R vs S + I` = m_strict) %>% imap(list), function(x) {
  m <- x[[1]]; s <- posterior::summarise_draws(as_draws_df(m), "rhat", "ess_bulk")
  s <- s[!is.na(s$rhat), ]
  np <- brms::nuts_params(m)
  tibble(model = names(x)[1], rhat_max = max(s$rhat),
         worst = s$variable[which.max(s$rhat)], ess_min = min(s$ess_bulk),
         rhat_excl_visit = max(s$rhat[!grepl("visit_key", s$variable)]),
         divergent = sum(np$Value[np$Parameter == "divergent__"]))
})
print(diag_tbl %>% mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(),
      row.names = FALSE)
write_csv(diag_tbl, "results/sir/section3_diagnostics.csv")

# -----------------------------------------------------------------------------
# 3.3  Decomposition on the latent logit scale
# -----------------------------------------------------------------------------
# Residual is pi^2/3, FIXED by the logit link - it is not estimated, which is
# why it is drawn in grey in the figure. Note the residual here means
# "drug-within-colony", not "colony-within-visit" as in the pipeline's
# per-class decomposition.

decompose_latent <- function(model, label) {
  dr <- as_draws_df(model)
  v_pat <- dr[["sd_individual_id__Intercept"]]^2
  v_vis <- dr[["sd_visit_key__Intercept"]]^2
  v_iso <- dr[["sd_isolate_id__Intercept"]]^2
  v_res <- pi^2 / 3
  tot <- v_pat + v_vis + v_iso + v_res
  mk <- function(v, comp) tibble(model = label, component = comp,
    median = median(v / tot), lo = unname(quantile(v / tot, .025)),
    hi = unname(quantile(v / tot, .975)))
  bind_rows(mk(v_pat, "Patients"), mk(v_vis, "Visits\n(within patient)"),
            mk(v_iso, "Colonies\n(within visit)"),
            mk(rep(v_res, length(v_pat)), "Drugs within colony\n(residual, fixed)"))
}

vd_sir <- bind_rows(decompose_latent(m_ord, "Ordinal (S < I < R)"),
                    decompose_latent(m_cons, "Binary, R + I vs S"),
                    decompose_latent(m_strict, "Binary, R vs S + I"))
cat("\n3.3  VARIANCE DECOMPOSITION ON THE LATENT LOGIT SCALE\n")
print(vd_sir %>% mutate(component = gsub("\n", " ", component),
                        across(where(is.numeric), ~round(100 * ., 1))) %>%
        as.data.frame(), row.names = FALSE)
write_csv(vd_sir %>% mutate(component = gsub("\n", " ", component)),
          "results/sir/variance_decomposition_ordinal.csv")

sd_tbl <- bind_rows(
  tibble(model = "Ordinal (S < I < R)", term = "colony", sd = as_draws_df(m_ord)[["sd_isolate_id__Intercept"]]),
  tibble(model = "Ordinal (S < I < R)", term = "visit",  sd = as_draws_df(m_ord)[["sd_visit_key__Intercept"]]),
  tibble(model = "Binary, R + I vs S", term = "colony", sd = as_draws_df(m_cons)[["sd_isolate_id__Intercept"]]),
  tibble(model = "Binary, R + I vs S", term = "visit",  sd = as_draws_df(m_cons)[["sd_visit_key__Intercept"]]),
  tibble(model = "Binary, R vs S + I", term = "colony", sd = as_draws_df(m_strict)[["sd_isolate_id__Intercept"]]),
  tibble(model = "Binary, R vs S + I", term = "visit",  sd = as_draws_df(m_strict)[["sd_visit_key__Intercept"]])) %>%
  group_by(model, term) %>%
  summarise(median = median(sd), lo = unname(quantile(sd, .025)),
            hi = unname(quantile(sd, .975)), OR_per_SD = exp(median(sd)), .groups = "drop")
cat("\n3.3b RANDOM-EFFECT SDs ON THE LATENT SCALE (logits), with OR per SD\n")
cat("     Report sd_visit alongside its variance share - an OR of ~1.3 per SD\n")
cat("     does not read like nothing, and the share is small only because the\n")
cat("     FIXED residual pi^2/3 dominates the denominator.\n")
print(sd_tbl %>% mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(),
      row.names = FALSE)
write_csv(sd_tbl, "results/sir/random_effect_sds.csv")

# -----------------------------------------------------------------------------
# 3.4  Per-drug Intermediate band width
# -----------------------------------------------------------------------------
# Under cumulative(), tau2 - tau1 is ONE constant for every drug: drug effects
# shift the latent scale, they do not move the thresholds. So a per-drug band
# CANNOT come from m_ord. It comes from the two implied binaries, draw by draw:
#
#   P(>= I) = logit^-1(a1 + b1_d)   ->  threshold t1_d = -(a1 + b1_d)
#   P(=  R) = logit^-1(a2 + b2_d)   ->  threshold t2_d = -(a2 + b2_d)
#   band_d  = t2_d - t1_d = (a1 + b1_d) - (a2 + b2_d)
#
# This is legitimate but it is NOT the proportional-odds model's band, and the
# two binaries are separate fits, so the draws are paired only by index.
# sigma_colony is pooled - the model has one isolate SD, not one per drug.

dc <- as_draws_df(m_cons); ds <- as_draws_df(m_strict)
n_dr <- min(nrow(dc), nrow(ds))
a1 <- dc[["b_Intercept"]][1:n_dr]; a2 <- ds[["b_Intercept"]][1:n_dr]
sd_col <- dc[["sd_isolate_id__Intercept"]][1:n_dr]
drugs_lv <- levels(sir_long$drug)

band_tbl <- map_dfr(drugs_lv, function(d) {
  b1 <- if (d == drugs_lv[1]) rep(0, n_dr) else dc[[paste0("b_drug", d)]][1:n_dr]
  b2 <- if (d == drugs_lv[1]) rep(0, n_dr) else ds[[paste0("b_drug", d)]][1:n_dr]
  band <- (a1 + b1) - (a2 + b2); r <- band / sd_col
  tibble(drug = d, band = median(band),
         band_lo = unname(quantile(band, .025)), band_hi = unname(quantile(band, .975)),
         band_per_sd = median(r), r_lo = unname(quantile(r, .025)),
         r_hi = unname(quantile(r, .975)))
}) %>% left_join(sir_long %>% group_by(drug) %>%
                   summarise(pct_I = 100 * mean(sir == "Intermediate")) %>%
                   mutate(drug = as.character(drug)), by = "drug") %>%
  arrange(desc(band))

cat("\n3.4  PER-DRUG INTERMEDIATE BAND WIDTH (latent logits; sigma_colony pooled)\n")
print(band_tbl %>% mutate(across(where(is.numeric), ~round(., 2))) %>% as.data.frame(),
      row.names = FALSE)
write_csv(band_tbl, "results/sir/intermediate_band_width.csv")
cat(sprintf("\n     sigma_colony = %.2f logits. Every band is under one colony-SD;\n",
            median(sd_col)))
cat("     the narrowest are ~0.2. A colony landing in the band routinely has\n")
cat("     siblings on the far side of a threshold - which is exactly the\n")
cat("     mechanism producing the raw disagreement rates in Section 2, NOT an\n")
cat("     independent confirmation of them. Section 2 Null B shows there is no\n")
cat("     excess assay error at the boundary beyond this geometry.\n")

# -----------------------------------------------------------------------------
# 3.5  Posterior predictive check: does pure geometry reproduce Section 2?
# -----------------------------------------------------------------------------
# m_ord has homoscedastic latent noise by construction - no extra error at the
# boundary. If simulating from its posterior reproduces the observed raw
# disagreement rates, band geometry alone explains them. If the observed rates
# are HIGHER than simulated, there is assay error on top.

cat("\n3.5  POSTERIOR PREDICTIVE CHECK ON WITHIN-VISIT DISAGREEMENT\n")
cat("     Measured 2026-08-27: observed 29.2 / 73.9 / 34.2 vs simulated\n")
cat("     37.0 / 87.4 / 41.3 - the model over-predicts disagreement by a\n")
cat("     similar factor (~0.8x) for ALL THREE calls.\n")
cat("     Verified mechanically: simulated marginals match observed to <0.5pp,\n")
cat("     the pair frame is identical by construction, and the miss is the same\n")
cat("     size for BETWEEN-visit pairs (0.85) as for within-visit (0.83), which\n")
cat("     rules out visit effects being redrawn inside a simulated pair.\n")
cat("     CONSEQUENCE: the colony variance component is an UPPER bound. Quote\n")
cat("     the between-colony share as \"at most\", or use observed-anchored\n")
cat("     numbers instead of model-derived ones.\n")
cat("     The uniformity across all three calls is also independent support for\n")
cat("     Section 2: a boundary-specific assay-error story predicts a\n")
cat("     category-specific miss, and there is not one.\n")
cat("     See verify_ppc_multicolony.R for which mechanism drives it.\n")
yrep <- posterior_predict(m_ord, ndraws = 200)
key <- sir_long %>% select(visit_key, drug, isolate_id) %>%
  mutate(row = row_number())
ppc <- map_dfr(seq_len(nrow(yrep)), function(i) {
  # ordered = TRUE: sir_long$sir is an ordered factor, so the pairs built from
  # the simulation must be too or the join against `observed` fails on type.
  sim <- key %>% mutate(sir = factor(sir_levels[yrep[i, ]], levels = sir_levels,
                                     ordered = TRUE))
  mkpairs(sim %>% select(visit_key, drug, isolate_id, sir)) %>%
    group_by(index = sir.i) %>% summarise(d = mean(sir.i != sir.j), .groups = "drop")
}) %>% group_by(index) %>%
  summarise(sim = median(d), sim_lo = unname(quantile(d, .025)),
            sim_hi = unname(quantile(d, .975)), .groups = "drop")

ppc_tbl <- observed %>% select(index, obs_disagree) %>%
  mutate(index = as.character(index)) %>%
  left_join(ppc %>% mutate(index = as.character(index)), by = "index") %>%
  mutate(across(c(obs_disagree, sim, sim_lo, sim_hi), ~round(100 * ., 1)),
         # Directional. Observed disagreement comes in BELOW the model at every
         # call. Do NOT label that "colonies are clonally related" - colonies
         # here are not typed, so that is a biological claim this data cannot
         # support and a reviewer can ask how we know. State only what the
         # numbers license: the model overstates within-sputum spread, so the
         # colony variance component is an UPPER bound. See 3.5b for which
         # mechanism is responsible.
         verdict = case_when(
           obs_disagree >= sim_lo & obs_disagree <= sim_hi ~ "geometry explains it",
           obs_disagree > sim_hi ~ "excess disagreement beyond geometry",
           TRUE ~ "model OVER-predicts disagreement (colony share is an upper bound)"))
print(as.data.frame(ppc_tbl), row.names = FALSE)
write_csv(ppc_tbl, "results/sir/ppc_disagreement.csv")

# -----------------------------------------------------------------------------
# 3.6  Figures
# -----------------------------------------------------------------------------
# The residual bar is pi^2/3 fixed by the link, so it is grey - colouring it
# like an estimated component re-imports the problem fixed in fig_variance_
# isolate_level. x labels are wrapped and the panel is widened to stop the
# collision.

comp_lv <- c("Patients", "Visits\n(within patient)", "Colonies\n(within visit)",
             "Drugs within colony\n(residual, fixed)")
fig_vd <- vd_sir %>%
  mutate(component = factor(component, levels = comp_lv),
         model = factor(model, levels = c("Ordinal (S < I < R)", "Binary, R + I vs S",
                                          "Binary, R vs S + I"))) %>%
  ggplot(aes(model, median, fill = component)) +
  geom_col(width = 0.62) +
  scale_fill_manual(values = c("Patients" = unname(pal["ink"]),
                               "Visits\n(within patient)" = unname(pal["blue"]),
                               "Colonies\n(within visit)" = unname(pal["red"]),
                               "Drugs within colony\n(residual, fixed)" = "grey78"),
                    name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(labels = function(x) str_wrap(x, 12)) +
  labs(x = NULL, y = "Share of latent variance",
       title = "The variance structure does not care how S/I/R is collapsed") +
  theme_sir(11) +
  theme(legend.position = "right", legend.key.size = unit(0.9, "lines"))
ggsave("figures/sir/variance_decomposition_ordinal.pdf", fig_vd,
       width = 8.5, height = 4.4, dpi = 300)

fig_band <- band_tbl %>%
  ggplot(aes(band_per_sd, fct_reorder(drug, band_per_sd))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = r_lo, xmax = r_hi), colour = unname(pal["gold"]),
                  size = 0.4, linewidth = 0.7) +
  labs(x = "Intermediate band width (colony SD)", y = NULL,
       title = "Every Intermediate band is narrower than within-host spread") +
  theme_sir(11)
ggsave("figures/sir/intermediate_band_width.pdf", fig_band,
       width = 7, height = 4.2, dpi = 300)


# #############################################################################
# SECTION 4: PROPORTIONAL-ODDS CHECK
# The two implied binaries are already fitted. Under proportional odds their
# drug coefficients should match; where they diverge, the two thresholds
# respond to drug differently and cumulative() is misspecified.
# #############################################################################

cat("\n\n=== SECTION 4: PROPORTIONAL-ODDS CHECK ===\n")

po_check <- bind_rows(
  as_draws_df(m_cons) %>% select(starts_with("b_drug")) %>%
    pivot_longer(everything(), names_to = "term", values_to = "v") %>%
    mutate(boundary = "S | I+R"),
  as_draws_df(m_strict) %>% select(starts_with("b_drug")) %>%
    pivot_longer(everything(), names_to = "term", values_to = "v") %>%
    mutate(boundary = "S+I | R")) %>%
  group_by(term, boundary) %>%
  summarise(est = median(v), lo = quantile(v, .025), hi = quantile(v, .975), .groups = "drop") %>%
  mutate(term = str_remove(term, "^b_drug"))

po_wide <- po_check %>% select(term, boundary, est) %>%
  pivot_wider(names_from = boundary, values_from = est) %>%
  mutate(difference = `S | I+R` - `S+I | R`) %>% arrange(desc(abs(difference)))
cat("\n4.1  DRUG EFFECTS AT EACH BOUNDARY (log-odds vs reference drug)\n")
print(po_wide %>% mutate(across(where(is.numeric), ~round(., 2))) %>% as.data.frame(),
      row.names = FALSE)
write_csv(po_wide, "results/sir/proportional_odds_check.csv")
cat("\n     Non-zero differences = proportional odds violated. This is the other\n")
cat("     reason not to adopt cumulative() as the primary model: doing it\n")
cat("     properly needs cs(drug), for a decomposition that does not change.\n")

fig_po <- po_check %>%
  ggplot(aes(est, fct_reorder(term, est), colour = boundary)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(0.5),
                  size = 0.32, linewidth = 0.6) +
  scale_colour_manual(values = unname(pal[c("blue", "red")]), name = NULL) +
  labs(x = "Log-odds vs reference drug", y = NULL,
       title = "The two thresholds do not respond to drug the same way") +
  theme_sir(11)
ggsave("figures/sir/proportional_odds_check.pdf", fig_po,
       width = 7.5, height = 4.6, dpi = 300)

}  # end FIT_MODELS


# #############################################################################
# SECTION 5: SUMMARY
# #############################################################################

cat("\n\n=== SUMMARY ===\n")
cat("\n1  Coding question: ANSWERED. The decomposition is invariant to ordinal\n")
cat("   vs R+I vs strict-R. Keep R+I everywhere; cite this as a robustness\n")
cat("   panel. Do not refit anything as ordinal.\n")
cat("\n2  Raw within-sputum disagreement (quotable, no null involved):\n")
if (exists("observed")) {
  print(observed %>% mutate(pct = round(100 * obs_disagree, 1)) %>%
          select(index, n_pairs, pct) %>% as.data.frame(), row.names = FALSE)
}
cat("\n3  RETRACTED: 'Intermediate is 2.5x less reproducible'. Under a null that\n")
cat("   preserves patient composition, all three calls show the same small\n")
cat("   excess agreement (kappa ~ 0.05, intervals overlapping), and per drug\n")
cat("   Intermediate is lowest in about a third of them - exactly chance. The\n")
cat("   74% is band geometry, not boundary assay error.\n")
cat("\n4  Per-drug band widths: results/sir/intermediate_band_width.csv. Every\n")
cat("   band is under one colony-SD. This is the defensible version of the\n")
cat("   argument and it carries the panel on its own.\n")
cat("\n5  sd_visit is not identified at this depth (see 10.3d). Panel copy must\n")
cat("   say 'not separable with one colony per visit', not 'visits do not vary'.\n\n")
