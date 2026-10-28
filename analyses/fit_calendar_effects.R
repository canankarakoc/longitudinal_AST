# =============================================================================
# analyses/fit_calendar_effects.R
#
# Calendar effects that the longitudinal claims have to survive:
#
#   1  SAMPLING DESIGN over calendar time: patients, cultures, colonies per
#      culture, single-colony share, interval, clinical status, sample type,
#      mucoid share, untested isolates - by half-year.
#   2  BREAKPOINTS. In 2017-2022 the only CLSI revision touching our agents for
#      P. aeruginosa is fluoroquinolones: ciprofloxacin and levofloxacin
#      lowered one dilution in M100 29th ed. (January 2019). Aminoglycoside and
#      piperacillin-tazobactam revisions for P. aeruginosa came in 2023, after
#      the window. Carbapenem, cephalosporin and aztreonam breakpoints did not
#      change in the window. When the laboratory adopted the 2019 FQ values is
#      NOT in the data (no MICs, no method field) and must come from the lab;
#      commercial panels lagged FDA recognition. So the adoption date is
#      estimated: for each drug, a colony-level logistic model with a patient
#      random intercept and a step at month c, scanned over c. A breakpoint
#      change predicts a step in levofloxacin only (others are negative
#      controls), WITHIN patient, with S falling.
#   3  FLUOROQUINOLONE TREND: raw cohort rise vs within-patient change; linear
#      calendar trend vs step vs both, adjusted for clinical status and
#      morphotype. 77 of 92 patients enter in 2017, so calendar time and
#      follow-up time are nearly collinear: trend and step are told apart by
#      shape, not by staggered entry.
#   4  SENSITIVITY of earlier analyses to period: within-culture discordance
#      (analysis 2, colony frame) before 2019, 2019, 2020 on; and the
#      within-patient FQ variance (analysis 1) with and without the step.
#
# Run from the repo root:  Rscript analyses/fit_calendar_effects.R
# Output: results/calendar/   figures/calendar/
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2); library(lme4)
})
out_res <- "results/calendar"
out_fig <- "figures/calendar"
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
B <- as.integer(Sys.getenv("CAL_B", "1000"))

SIR <- c("Susceptible", "Intermediate", "Resistant")
NS  <- c("Resistant", "Intermediate")
DRUGS <- c(levofloxacin = "Fluoroquinolone", tobramycin = "Aminoglycoside", amikacin = "Aminoglycoside",
           gentamicin = "Aminoglycoside", meropenem = "Carbapenem", cefepime = "Cephalosporin",
           ceftazidime = "Cephalosporin", piperacillin_tazobactam = "Penicillin", aztreonam = "Monobactam")
CLASS_DRUGS <- list(
  Aminoglycoside  = c("amikacin", "tobramycin", "gentamicin"),
  Carbapenem      = c("meropenem", "imipenem"),
  Fluoroquinolone = c("ciprofloxacin", "levofloxacin"),
  Cephalosporin   = c("cefepime", "ceftazidime"))
CLSI_FQ <- as.Date("2019-01-01")
T0 <- as.Date("2017-01-01")

iso <- readRDS("data/pa_isolates.rds") %>%
  mutate(id = as.character(individual_id), date = as.Date(date_collected_from_subject),
         half = sprintf("%s-H%d", format(date, "%Y"), ifelse(as.integer(format(date, "%m")) <= 6, 1, 2)),
         period = case_when(date < CLSI_FQ ~ "2017-2018", date < as.Date("2020-01-01") ~ "2019",
                            TRUE ~ "2020-2022"),
         status = coalesce(clinical_status, "unknown"),
         mucoid = pseudomonas_aeruginosa_phenotype == "mucoid")
iso$tested_any <- rowSums(sapply(names(DRUGS), function(d) iso[[d]] %in% SIR)) > 0


# #############################################################################
# 1  SAMPLING DESIGN BY HALF-YEAR
# #############################################################################

cultures <- iso %>% group_by(id, date, half) %>%
  summarise(colonies_all = n(), colonies_tested = sum(tested_any),
            status = first(status), sample_type = first(sample_type), .groups = "drop") %>%
  arrange(id, date) %>% group_by(id) %>% mutate(gap = as.numeric(date - lag(date))) %>% ungroup()

design <- cultures %>% group_by(half) %>%
  summarise(patients = n_distinct(id), cultures = n(),
            cultures_per_patient = cultures / patients,
            colonies_tested_mean = mean(colonies_tested[colonies_tested > 0]),
            pct_single_tested = 100 * mean(colonies_tested[colonies_tested > 0] == 1),
            pct_cultures_untested = 100 * mean(colonies_tested == 0),
            median_gap_days = median(gap, na.rm = TRUE),
            pct_APE_IV = 100 * mean(status == "APE IV Antibiotics"),
            pct_APE_out = 100 * mean(status == "APE Outpatient"),
            pct_stable = 100 * mean(status == "Clinically Stable"),
            pct_sputum = 100 * mean(sample_type %in% "Sputum"), .groups = "drop") %>%
  left_join(iso %>% filter(tested_any) %>% group_by(half) %>%
              summarise(pct_mucoid_tested = 100 * mean(mucoid, na.rm = TRUE),
                        pct_levo_tested = 100 * mean(levofloxacin %in% SIR),
                        pct_aztreonam_tested = 100 * mean(aztreonam %in% SIR), .groups = "drop"),
            by = "half")
write_csv(design, file.path(out_res, "design_by_half_year.csv"))
cat("\nSampling design by half-year:\n"); print(as.data.frame(design), digits = 3, row.names = FALSE)


# #############################################################################
# 2  CHANGE-POINT SCAN, every drug
# #############################################################################

col_long <- map_dfr(names(DRUGS), function(d) {
  iso %>% filter(.data[[d]] %in% SIR) %>%
    transmute(id, date, half, status, mucoid, drug = d, class = DRUGS[[d]], sir = .data[[d]],
              ns = as.integer(sir %in% NS), r = as.integer(sir == "Resistant"),
              s = as.integer(sir == "Susceptible"),
              t = as.numeric(date - T0) / 365.25)
})

composition <- col_long %>% group_by(drug, half) %>%
  summarise(n = n(), S = mean(sir == "Susceptible"), I = mean(sir == "Intermediate"),
            R = mean(sir == "Resistant"), .groups = "drop")
write_csv(composition, file.path(out_res, "sir_composition_by_half_year.csv"))

CANDIDATES <- seq(as.Date("2017-07-01"), as.Date("2021-12-01"), by = "month")
ctrl <- glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)
fit_ll <- function(f, d) {
  m <- try(suppressMessages(glmer(f, data = d, family = binomial, nAGQ = 0, control = ctrl)), silent = TRUE)
  if (inherits(m, "try-error")) NA_real_ else as.numeric(logLik(m))
}

scan <- map_dfr(names(DRUGS), function(dr) {
  d <- col_long %>% filter(drug == dr)
  ll0 <- fit_ll(ns ~ 1 + (1 | id), d)
  map_dfr(CANDIDATES, function(cc) {
    d$step <- as.integer(d$date >= cc)
    if (min(table(d$step)) < 50) return(NULL)
    tibble(drug = dr, change = cc, ll = fit_ll(ns ~ step + (1 | id), d), ll0 = ll0)
  })
}) %>% mutate(lr = 2 * (ll - ll0))
write_csv(scan, file.path(out_res, "changepoint_scan.csv"))

# best change point per drug; p for the MAX over the scan by within-patient
# permutation of dates is expensive, so report the LR at the best c and at the
# CLSI date, and calibrate the max-LR against a parametric bootstrap below
best <- scan %>% group_by(drug) %>% slice_max(lr, n = 1, with_ties = FALSE) %>% ungroup()

# Max-LR null: simulate from the no-step model (patient intercepts kept),
# rescan on a coarse grid of quarters. 100 sims per drug.
NULL_SIMS <- as.integer(Sys.getenv("CAL_NULL", "100"))
Q_CAND <- seq(as.Date("2017-07-01"), as.Date("2021-10-01"), by = "quarter")
maxlr_null <- map_dfr(names(DRUGS), function(dr) {
  d <- col_long %>% filter(drug == dr)
  m0 <- suppressMessages(glmer(ns ~ 1 + (1 | id), data = d, family = binomial, nAGQ = 0, control = ctrl))
  p0 <- fitted(m0)
  mx <- replicate(NULL_SIMS, {
    d$ns <- rbinom(nrow(d), 1, p0)
    l0 <- fit_ll(ns ~ 1 + (1 | id), d)
    max(map_dbl(Q_CAND, function(cc) {
      d$step <- as.integer(d$date >= cc)
      if (min(table(d$step)) < 50) return(NA_real_)
      2 * (fit_ll(ns ~ step + (1 | id), d) - l0)
    }), na.rm = TRUE)
  })
  tibble(drug = dr, null_maxlr_95 = quantile(mx, 0.95), null_maxlr = list(mx))
})
best <- best %>% left_join(maxlr_null %>% select(drug, null_maxlr_95), by = "drug") %>%
  mutate(p_maxlr = map2_dbl(drug, lr, function(dr, l) {
    mx <- maxlr_null$null_maxlr[[which(maxlr_null$drug == dr)]]; (1 + sum(mx >= l)) / (1 + length(mx))
  }))

# step size at the best date and at the CLSI date, Laplace fits, with S and R
step_fit <- function(dr, cc, outcome) {
  d <- col_long %>% filter(drug == dr) %>% mutate(step = as.integer(date >= cc))
  d$y <- d[[outcome]]
  m <- suppressMessages(glmer(y ~ step + (1 | id), data = d, family = binomial, control = ctrl))
  s <- summary(m)$coefficients["step", ]
  tibble(drug = dr, change = cc, outcome = outcome, log_or = s[["Estimate"]], se = s[["Std. Error"]],
         or = exp(log_or), or_lo = exp(log_or - 1.96 * se), or_hi = exp(log_or + 1.96 * se))
}
steps <- bind_rows(
  pmap_dfr(best %>% select(drug, change), function(drug, change)
    map_dfr(c("ns", "s", "r"), ~ step_fit(drug, change, .x)) %>% mutate(at = "best change point")),
  map_dfr(names(DRUGS), function(dr)
    map_dfr(c("ns", "s", "r"), ~ step_fit(dr, CLSI_FQ, .x)) %>% mutate(at = "CLSI 2019-01-01")))
write_csv(best %>% select(-any_of("null_maxlr")), file.path(out_res, "changepoint_best.csv"))
# Step or drift? For every drug: patient-only, linear trend, step at the drug's
# best date, step at the CLSI date, trend + best step (AIC; patient intercept).
shape <- map_dfr(names(DRUGS), function(dr) {
  d <- col_long %>% filter(drug == dr) %>%
    mutate(step_best = as.integer(date >= best$change[best$drug == dr]),
           step_clsi = as.integer(date >= CLSI_FQ))
  fs <- list(none = ns ~ 1 + (1 | id), trend = ns ~ t + (1 | id), step_best = ns ~ step_best + (1 | id),
             step_clsi = ns ~ step_clsi + (1 | id), trend_step = ns ~ t + step_best + (1 | id))
  ms <- map(fs, ~ suppressMessages(glmer(.x, data = d, family = binomial, control = ctrl)))
  cf <- summary(ms$trend_step)$coefficients
  tibble(drug = dr, best_change = best$change[best$drug == dr],
         dAIC_trend = AIC(ms$trend) - AIC(ms$none), dAIC_step_best = AIC(ms$step_best) - AIC(ms$none),
         dAIC_step_clsi = AIC(ms$step_clsi) - AIC(ms$none), dAIC_trend_step = AIC(ms$trend_step) - AIC(ms$none),
         trend_OR_given_step = exp(cf["t", 1]), trend_p_given_step = 2 * pnorm(-abs(cf["t", 3])),
         step_OR_given_trend = exp(cf["step_best", 1]), step_p_given_trend = 2 * pnorm(-abs(cf["step_best", 3])))
})
write_csv(shape, file.path(out_res, "step_vs_trend_by_drug.csv"))
cat("\nStep or drift, per drug (AIC relative to patient-only):\n"); print(as.data.frame(shape), digits = 3, row.names = FALSE)

write_csv(steps, file.path(out_res, "changepoint_steps.csv"))
cat("\nBest change point per drug:\n"); print(as.data.frame(best), digits = 3, row.names = FALSE)


# #############################################################################
# 3  FLUOROQUINOLONE: trend vs step vs both, adjusted
# #############################################################################

levo <- col_long %>% filter(drug == "levofloxacin")
c_star <- best$change[best$drug == "levofloxacin"]
levo <- levo %>% mutate(step_best = as.integer(date >= c_star), step_clsi = as.integer(date >= CLSI_FQ),
                        status = factor(status), mucoid = factor(mucoid))
fq_models <- list(
  `patient only`                = ns ~ 1 + (1 | id),
  `linear calendar trend`       = ns ~ t + (1 | id),
  `step at best date`           = ns ~ step_best + (1 | id),
  `step at CLSI date`           = ns ~ step_clsi + (1 | id),
  `trend + step (best)`         = ns ~ t + step_best + (1 | id),
  `trend + step (best), adjusted` = ns ~ t + step_best + status + mucoid + (1 | id),
  `trend, adjusted`             = ns ~ t + status + mucoid + (1 | id))
fq_fits <- imap(fq_models, function(f, nm) suppressMessages(glmer(f, data = levo, family = binomial, control = ctrl)))
fq_tab <- imap_dfr(fq_fits, function(m, nm) {
  cf <- summary(m)$coefficients
  get <- function(term) if (term %in% rownames(cf))
    sprintf("%.2f (%.2f–%.2f)", exp(cf[term, 1]), exp(cf[term, 1] - 1.96 * cf[term, 2]),
            exp(cf[term, 1] + 1.96 * cf[term, 2])) else ""
  tibble(model = nm, AIC = AIC(m), `OR per year` = get("t"),
         `OR step` = paste0(get("step_best"), get("step_clsi")),
         sd_patient = sqrt(VarCorr(m)$id[1]))
}) %>% mutate(dAIC = AIC - min(AIC))
write_csv(fq_tab, file.path(out_res, "fq_trend_models.csv"))
cat("\nLevofloxacin models:\n"); print(as.data.frame(fq_tab), digits = 3, row.names = FALSE)

# raw cohort mean vs within-patient: patient-weighted yearly mean, and the
# yearly mean of each patient's deviation from their own overall mean
fq_raw <- levo %>% mutate(year = format(date, "%Y")) %>%
  group_by(id) %>% mutate(dev = ns - mean(ns)) %>% ungroup() %>%
  group_by(year, id) %>% summarise(ns = mean(ns), dev = mean(dev), .groups = "drop") %>%
  group_by(year) %>% summarise(patients = n(), cohort_mean_NS = mean(ns),
                               within_patient_deviation = mean(dev), .groups = "drop")
write_csv(fq_raw, file.path(out_res, "fq_raw_vs_within.csv"))
cat("\nLevofloxacin by year, cohort vs within-patient:\n"); print(as.data.frame(fq_raw), digits = 3)

# who leaves: mean levo NS of patients by the year of their last culture
attrition <- levo %>% group_by(id) %>% summarise(last = format(max(date), "%Y"), ns = mean(ns)) %>%
  group_by(last) %>% summarise(patients = n(), mean_NS = mean(ns), .groups = "drop")
write_csv(attrition, file.path(out_res, "fq_attrition.csv"))


# #############################################################################
# 4  SENSITIVITY OF EARLIER ANALYSES TO PERIOD
# #############################################################################

colony_frame <- imap_dfr(CLASS_DRUGS, function(drugs, cl) {
  tests <- rowSums(sapply(drugs, function(x) iso[[x]] %in% SIR))
  ns    <- rowSums(sapply(drugs, function(x) iso[[x]] %in% NS)) > 0
  tibble(id = iso$id, date = iso$date, period = iso$period, tests, ns) %>%
    filter(tests > 0) %>% group_by(id, date, period) %>%
    summarise(m = n(), k = sum(ns), .groups = "drop") %>% mutate(class = cl)
})

# (a) within-culture discordance by period, culture-weighted, patient bootstrap
disc_stat <- function(d) d %>% filter(m >= 2) %>%
  mutate(dd = 2 * k * (m - k) / (m * (m - 1))) %>%
  group_by(class, period) %>% summarise(disc = mean(dd), cultures = n(), .groups = "drop")
disc_pt <- disc_stat(colony_frame)
by_id <- split(colony_frame, colony_frame$id); ids <- names(by_id)
disc_bt <- map_dfr(seq_len(B), function(b) disc_stat(bind_rows(by_id[sample(ids, length(ids), TRUE)])))
disc <- disc_pt %>% left_join(disc_bt %>% group_by(class, period) %>%
                                summarise(lo = quantile(disc, 0.025), hi = quantile(disc, 0.975), .groups = "drop"),
                              by = c("class", "period"))
write_csv(disc, file.path(out_res, "discordance_by_period.csv"))

# (b) analysis-1 within-patient variance, FQ colony frame, centred on the
# patient's mean (as analysis 1) vs on the patient's mean within each side of
# the change point. Pooled-E moment estimator, m >= 2 cultures.
vw <- function(d, centre_by) {
  d <- d %>% filter(m >= 2) %>% mutate(q = k / m)
  mu <- mean(d$q); Einv <- mean(1 / d$m)
  V_tot <- (var(d$q) - mu * (1 - mu) * Einv) / (1 - Einv); E <- mu * (1 - mu) - V_tot
  w <- d %>% group_by(across(all_of(centre_by))) %>% filter(n() >= 2) %>%
    summarise(n = n(), ss = sum((q - mean(q))^2), s = mean(E / m), .groups = "drop")
  tibble(V_tot = V_tot, V_w = (sum(w$ss) - sum((w$n - 1) * w$s)) / sum(w$n - 1),
         share = V_w / V_tot, R_wp_1 = ifelse(V_w > 0, V_w / (V_w + E), 0))
}
fq_cf <- colony_frame %>% filter(class == "Fluoroquinolone") %>%
  mutate(side = ifelse(date >= c_star, "after", "before"))
vw_tab <- bind_rows(
  vw(fq_cf, "id") %>% mutate(centring = "patient (analysis 1)"),
  vw(fq_cf, c("id", "side")) %>% mutate(centring = paste("patient x side of", c_star)))
fq_by_id <- split(fq_cf, fq_cf$id)
vw_bt <- map_dfr(seq_len(B), function(b) {
  pick <- sample(names(fq_by_id), replace = TRUE)
  db <- bind_rows(map2(fq_by_id[pick], seq_along(pick), ~ mutate(.x, id = as.character(.y))))
  bind_rows(vw(db, "id") %>% mutate(centring = "patient (analysis 1)"),
            vw(db, c("id", "side")) %>% mutate(centring = paste("patient x side of", c_star))) %>% mutate(rep = b)
})
vw_tab <- vw_tab %>% left_join(vw_bt %>% group_by(centring) %>%
                                 summarise(share_lo = quantile(share, 0.025), share_hi = quantile(share, 0.975),
                                           .groups = "drop"), by = "centring")
write_csv(vw_tab, file.path(out_res, "fq_within_patient_share_by_centring.csv"))


# #############################################################################
# FIGURE: drawn by analyses/figures_new.R from the CSVs written here, at print
# size (analyses/fig_style.R). This script no longer draws.
# #############################################################################

# #############################################################################
# REPORT
# #############################################################################

md_tab <- function(df) {
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
fmt <- function(df) df %>% mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))
writeLines(c(
  "# Calendar effects", "",
  "## Sampling design by half-year", "", md_tab(fmt(design)),
  "", "## Change-point scan: best step date per drug", "",
  sprintf("Max-LR null: %d parametric simulations from the no-step model, quarterly grid.", NULL_SIMS), "",
  md_tab(best %>% transmute(drug, best_change = format(change, "%Y-%m"), LR = sprintf("%.1f", lr),
                            `null 95% max-LR` = sprintf("%.1f", null_maxlr_95), p = sprintf("%.3f", p_maxlr))),
  "", "## Step or drift: AIC relative to patient-only, per drug", "",
  md_tab(shape %>% transmute(drug, best_change = format(best_change, "%Y-%m"),
                             `dAIC trend` = sprintf("%.1f", dAIC_trend), `dAIC step (best)` = sprintf("%.1f", dAIC_step_best),
                             `dAIC step (CLSI)` = sprintf("%.1f", dAIC_step_clsi), `dAIC trend+step` = sprintf("%.1f", dAIC_trend_step),
                             `trend OR/yr given step` = sprintf("%.2f (p %.3f)", trend_OR_given_step, trend_p_given_step),
                             `step OR given trend` = sprintf("%.2f (p %.3f)", step_OR_given_trend, step_p_given_trend))),
  "", "## Step odds ratios (patient random intercept)", "",
  md_tab(steps %>% transmute(at, drug, change = format(change, "%Y-%m"), outcome,
                             OR = sprintf("%.2f (%.2f–%.2f)", or, or_lo, or_hi))),
  "", "## Levofloxacin: trend vs step", "", md_tab(fmt(fq_tab)),
  "", "## Levofloxacin by year: cohort mean vs within-patient deviation", "", md_tab(fmt(fq_raw)),
  "", "## Levofloxacin NS by year of a patient's last culture (attrition)", "", md_tab(fmt(attrition)),
  "", "## Within-culture discordance by period (colony frame, culture-weighted)", "",
  md_tab(disc %>% transmute(class, period, cultures, disc = sprintf("%.2f (%.2f–%.2f)", disc, lo, hi))),
  "", "## Analysis-1 within-patient share, fluoroquinolone, by centring", "",
  md_tab(vw_tab %>% transmute(centring, V_tot = sprintf("%.3f", V_tot), V_w = sprintf("%.3f", V_w),
                              share = sprintf("%.2f (%.2f–%.2f)", share, share_lo, share_hi),
                              `R_wp at 1 colony` = sprintf("%.2f", R_wp_1)))),
  file.path(out_res, "calendar_effects.md"))
cat("\nwrote", out_res, "and", out_fig, "\n")
