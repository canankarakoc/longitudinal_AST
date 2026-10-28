# =============================================================================
# analyses/fit_calendar_ou_null.R
#
# Companion to fit_calendar_effects.R. Its max-LR null simulated each drug
# from a patient-intercept model with INDEPENDENT cultures, so any ordinary
# within-patient persistence reads as a "step": every drug but gentamicin beat
# that null. Here the null keeps the dynamics.
#
#   A  Class-level change-point scan (colony frame: a colony is NS if any
#      agent in the class reads R or I; culture-level binomial GLMM, patient
#      intercept, step at quarter c). Null: the analysis-5 calibrated model
#      ("selection" version: patient level + OU within-patient process with
#      tau ~ 130-240 d + heterogeneity-dependent depth) simulated on the
#      cohort's own dates, rescanned. The calibration was fitted on data that
#      include any step, so this null is if anything generous to "no step".
#      Also the LR at the CLSI date (2019-01-01) against the same null.
#   B  Are the 2020-on drops (tobramycin, cefepime) the shift to clinically
#      stable visits? Step at 2020-03-01 with and without clinical status.
#
# Run from the repo root:  Rscript analyses/fit_calendar_ou_null.R
# Output: results/calendar/ou_null_*.csv, status_adjusted_2020_steps.csv
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(lme4)
})
out_res <- "results/calendar"
set.seed(42)
NSIM <- as.integer(Sys.getenv("OU_NSIM", "200"))

CLASSES <- c("Aminoglycoside", "Carbapenem", "Fluoroquinolone", "Cephalosporin")
CLASS_DRUGS <- list(
  Aminoglycoside  = c("amikacin", "tobramycin", "gentamicin"),
  Carbapenem      = c("meropenem", "imipenem"),
  Fluoroquinolone = c("ciprofloxacin", "levofloxacin"),
  Cephalosporin   = c("cefepime", "ceftazidime"))
SIR <- c("Susceptible", "Intermediate", "Resistant")
NS  <- c("Resistant", "Intermediate")
CLSI_FQ <- as.Date("2019-01-01")
Q_CAND <- seq(as.Date("2017-07-01"), as.Date("2021-10-01"), by = "quarter")
ctrl <- glmerControl(optimizer = "bobyqa", calc.derivs = FALSE)

iso <- readRDS("data/pa_isolates.rds")
cult <- imap_dfr(CLASS_DRUGS, function(drugs, cl) {
  tests <- rowSums(sapply(drugs, function(x) iso[[x]] %in% SIR))
  ns    <- rowSums(sapply(drugs, function(x) iso[[x]] %in% NS)) > 0
  tibble(id = as.character(iso$individual_id), date = as.Date(iso$date_collected_from_subject),
         status = coalesce(iso$clinical_status, "unknown"), tests, ns) %>%
    filter(tests > 0) %>% group_by(id, date) %>%
    summarise(m = n(), k = sum(ns), status = first(status), .groups = "drop") %>% mutate(class = cl)
}) %>% arrange(class, id, date)

pars <- read_csv("results/power_recalibrated/calibration_parameters.csv", show_col_types = FALSE) %>%
  filter(version == "selection")


# #############################################################################
# A  SCAN WITH A DYNAMICS-PRESERVING NULL
# #############################################################################

ll_step <- function(d, cc) {
  d$step <- as.integer(d$date >= cc)
  if (min(table(d$step)) < 30) return(NA_real_)
  m <- try(suppressMessages(glmer(cbind(k, m - k) ~ step + (1 | id), data = d, family = binomial,
                                  nAGQ = 0, control = ctrl)), silent = TRUE)
  if (inherits(m, "try-error")) NA_real_ else as.numeric(logLik(m))
}
ll_null <- function(d) as.numeric(logLik(suppressMessages(
  glmer(cbind(k, m - k) ~ 1 + (1 | id), data = d, family = binomial, nAGQ = 0, control = ctrl))))
scan_lr <- function(d) {
  l0 <- ll_null(d)
  lr <- map_dbl(Q_CAND, ~ 2 * (ll_step(d, .x) - l0))
  c(max = max(lr, na.rm = TRUE), at = Q_CAND[which.max(lr)], clsi = 2 * (ll_step(d, CLSI_FQ) - l0))
}

# simulate the calibrated model on this class's own dates and colony counts
sim_class <- function(d, p) {
  d <- d %>% group_by(id) %>% mutate(gap = as.numeric(date - lag(date))) %>% ungroup()
  ids <- unique(d$id)
  mu <- setNames(rnorm(length(ids), p$mu0, p$sd_pat), ids)
  x <- numeric(nrow(d))
  for (i in seq_len(nrow(d))) {
    if (is.na(d$gap[i])) x[i] <- rnorm(1, 0, p$sd_w)
    else { r <- exp(-d$gap[i] / p$tau_days); x[i] <- r * x[i - 1] + p$sd_w * sqrt(1 - r^2) * rnorm(1) }
  }
  th <- plogis(mu[d$id] + x)
  # heterogeneity-dependent depth, share of multi-colony cultures as observed
  h <- 4 * th * (1 - th); f_multi <- mean(d$m >= 2); md <- d$m[d$m >= 2]
  a <- uniroot(function(a) mean(plogis(a + p$b * h)) - f_multi, c(-60, 60))$root
  multi <- runif(nrow(d)) < plogis(a + p$b * h)
  d$m <- ifelse(multi, sample(md, nrow(d), replace = TRUE), 1L)
  d$k <- rbinom(nrow(d), d$m, th)
  d
}

ou <- map_dfr(CLASSES, function(cl) {
  d <- cult %>% filter(class == cl); p <- pars %>% filter(class == cl)
  obs <- scan_lr(d)
  null <- t(replicate(NSIM, scan_lr(sim_class(d, p))))
  cat(sprintf("  %-16s observed max-LR %.1f at %s, CLSI-date LR %.1f\n", cl, obs[["max"]],
              format(as.Date(obs[["at"]], origin = "1970-01-01"), "%Y-%m"), obs[["clsi"]]))
  tibble(class = cl, max_lr = obs[["max"]], at = as.Date(obs[["at"]], origin = "1970-01-01"),
         clsi_lr = obs[["clsi"]],
         null_max_95 = quantile(null[, "max"], 0.95, na.rm = TRUE),
         p_max = (1 + sum(null[, "max"] >= obs[["max"]], na.rm = TRUE)) / (1 + sum(!is.na(null[, "max"]))),
         null_clsi_95 = quantile(null[, "clsi"], 0.95, na.rm = TRUE),
         p_clsi = (1 + sum(null[, "clsi"] >= obs[["clsi"]], na.rm = TRUE)) / (1 + sum(!is.na(null[, "clsi"]))))
})
write_csv(ou, file.path(out_res, "ou_null_changepoint.csv"))
cat("\nClass-level scan against the calibrated-dynamics null:\n"); print(as.data.frame(ou), digits = 3)


# #############################################################################
# B  2020 STEPS, WITH AND WITHOUT CLINICAL STATUS
# #############################################################################

col_long <- map_dfr(c("tobramycin", "cefepime", "levofloxacin", "meropenem"), function(dr) {
  iso %>% filter(.data[[dr]] %in% SIR) %>%
    transmute(id = as.character(individual_id), date = as.Date(date_collected_from_subject),
              status = factor(coalesce(clinical_status, "unknown")), drug = dr,
              ns = as.integer(.data[[dr]] %in% NS))
})
adj <- map_dfr(unique(col_long$drug), function(dr) {
  d <- col_long %>% filter(drug == dr) %>% mutate(step = as.integer(date >= as.Date("2020-03-01")))
  map_dfr(c(FALSE, TRUE), function(a) {
    f <- if (a) ns ~ step + status + (1 | id) else ns ~ step + (1 | id)
    cf <- summary(suppressMessages(glmer(f, data = d, family = binomial, control = ctrl)))$coefficients["step", ]
    tibble(drug = dr, adjusted_for_status = a,
           OR = exp(cf[["Estimate"]]), lo = exp(cf[["Estimate"]] - 1.96 * cf[["Std. Error"]]),
           hi = exp(cf[["Estimate"]] + 1.96 * cf[["Std. Error"]]))
  })
})
status_ns <- col_long %>% group_by(drug, status) %>% summarise(n = n(), NS = mean(ns), .groups = "drop")
write_csv(adj, file.path(out_res, "status_adjusted_2020_steps.csv"))
write_csv(status_ns, file.path(out_res, "ns_by_clinical_status.csv"))
cat("\nStep at 2020-03, with and without clinical status:\n"); print(as.data.frame(adj), digits = 3)
cat("\nNS by clinical status:\n"); print(as.data.frame(status_ns), digits = 3)
