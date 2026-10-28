# =============================================================================
# analyses/fit_power_frontiers.R
#
# Power under the recalibrated model (analyses/fit_power_calibration.R,
# "selection" version), with the levers separated and costs held fixed.
#
# Latent model per class: patient level mu_i ~ N(mu0, sd_pat^2); stationary OU
# x_i(t) with SD sd_w and memory tau; theta = plogis(mu_i + x_i(t)); m colonies
# picked at RANDOM per visit (a designed study chooses its depth; fractional m
# = each visit floor(m) or ceiling(m)). Regular visit schedule.
#
# Estimands (cohort level, alpha = 0.05, two-sided):
#   step      every patient's latent level shifts by DELTA (logit) halfway
#             through follow-up (e.g. a new therapy). Detector: per-patient
#             post-minus-pre mean NS fraction, one-sample t across patients.
#   coupling  class X drives class Y in CONTINUOUS time,
#             dy = (-y / tau + kappa x) dt + noise, so the lag coefficient at
#             interval dt is beta(dt) = kappa dt exp(-dt / tau): a denser
#             schedule sees a smaller per-step effect, which the 14.3 grid
#             (beta fixed per step) did not. Exact discretisation of the 2-D
#             OU. Detector: within-patient (fixed-effect) cross-lagged OLS,
#             py ~ lag py + lag px, cluster-robust SE by patient.
# Effect sizes are set so the cohort's own design (N = 90, L = 12, dt = 84,
# m = 1.5) has power ~ 0.5, off floor and ceiling, then held fixed.
#
# Levers, from the baseline N = 90, L = 12, dt = 84 d, m = 1.5:
#   interval   dt varies, L fixed: follow-up T = 11 dt shrinks with dt
#   visits     L varies, T fixed at 924 d: dt = T / (L - 1)
#   follow-up  L varies, dt fixed at 84 d: T = 84 (L - 1)
#   depth      m varies
# each (a) with N = 90 (what the lever buys) and (b) at FIXED COST, N set by
# budget = 90 x cost(L = 12, m = 1.5), cost per patient = c_p + L (c_v + m c_c).
# Cost weights are assumptions: base c_c = 1 (one colony AST), c_v = 2 (a
# culture), c_p = 20 (enrolling a patient); "visits dear" c_v = 10.
# The interval lever costs nothing per patient; its price is calendar time.
#
# Sensitivity: tau at the ends of its profile (60, 480 d); cephalosporin
# parameters in place of aminoglycoside.
#
# Run from the repo root:  Rscript analyses/fit_power_frontiers.R
# Output: results/power_recalibrated/   figures/power_recalibrated/
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2)
})
out_res <- "results/power_recalibrated"
out_fig <- "figures/power_recalibrated"
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
REPS <- as.integer(Sys.getenv("POWER_REPS", "1000"))

cal <- read_csv(file.path(out_res, "calibration_parameters.csv"), show_col_types = FALSE) %>%
  filter(version == "selection")
par_of <- function(cl, tau = NULL) {
  r <- cal %>% filter(class == cl)
  list(mu0 = r$mu0, sd_pat = r$sd_pat, sd_w = r$sd_w, tau = if (is.null(tau)) r$tau_days else tau)
}
BASE <- list(N = 90, L = 12, dt = 84, m = 1.5)
COSTS <- list(base = c(c_c = 1, c_v = 2, c_p = 20), `visits dear` = c(c_c = 1, c_v = 10, c_p = 20))
cost_pp <- function(L, m, cs) cs[["c_p"]] + L * (cs[["c_v"]] + m * cs[["c_c"]])


# #############################################################################
# SIMULATORS
# #############################################################################

draw_m <- function(N, L, m) {
  lo <- floor(m); matrix(lo + (runif(N * L) < m - lo), N, L)
}

sim_ou <- function(N, L, dt, sd_w, tau) {
  rho <- exp(-dt / tau)
  X <- matrix(0, N, L); X[, 1] <- rnorm(N, 0, sd_w)
  for (j in seq_len(L)[-1]) X[, j] <- rho * X[, j - 1] + sd_w * sqrt(1 - rho^2) * rnorm(N)
  X
}

# coupled 2-D OU, exact discretisation; x drives y at rate kappa per day
sim_coupled <- function(N, L, dt, sd_w, tau, kappa) {
  a <- 1 / tau; q <- 2 * a * sd_w^2
  E <- exp(-a * dt) * matrix(c(1, kappa * dt, 0, 1), 2)          # e^{A dt}
  e2 <- exp(-2 * a * dt)
  I0 <- (1 - e2) / (2 * a)
  I1 <- (1 - e2 * (1 + 2 * a * dt)) / (4 * a^2)
  I2 <- (2 - e2 * (2 + 4 * a * dt + 4 * a^2 * dt^2)) / (8 * a^3)
  Q  <- q * matrix(c(I0, kappa * I1, kappa * I1, I0 + kappa^2 * I2), 2)
  P  <- q * matrix(c(1 / (2 * a), kappa / (4 * a^2), kappa / (4 * a^2),
                     1 / (2 * a) + kappa^2 / (4 * a^3)), 2)          # stationary
  cq <- chol(Q); cp <- chol(P)
  Z <- matrix(rnorm(2 * N), N) %*% cp
  X <- Y <- matrix(0, N, L); X[, 1] <- Z[, 1]; Y[, 1] <- Z[, 2]
  for (j in seq_len(L)[-1]) {
    eps <- matrix(rnorm(2 * N), N) %*% cq
    X[, j] <- E[1, 1] * X[, j - 1] + eps[, 1]
    Y[, j] <- E[2, 1] * X[, j - 1] + E[2, 2] * Y[, j - 1] + eps[, 2]
  }
  list(X = X, Y = Y)
}

observe <- function(lat, mu, m_mat) {
  K <- matrix(rbinom(length(m_mat), m_mat, plogis(mu + lat)), nrow(m_mat))
  K / m_mat
}


# #############################################################################
# DETECTORS
# #############################################################################

p_step <- function(Q, L) {
  pre <- seq_len(floor(L / 2)); post <- setdiff(seq_len(L), pre)
  d <- rowMeans(Q[, post, drop = FALSE]) - rowMeans(Q[, pre, drop = FALSE])
  s <- sd(d); if (s == 0) return(1)
  2 * pt(-abs(mean(d) / (s / sqrt(length(d)))), length(d) - 1)
}

p_coupling <- function(QX, QY) {
  N <- nrow(QY); L <- ncol(QY)
  y <- QY[, -1]; ly <- QY[, -L]; lx <- QX[, -L]
  dm <- function(M) M - rowMeans(M)                  # within-patient demeaning
  yv <- as.vector(dm(y)); X <- cbind(as.vector(dm(ly)), as.vector(dm(lx)))
  id <- rep(seq_len(N), times = L - 1)
  XtX <- crossprod(X); if (rcond(XtX) < 1e-10) return(NA_real_)
  bread <- solve(XtX); b <- bread %*% crossprod(X, yv)
  e <- yv - X %*% b
  sc <- rowsum(X * as.vector(e), id)                 # per-patient scores
  V <- bread %*% crossprod(sc) %*% bread * N / (N - 1)
  2 * pt(-abs(b[2] / sqrt(V[2, 2])), N - 1)
}

power_step <- function(N, L, dt, m, pars, delta, reps = REPS) {
  pv <- replicate(reps, {
    mu <- rnorm(N, pars$mu0, pars$sd_pat)
    lat <- sim_ou(N, L, dt, pars$sd_w, pars$tau)
    lat[, -seq_len(floor(L / 2))] <- lat[, -seq_len(floor(L / 2))] + delta
    p_step(observe(lat, mu, draw_m(N, L, m)), L)
  })
  mean(pv < 0.05)
}

power_coupling <- function(N, L, dt, m, pars, kappa, reps = REPS) {
  pv <- replicate(reps, {
    lat <- sim_coupled(N, L, dt, pars$sd_w, pars$tau, kappa)
    mm <- draw_m(N, L, m)                            # one set of colonies, both classes
    QX <- observe(lat$X, rnorm(N, pars$mu0, pars$sd_pat), mm)
    QY <- observe(lat$Y, rnorm(N, pars$mu0, pars$sd_pat), mm)
    p_coupling(QX, QY)
  })
  mean(pv < 0.05, na.rm = TRUE)
}


# #############################################################################
# 1  EFFECT SIZES: power ~ 0.5 at the cohort's own design (aminoglycoside)
# #############################################################################

PAR <- par_of("Aminoglycoside")
cat(sprintf("Aminoglycoside: mu0 %.2f, sd_pat %.2f, sd_w %.2f, tau %.0f d\n",
            PAR$mu0, PAR$sd_pat, PAR$sd_w, PAR$tau))

fpr <- tibble(
  estimand = c("step", "coupling"),
  false_positive = c(power_step(90, 12, 84, 1.5, PAR, 0, 2000),
                     power_coupling(90, 12, 84, 1.5, PAR, 0, 2000)))
cat("\nFalse-positive rate at the baseline design (2000 reps):\n"); print(fpr)

pick_effect <- function(f, grid) {
  pw <- map_dbl(grid, f)
  cat("   ", paste(sprintf("%.4g -> %.2f", grid, pw), collapse = " | "), "\n")
  grid[which.min(abs(pw - 0.5))]
}
cat("\nstep delta (logit):\n")
DELTA <- pick_effect(function(d) power_step(90, 12, 84, 1.5, PAR, d, 1000),
                     c(0.25, 0.5, 0.75, 1, 1.5, 2))
cat("coupling beta at 84 d (kappa = beta / (84 exp(-84 / tau))):\n")
k_of <- function(beta84, tau = PAR$tau) beta84 / (84 * exp(-84 / tau))
B84 <- pick_effect(function(b) power_coupling(90, 12, 84, 1.5, PAR, k_of(b), 1000),
                   c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4))
KAPPA <- k_of(B84)
effects <- tibble(delta_logit = DELTA, beta_at_84d = B84, kappa_per_day = KAPPA, tau = PAR$tau)
write_csv(bind_cols(effects, fpr %>% pivot_wider(names_from = estimand, values_from = false_positive,
                                                 names_prefix = "fpr_")),
          file.path(out_res, "power_effect_sizes.csv"))


# #############################################################################
# 2  LEVERS: free (N = 90) and at fixed cost
# #############################################################################

T0 <- (BASE$L - 1) * BASE$dt
LEVELS_L <- c(4, 6, 8, 12, 16, 24, 32, 48)
levers <- bind_rows(
  tibble(lever = "interval (L fixed)", dt = c(14, 28, 42, 56, 84, 126, 168),
         L = BASE$L, m = BASE$m, x = dt),
  tibble(lever = "visits (follow-up fixed)", L = LEVELS_L, dt = T0 / (L - 1), m = BASE$m, x = L),
  tibble(lever = "follow-up (interval fixed)", L = LEVELS_L, dt = BASE$dt, m = BASE$m, x = L),
  tibble(lever = "depth", m = c(1, 1.5, 2, 3, 4, 6, 8), L = BASE$L, dt = BASE$dt, x = m)) %>%
  mutate(followup_days = (L - 1) * dt)

designs <- bind_rows(
  levers %>% mutate(budget = "N fixed at 90", N = 90),
  imap_dfr(COSTS, function(cs, nm) {
    B <- 90 * cost_pp(BASE$L, BASE$m, cs)
    levers %>% mutate(budget = paste("fixed cost:", nm), N = floor(B / cost_pp(L, m, cs)))
  }))

run_designs <- function(des, pars, label) {
  t0 <- Sys.time()
  out <- des %>% mutate(
    power_step = pmap_dbl(list(N, L, dt, m), function(N, L, dt, m) power_step(N, L, dt, m, pars, DELTA)),
    power_coupling = pmap_dbl(list(N, L, dt, m), function(N, L, dt, m)
      power_coupling(N, L, dt, m, pars, KAPPA)),   # kappa is the physical rate: held fixed
    params = label)
  cat(sprintf("  %-34s %d cells  %.0f s\n", label, nrow(des), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  out
}

main <- run_designs(designs, PAR, "Aminoglycoside (fitted tau)")
sens_des <- designs %>% filter(budget %in% c("N fixed at 90", "fixed cost: base"))
sens <- bind_rows(
  run_designs(sens_des, par_of("Aminoglycoside", 60),  "Aminoglycoside, tau = 60 d"),
  run_designs(sens_des, par_of("Aminoglycoside", 480), "Aminoglycoside, tau = 480 d"),
  run_designs(sens_des, par_of("Cephalosporin"),       "Cephalosporin (fitted)"))
grid_all <- bind_rows(main, sens) %>%
  mutate(se_step = sqrt(power_step * (1 - power_step) / REPS),
         se_coupling = sqrt(power_coupling * (1 - power_coupling) / REPS))
write_csv(grid_all, file.path(out_res, "power_frontiers.csv"))


# #############################################################################
# 3  MATCHED 2x LEVERS (the Table 2 question, recalibrated)
# #############################################################################

matched <- tribble(
  ~lever,                          ~N,  ~L, ~dt, ~m,
  "baseline (as run)",             90,  12, 84,  1.5,
  "cohort 90 -> 180",             180,  12, 84,  1.5,
  "visits 12 -> 24, same follow-up", 90, 24, 40,  1.5,
  "follow-up doubled (24 visits, 84 d)", 90, 24, 84, 1.5,
  "depth 1.5 -> 3 colonies",       90,  12, 84,  3,
  "interval 84 -> 42 d, same visits", 90, 12, 42, 1.5)
matched_cost <- matched %>% mutate(
  N = floor(90 * cost_pp(12, 1.5, COSTS$base) / cost_pp(L, m, COSTS$base)),
  N = ifelse(lever == "cohort 90 -> 180", NA, N))
tab2 <- bind_rows(matched %>% mutate(budget = "free"),
                  matched_cost %>% filter(!is.na(N)) %>% mutate(budget = "fixed cost (base)")) %>%
  mutate(power_step = pmap_dbl(list(N, L, dt, m), function(N, L, dt, m) power_step(N, L, dt, m, PAR, DELTA)),
         power_coupling = pmap_dbl(list(N, L, dt, m), function(N, L, dt, m) power_coupling(N, L, dt, m, PAR, KAPPA)),
         followup_days = (L - 1) * dt) %>%
  group_by(budget) %>%
  mutate(gain_step = power_step - power_step[lever == "baseline (as run)"],
         gain_coupling = power_coupling - power_coupling[lever == "baseline (as run)"]) %>% ungroup()
write_csv(tab2, file.path(out_res, "power_matched_levers.csv"))


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
f2 <- function(x) sprintf("%.2f", x)
sg <- function(x) sprintf("%+.2f", x)
lever_tab <- grid_all %>%
  transmute(params, budget, lever, N, L, dt = round(dt), m, followup_days = round(followup_days),
            step = f2(power_step), coupling = f2(power_coupling))

writeLines(c(
  "# Power, recalibrated",
  "", sprintf("%d reps per cell (MC SE <= %.3f). Aminoglycoside calibration unless stated.", REPS, sqrt(0.25 / REPS)),
  sprintf("Effect sizes: step = %.2g logit; coupling beta = %.2g at 84 d (kappa = %.2e /day, tau = %.0f d).",
          DELTA, B84, KAPPA, PAR$tau),
  sprintf("False-positive rate at baseline (2000 reps): step %.3f, coupling %.3f.", fpr$false_positive[1], fpr$false_positive[2]),
  "", "## Matched 2x levers (Table 2, recalibrated)", "",
  md_tab(tab2 %>% transmute(budget, lever, N, L, dt = round(dt), m, followup_days = round(followup_days),
                            step = f2(power_step), `gain step` = sg(gain_step),
                            coupling = f2(power_coupling), `gain coupling` = sg(gain_coupling))),
  "", "## All lever cells", "", md_tab(lever_tab)),
  file.path(out_res, "power_recalibrated.md"))
cat("\nwrote", out_res, "and", out_fig, "\n")
