# =============================================================================
# analyses/fit_power_calibration.R
#
# Recalibrates the generative model behind the power simulation (pipeline 14.3,
# Table 2, Fig 5). 14.3 set the latent AR(1) rho to r_obs(dt) / 0.74, where
# 0.74 is the POOLED reliability; the patient-centred autocorrelation is
# attenuated by the within-patient reliability (~0.20, analysis 1), which would
# imply rho > 1. It also fixed sd_pat = sd_inn = 1 without calibration, started
# each series off its stationary distribution, and held the per-step coupling
# fixed as the interval changed.
#
# Model (per class, colony frame: a colony is NS if any agent in the class
# reads R or I), on the logit scale:
#   mu_i      ~ N(mu0, sd_pat^2)                    patient level
#   x_i(t)    stationary Ornstein-Uhlenbeck, SD sd_w, corr exp(-dt / tau)
#   theta_ij  = plogis(mu_i + x_i(t_ij))
#   k_ij      ~ Binomial(m_ij, theta_ij)            colonies tested
# simulated on the cohort's OWN design: every patient's culture dates and
# colonies per culture.
#
# Two versions:
#   base       depth m_ij as observed, independent of theta (the 14.3 model)
#   selection  depth depends on the culture's heterogeneity: a culture gets
#              >= 2 colonies with probability plogis(a + b * 4 theta(1 - theta)),
#              a set so the share of multi-colony cultures matches the cohort,
#              and multi-colony depths drawn from the cohort's own. The lab
#              picks one colony per visible morphotype (README issue 1), so
#              mixed cultures are the ones that get depth. b = 0 is base.
#
# Fitted by the method of simulated moments to eight cohort summaries:
#   mean        mean NS fraction per culture
#   var_pat     variance of patient means (patients with >= 2 cultures)
#   var_within  pooled within-patient variance of the culture NS fraction
#   disc        within-culture discordance, P(two colonies differ), m >= 2
#   r_0_42 .. r_180_inf  patient-centred consecutive-culture correlation in
#               four gap bins (as Fig 4, colony frame, per class)
# weighted by patient-bootstrap SDs of the observed summaries. mu0 is set
# analytically from the mean given the total latent SD; sd_pat, sd_w and tau
# by grid search with common random numbers, then refined.
#
# The fit is the check: if one binomial-picking model cannot reproduce the
# within-culture discordance and the short-gap autocorrelation together, the
# misfit is reported, not hidden (analysis 1 found exactly that tension).
#
# Run from the repo root:  Rscript analyses/fit_power_calibration.R
# Output: results/power_recalibrated/calibration_*.csv
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
})
out_res <- "results/power_recalibrated"
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
R_GRID <- as.integer(Sys.getenv("CAL_R", "40"))     # simulated cohorts per objective
B_OBS  <- as.integer(Sys.getenv("CAL_B", "1000"))   # bootstrap for observed SDs

CLASSES <- c("Aminoglycoside", "Carbapenem", "Fluoroquinolone", "Cephalosporin")
CLASS_DRUGS <- list(
  Aminoglycoside  = c("amikacin", "tobramycin", "gentamicin"),
  Carbapenem      = c("meropenem", "imipenem"),
  Fluoroquinolone = c("ciprofloxacin", "levofloxacin"),
  Cephalosporin   = c("cefepime", "ceftazidime"))
SIR <- c("Susceptible", "Intermediate", "Resistant")
NS  <- c("Resistant", "Intermediate")
GAP_BREAKS <- c(0, 42, 90, 180, Inf)
STATS <- c("mean", "var_pat", "var_within", "disc", "r_0_42", "r_42_90", "r_90_180", "r_180_inf")

iso <- readRDS("data/pa_isolates.rds")
cult <- imap_dfr(CLASS_DRUGS, function(drugs, cl) {
  tests <- rowSums(sapply(drugs, function(x) iso[[x]] %in% SIR))
  ns    <- rowSums(sapply(drugs, function(x) iso[[x]] %in% NS)) > 0
  tibble(id = as.character(iso$individual_id),
         date = as.Date(iso$date_collected_from_subject), tests, ns) %>%
    filter(tests > 0) %>% group_by(id, date) %>%
    summarise(m = n(), k = sum(ns), .groups = "drop") %>% mutate(class = cl)
}) %>% arrange(class, id, date)


# #############################################################################
# DESIGN MATRICES: patients x culture index (NA-padded)
# #############################################################################

to_design <- function(d) {
  d <- d %>% group_by(id) %>% mutate(j = row_number(),
                                     gap = as.numeric(date - lag(date))) %>% ungroup()
  ids <- unique(d$id); P <- length(ids); J <- max(d$j)
  M <- matrix(NA_integer_, P, J); K <- M; G <- matrix(NA_real_, P, J)
  ix <- cbind(match(d$id, ids), d$j)
  M[ix] <- d$m; K[ix] <- d$k; G[ix] <- d$gap
  list(ids = ids, M = M, K = K, G = G, P = P, J = J,
       f_multi = mean(d$m >= 2), multi_depths = sort(d$m[d$m >= 2]))
}

# summaries from a K (NS colonies) and M (colonies) matrix pair
summ <- function(K, M, G) {
  Q <- K / M
  n_i <- rowSums(!is.na(Q)); pm <- rowMeans(Q, na.rm = TRUE)
  Qc <- Q - pm
  multi <- !is.na(M) & M >= 2
  disc <- 2 * K[multi] * (M[multi] - K[multi]) / (M[multi] * (M[multi] - 1))
  x <- Qc[, -ncol(Qc), drop = FALSE]; y <- Qc[, -1, drop = FALSE]; g <- G[, -1, drop = FALSE]
  ok <- !is.na(x) & !is.na(y)
  bin <- cut(g[ok], GAP_BREAKS)
  r <- tapply(seq_along(bin), bin, function(ix) cor(x[ok][ix], y[ok][ix]))
  c(mean = mean(Q, na.rm = TRUE),
    var_pat = var(pm[n_i >= 2]),
    var_within = sum(Qc^2, na.rm = TRUE) / sum(n_i[n_i >= 2] - 1),
    disc = mean(disc),
    setNames(as.numeric(r), c("r_0_42", "r_42_90", "r_90_180", "r_180_inf")))
}


# #############################################################################
# SIMULATOR on a design, with common random numbers
# #############################################################################

make_crn <- function(des, R) {
  list(zp = matrix(rnorm(des$P * R), des$P, R),
       zx = array(rnorm(des$P * des$J * R), c(des$P, des$J, R)),
       u  = array(runif(des$P * des$J * R), c(des$P, des$J, R)),
       u_multi = array(runif(des$P * des$J * R), c(des$P, des$J, R)),
       u_depth = array(runif(des$P * des$J * R), c(des$P, des$J, R)))
}

simulate_K <- function(des, mu0, sd_pat, sd_w, tau, crn, r, b = 0) {
  P <- des$P; J <- des$J
  mu <- mu0 + sd_pat * crn$zp[, r]
  X <- matrix(NA_real_, P, J)
  X[, 1] <- sd_w * crn$zx[, 1, r]
  for (j in seq_len(J)[-1]) {
    rho <- exp(-des$G[, j] / tau)
    X[, j] <- rho * X[, j - 1] + sd_w * sqrt(1 - rho^2) * crn$zx[, j, r]
  }
  TH <- plogis(mu + X)
  M <- des$M
  if (b > 0) {
    present <- !is.na(M)
    h <- 4 * TH[present] * (1 - TH[present])
    a <- uniroot(function(a) mean(plogis(a + b * h)) - des$f_multi, c(-60, 60))$root
    multi <- crn$u_multi[, , r][present] < plogis(a + b * h)
    md <- des$multi_depths
    M[present] <- ifelse(multi, md[pmin(length(md), 1 + floor(crn$u_depth[, , r][present] * length(md)))], 1L)
  }
  K <- qbinom(crn$u[, , r], M, TH)            # inverse-CDF: CRN across parameters
  K[is.na(M)] <- NA
  attr(K, "M") <- M
  K
}

# mu0 such that E[plogis(mu0 + N(0, s^2))] = target, by Gauss-Hermite
gh <- local({                      # Golub-Welsch, 40 nodes: E[f(Z)] = sum(w f(x))
  n <- 40; i <- seq_len(n - 1); Jm <- matrix(0, n, n)
  Jm[cbind(i, i + 1)] <- sqrt(i / 2); Jm[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(Jm, symmetric = TRUE)
  list(x = e$values * sqrt(2), w = e$vectors[1, ]^2)
})
solve_mu0 <- function(target, s) {
  uniroot(function(m0) sum(gh$w * plogis(m0 + s * gh$x)) - target, c(-15, 15))$root
}

sim_stats <- function(des, pars, crn, R = dim(crn$zx)[3]) {
  b <- if (is.null(pars$b)) 0 else pars$b
  ss <- vapply(seq_len(R), function(r) {
    K <- simulate_K(des, pars$mu0, pars$sd_pat, pars$sd_w, pars$tau, crn, r, b)
    summ(K, attr(K, "M"), des$G)
  }, numeric(length(STATS)))
  rowMeans(ss, na.rm = TRUE)
}


# #############################################################################
# FIT PER CLASS
# #############################################################################

GRID <- expand_grid(sd_pat = c(0.5, 1, 1.5, 2, 3, 4, 6, 8),
                    sd_w   = c(0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3),
                    tau    = c(14, 30, 60, 120, 240, 480, 1000, 2000))
B_GRID <- c(2, 5, 10, 20)

fit_class <- function(cl, version) {
  d <- cult %>% filter(class == cl)
  des <- to_design(d)
  obs <- summ(des$K, des$M, des$G)
  set.seed(1)                                   # same bootstrap for both versions
  bt <- t(replicate(B_OBS, {
    pick <- sample(des$P, des$P, replace = TRUE)
    summ(des$K[pick, , drop = FALSE], des$M[pick, , drop = FALSE], des$G[pick, , drop = FALSE])
  }))
  sdv <- apply(bt, 2, sd, na.rm = TRUE)
  set.seed(2)
  crn <- make_crn(des, R_GRID)
  obj <- function(sd_pat, sd_w, tau, b) {
    mu0 <- solve_mu0(obs[["mean"]], sqrt(sd_pat^2 + sd_w^2))
    s <- sim_stats(des, list(mu0 = mu0, sd_pat = sd_pat, sd_w = sd_w, tau = tau, b = b), crn)
    sum(((s - obs) / sdv)[STATS[-1]]^2)        # mean is matched by construction
  }
  t0 <- Sys.time()
  grid <- if (version == "base") GRID %>% mutate(b = 0) else expand_grid(GRID, b = B_GRID)
  g <- grid %>% mutate(obj = pmap_dbl(list(sd_pat, sd_w, tau, b), obj))
  best <- g %>% slice_min(obj, n = 1, with_ties = FALSE)
  if (version == "base") {
    ref <- optim(log(c(best$sd_pat, best$sd_w, best$tau)),
                 function(p) obj(exp(p[1]), exp(p[2]), exp(p[3]), 0),
                 method = "Nelder-Mead", control = list(maxit = 200, reltol = 1e-4))
    pars <- list(sd_pat = exp(ref$par[1]), sd_w = exp(ref$par[2]), tau = exp(ref$par[3]), b = 0)
  } else {
    ref <- optim(log(c(best$sd_pat, best$sd_w, best$tau, best$b)),
                 function(p) obj(exp(p[1]), exp(p[2]), exp(p[3]), exp(p[4])),
                 method = "Nelder-Mead", control = list(maxit = 300, reltol = 1e-4))
    pars <- list(sd_pat = exp(ref$par[1]), sd_w = exp(ref$par[2]), tau = exp(ref$par[3]),
                 b = exp(ref$par[4]))
  }
  pars$mu0 <- solve_mu0(obs[["mean"]], sqrt(pars$sd_pat^2 + pars$sd_w^2))
  set.seed(3)
  fit_s <- sim_stats(des, pars, make_crn(des, 200))   # fresh random numbers
  cat(sprintf("  %-9s %-16s grid best %.1f -> refined %.1f  (%.0f s)\n", version, cl, best$obj,
              ref$value, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  prof <- bind_rows(
    g %>% group_by(tau)  %>% slice_min(obj, n = 1, with_ties = FALSE) %>% ungroup() %>% mutate(profile = "tau"),
    g %>% group_by(sd_w) %>% slice_min(obj, n = 1, with_ties = FALSE) %>% ungroup() %>% mutate(profile = "sd_w"))
  list(class = cl, version = version, des = des, obs = obs, sd = sdv, pars = pars, fit = fit_s,
       obj = ref$value, grid = g %>% mutate(class = cl, version = version),
       prof = prof %>% mutate(class = cl, version = version))
}

fits <- unlist(map(c("base", "selection"), function(v) {
  f <- map(CLASSES, fit_class, version = v); names(f) <- paste(v, CLASSES); f
}), recursive = FALSE)


# #############################################################################
# OUTPUT
# #############################################################################

pars_tab <- map_dfr(fits, ~ tibble(version = .x$version, class = .x$class, mu0 = .x$pars$mu0,
                                   sd_pat = .x$pars$sd_pat, sd_w = .x$pars$sd_w, tau_days = .x$pars$tau,
                                   b = .x$pars$b, objective = .x$obj,
                                   patients = .x$des$P, cultures = sum(!is.na(.x$des$M))))
fit_tab <- map_dfr(fits, ~ tibble(version = .x$version, class = .x$class, stat = STATS, observed = .x$obs,
                                  boot_sd = .x$sd, fitted = .x$fit) %>%
                     mutate(z = (fitted - observed) / boot_sd))
write_csv(pars_tab, file.path(out_res, "calibration_parameters.csv"))
write_csv(fit_tab, file.path(out_res, "calibration_fit.csv"))
write_csv(map_dfr(fits, "prof"), file.path(out_res, "calibration_profiles.csv"))
write_csv(map_dfr(fits, "grid"), file.path(out_res, "calibration_grid.csv"))
saveRDS(map(fits, ~ .x[c("class", "version", "des", "obs", "sd", "pars")]),
        file.path(out_res, "calibration.rds"))

cat("\nParameters:\n"); print(as.data.frame(pars_tab), digits = 3, row.names = FALSE)
cat("\nFit (z = (fitted - observed) / bootstrap SD):\n")
print(as.data.frame(fit_tab %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
