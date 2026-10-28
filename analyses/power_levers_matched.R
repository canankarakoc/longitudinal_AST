if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
# Power levers, MATCHED at 2x and off the ceiling.
#
# The Section 14.3 grid cannot answer this: it has m = 3,10 only, dt = 14/28/56/84
# only, n_pat = 30/90/200 only, and beta = 0.2 where visits and colonies both
# saturate at ~1.00 so the ranking is invisible. Five cells re-run here:
#
#   baseline  n=92, L=12, m=2 effective tests (= 1.7 colonies), dt=83 d
#   each lever doubled, one at a time, beta = 0.1
#
# m is EFFECTIVE independent tests, not drug-tests: 3 AG drugs on one colony is
# 1.4 independent tests (rho 0.57, deff 2.14), so 1.7 colonies ~ 2.4 -> use 2,
# and 3.4 colonies -> 4. The fold change is exactly 2x, which is the point.
suppressMessages({library(tidyverse); library(lme4)})
set.seed(42)
REPS <- 300

sim_cohort <- function(n_pat, L, rho, beta, m, sd_pat = 1.0, sd_inn = 1.0) {
  purrr::map_dfr(seq_len(n_pat), function(i) {
    mux <- rnorm(1, 0, sd_pat); muy <- rnorm(1, 0, sd_pat)
    x <- y <- numeric(L); x[1] <- rnorm(1); y[1] <- rnorm(1)
    for (t in 2:L) {
      x[t] <- rho * x[t-1] + rnorm(1, 0, sd_inn)
      y[t] <- rho * y[t-1] + beta * x[t-1] + rnorm(1, 0, sd_inn)
    }
    tibble(id = i, t = seq_len(L),
           kx = rbinom(L, m, plogis(mux + x)), ky = rbinom(L, m, plogis(muy + y)), m = m)
  })
}
detect_xlag <- function(d) {
  dd <- d %>% group_by(id) %>%
    mutate(py = ky/m, px = kx/m, lpy = lag(py), lpx = lag(px)) %>%
    ungroup() %>% filter(!is.na(lpy))
  fit <- try(suppressMessages(lmer(py ~ lpy + lpx + (1|id), data = dd)), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  s <- summary(fit)$coefficients
  if (!"lpx" %in% rownames(s)) return(NA_real_)
  2 * pnorm(-abs(s["lpx", "t value"]))
}
# latent AR(1) implied by sampling interval: observed r(dt) = a exp(-dt/tau),
# divided by reliability because the binomial step is applied separately
A_HAT <- 0.343; TAU <- 86; RELY <- 0.70
rho_at <- function(dt) min(max(A_HAT * exp(-dt/TAU) / RELY, 0.02), 0.95)

BETA <- 0.1
cells <- tribble(
  ~lever,             ~n_pat, ~L, ~m, ~dt,   ~cost,
  "baseline",             92, 12,  2, 83,    "as run today",
  "More colonies",        92, 12,  4, 83,    "1.7 -> 3.4 colonies: a bench protocol change",
  "Shorter interval",     92, 12,  2, 41,    "83 -> 41 days: fortnightly sputum, every patient",
  "Bigger cohort",       184, 12,  2, 83,    "92 -> 184 patients: a multi-centre study",
  "More visits",          92, 24,  2, 83,    "12 -> 24 visits: five more years of follow-up")

cat("all levers doubled; beta =", BETA, ";", REPS, "reps/cell\n\n")
res <- purrr::pmap_dfr(cells %>% select(-cost), function(lever, n_pat, L, m, dt) {
  r <- rho_at(dt)
  p <- replicate(REPS, detect_xlag(sim_cohort(n_pat, L, r, BETA, m)))
  pw <- mean(p < 0.05, na.rm = TRUE); nok <- sum(!is.na(p))
  cat(sprintf("  %-18s n=%3d L=%2d m=%d dt=%2d rho=%.2f -> power %.3f\n",
              lever, n_pat, L, m, dt, r, pw))
  tibble(lever, n_pat, L, m, dt, rho = r, power = pw,
         power_se = sqrt(pw*(1-pw)/nok))
}) %>% left_join(cells %>% select(lever, cost), by = "lever")

base <- res$power[res$lever == "baseline"]
res <- res %>% mutate(baseline = base, gain = power - base)
cat("\n=== MATCHED 2x LEVERS (beta = 0.1) ===\n")
print(res %>% select(lever, cost, power, power_se, gain) %>%
        mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(), row.names = FALSE)
write_csv(res, "results/design/lever_comparison_matched.csv")
