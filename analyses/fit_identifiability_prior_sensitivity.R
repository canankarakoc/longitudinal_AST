# =============================================================================
# analyses/fit_identifiability_prior_sensitivity.R
#
# Two supplementary analyses for the AST measurement manuscript. Both refit
# models; neither touches the pipeline's cached fits.
#
#   1. Identifiability of the between-visit SD (Fig. S1). Re-runs AMR_pipeline.R
#      10.3d, whose results were recorded only in a code comment: simulate
#      outcomes on the observed drug x isolate design with a known sd_visit,
#      refit, see what comes back. Replicate 1 uses set.seed(1), exactly as in
#      10.3d, so it should reproduce the recorded values
#      (0 -> 0.267, 0.3 -> 0.215, 0.5 -> 0.401, 0.8 -> 0.666); replicates 2-5
#      add new simulated data sets so one draw of the data is not the result.
#
#   2. Prior sensitivity on the variance SDs (Fig. S5). Exponential(1) is not
#      innocuous when a component sits at zero. Refit the aminoglycoside
#      drug x isolate model and the four Table 1 models under:
#        exp1      Exponential(1) on every SD            (the published prior)
#        halft     half-Student-t(3, 0, 2.5) on every SD (brms default)
#        exp02     Exponential(0.2) on every SD          (weaker, mean 5)
#        exp5_vis  Exponential(5) on the visit SD only   (strong pull to zero)
#      exp1 is read from the existing cached fits rather than refitted.
#
# Run from the repo root:  Rscript analyses/fit_identifiability_prior_sensitivity.R
# Sequential by design (README: two brms jobs at once exhaust 24 GB).
# Results are appended row by row, so an interrupted run keeps what it did.
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(brms)
})
options(mc.cores = 4)

out_sim   <- "results/sensitivity/sd_visit_recovery.csv"
out_prior <- "results/sensitivity/prior_sensitivity.csv"
SIR <- c("Susceptible", "Intermediate", "Resistant")
NS  <- c("Resistant", "Intermediate")

iso <- readRDS("data/pa_isolates.rds")

# Same frame as AMR_pipeline.R 10.3
ag_long <- iso %>%
  mutate(visit_date = as.Date(date_collected_from_subject),
         isolate_id = paste0("iso", row_number()),
         visit_key  = paste(individual_id, visit_date, sep = "_")) %>%
  dplyr::select(individual_id, visit_key, isolate_id, amikacin, tobramycin, gentamicin) %>%
  pivot_longer(c(amikacin, tobramycin, gentamicin), names_to = "drug", values_to = "res") %>%
  filter(res %in% SIR) %>%
  mutate(nonsusc = as.integer(res %in% NS))
stopifnot(nrow(ag_long) == 4775)

q <- function(x) c(median = median(x), lo = unname(quantile(x, 0.025)),
                   hi = unname(quantile(x, 0.975)))
append_csv <- function(df, path) write_csv(df, path, append = file.exists(path))

# #############################################################################
# 1. IDENTIFIABILITY OF sd_visit
# #############################################################################

SD_TRUE <- c(0, 0.3, 0.5, 0.8)
REPS    <- 1:5
skel    <- ag_long %>% dplyr::select(individual_id, visit_key, isolate_id, drug)

simulate_on_skeleton <- function(sd_v, seed, sd_pat = 3.75, sd_iso = 2.03) {
  set.seed(seed)   # seed 1 = the 10.3d run
  u_p <- distinct(skel, individual_id) %>% mutate(u = rnorm(n(), 0, sd_pat))
  u_v <- distinct(skel, visit_key)     %>% mutate(u = rnorm(n(), 0, sd_v))
  u_i <- distinct(skel, isolate_id)    %>% mutate(u = rnorm(n(), 0, sd_iso))
  skel %>%
    left_join(u_p, by = "individual_id") %>% rename(up = u) %>%
    left_join(u_v, by = "visit_key")     %>% rename(uv = u) %>%
    left_join(u_i, by = "isolate_id")    %>% rename(ui = u) %>%
    mutate(nonsusc = rbinom(n(), 1, plogis(up + uv + ui)))
}

done_sim <- if (file.exists(out_sim)) read_csv(out_sim, show_col_types = FALSE) else NULL
template <- NULL
cat("\n=== 1. sd_visit recovery:", length(SD_TRUE), "values x", length(REPS), "replicates ===\n")
for (r in REPS) for (sv in SD_TRUE) {
  if (!is.null(done_sim) && any(done_sim$rep == r & done_sim$sd_visit_true == sv)) next
  d <- simulate_on_skeleton(sv, r)
  t0 <- Sys.time()
  fit <- if (is.null(template)) {
    brm(nonsusc ~ drug + (1 | individual_id) + (1 | visit_key) + (1 | isolate_id),
        data = d, family = bernoulli(),
        prior = c(prior(normal(0, 1.5), class = "b"), prior(normal(0, 2), class = "Intercept"),
                  prior(exponential(1), class = "sd")),
        chains = 2, cores = 2, iter = 2000, warmup = 1000,
        control = list(adapt_delta = 0.95), seed = 42, refresh = 0, silent = 2)
  } else {
    update(template, newdata = d, recompile = FALSE, chains = 2, cores = 2,
           iter = 2000, warmup = 1000, control = list(adapt_delta = 0.95),
           seed = 42, refresh = 0, silent = 2)
  }
  if (is.null(template)) template <- fit
  dr <- as_draws_df(fit)
  np <- nuts_params(fit)
  row <- tibble(rep = r, sd_visit_true = sv,
                est = median(dr$sd_visit_key__Intercept),
                lo  = quantile(dr$sd_visit_key__Intercept, 0.025),
                hi  = quantile(dr$sd_visit_key__Intercept, 0.975),
                sd_patient_est = median(dr$sd_individual_id__Intercept),
                sd_colony_est  = median(dr$sd_isolate_id__Intercept),
                rhat_max = max(rhat(fit), na.rm = TRUE),
                divergent = sum(np$Value[np$Parameter == "divergent__"]),
                minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  append_csv(row, out_sim)
  cat(sprintf("  rep %d  true %.1f -> %.3f [%.3f, %.3f]   (%.1f min)\n",
              r, sv, row$est, row$lo, row$hi, row$minutes))
  rm(fit, dr); gc(verbose = FALSE)
}
rm(template); gc(verbose = FALSE)

# #############################################################################
# 2. PRIOR SENSITIVITY
# #############################################################################

priors_sd <- function(label, visit_group) {
  switch(label,
    halft    = prior(student_t(3, 0, 2.5), class = "sd"),
    exp02    = prior(exponential(0.2), class = "sd"),
    exp5_vis = c(prior(exponential(1), class = "sd"),
                 set_prior("exponential(5)", class = "sd", group = visit_group)))
}
PRIORS <- c("halft", "exp02", "exp5_vis")

summarise_iso <- function(fit, prior_label) {
  dr <- as_draws_df(fit)
  v <- tibble(p = dr$sd_individual_id__Intercept^2, v = dr$sd_visit_key__Intercept^2,
              c = dr$sd_isolate_id__Intercept^2, r = pi^2 / 3) %>% mutate(t = p + v + c + r)
  np <- nuts_params(fit)
  bind_rows(
    tibble(component = c("sd_patient", "sd_visit", "sd_colony"),
           stat = "sd",
           !!!as.data.frame(rbind(q(sqrt(v$p)), q(sqrt(v$v)), q(sqrt(v$c))))),
    tibble(component = c("share_patient", "share_visit", "share_colony", "share_residual"),
           stat = "share",
           !!!as.data.frame(rbind(q(v$p / v$t), q(v$v / v$t), q(v$c / v$t), q(v$r / v$t))))) %>%
    mutate(model = "AG drug x isolate", class = "Aminoglycoside", prior = prior_label,
           rhat_max = max(rhat(fit), na.rm = TRUE),
           divergent = sum(np$Value[np$Parameter == "divergent__"]), .before = 1)
}

summarise_t1 <- function(fit, cls, prior_label) {
  dr <- as_draws_df(fit)
  p <- dr$sd_individual_id__Intercept^2
  v <- dr[["sd_individual_id:visit_number__Intercept"]]^2
  t <- p + v + pi^2 / 3
  np <- nuts_params(fit)
  bind_rows(
    tibble(component = c("sd_patient", "sd_visit"), stat = "sd",
           !!!as.data.frame(rbind(q(sqrt(p)), q(sqrt(v))))),
    tibble(component = c("share_patient", "share_visit", "share_colony"), stat = "share",
           !!!as.data.frame(rbind(q(p / t), q(v / t), q((pi^2 / 3) / t))))) %>%
    mutate(model = "Table 1", class = cls, prior = prior_label,
           rhat_max = max(rhat(fit), na.rm = TRUE),
           divergent = sum(np$Value[np$Parameter == "divergent__"]), .before = 1)
}

done_prior <- if (file.exists(out_prior)) read_csv(out_prior, show_col_types = FALSE) else NULL
have <- function(model, cls, pr) !is.null(done_prior) &&
  any(done_prior$model == model & done_prior$class == cls & done_prior$prior == pr)

cat("\n=== 2a. prior sensitivity, AG drug x isolate model ===\n")
if (!have("AG drug x isolate", "Aminoglycoside", "exp1")) {
  base <- readRDS("models/sensitivity/m_var_isolate.rds")
  append_csv(summarise_iso(base, "exp1"), out_prior); rm(base); gc(verbose = FALSE)
  cat("  exp1 (cached fit) summarised\n")
}
for (pr in PRIORS) {
  if (have("AG drug x isolate", "Aminoglycoside", pr)) next
  t0 <- Sys.time()
  fit <- brm(nonsusc ~ drug + (1 | individual_id) + (1 | visit_key) + (1 | isolate_id),
             data = ag_long, family = bernoulli(),
             prior = c(prior(normal(0, 1.5), class = "b"), prior(normal(0, 2), class = "Intercept"),
                       priors_sd(pr, "visit_key")),
             chains = 4, cores = 4, iter = 4000, warmup = 2000,
             control = list(adapt_delta = 0.95), seed = 42, refresh = 0, silent = 2)
  append_csv(summarise_iso(fit, pr), out_prior)
  cat(sprintf("  %-8s done (%.1f min)\n", pr, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  rm(fit); gc(verbose = FALSE)
}

cat("\n=== 2b. prior sensitivity, Table 1 models ===\n")
T1 <- list(Aminoglycoside  = list(drugs = c("amikacin", "tobramycin", "gentamicin"), tag = "ag"),
           Carbapenem      = list(drugs = c("meropenem", "imipenem"),               tag = "carb"),
           Fluoroquinolone = list(drugs = c("ciprofloxacin", "levofloxacin"),       tag = "fq"),
           Cephalosporin   = list(drugs = c("cefepime", "ceftazidime"),             tag = "ceph"))
for (cls in names(T1)) {
  dr_cols <- T1[[cls]]$drugs
  # same rows and outcome as AMR_common.R fit_variance_decomp()
  d <- iso %>% mutate(
    .tested = rowSums(sapply(dr_cols, function(x) iso[[x]] %in% SIR)),
    R = if_else(.tested == 0, NA_integer_,
                as.integer(rowSums(sapply(dr_cols, function(x) iso[[x]] %in% NS)) > 0))) %>%
    filter(!is.na(R))
  if (!have("Table 1", cls, "exp1")) {
    base <- readRDS(sprintf("models/noise/vardecomp_%s.rds", T1[[cls]]$tag))
    stopifnot(nobs(base) == nrow(d))
    append_csv(summarise_t1(base, cls, "exp1"), out_prior); rm(base); gc(verbose = FALSE)
  }
  for (pr in PRIORS) {
    if (have("Table 1", cls, pr)) next
    t0 <- Sys.time()
    fit <- brm(R ~ 1 + (1 | individual_id) + (1 | individual_id:visit_number),
               data = d, family = bernoulli(),
               prior = c(prior(normal(0, 2), class = "Intercept"),
                         priors_sd(pr, "individual_id:visit_number")),
               chains = 4, cores = 4, iter = 4000, warmup = 2000,
               control = list(adapt_delta = 0.95), seed = 42, refresh = 0, silent = 2)
    append_csv(summarise_t1(fit, cls, pr), out_prior)
    cat(sprintf("  %-15s %-8s done (%.1f min)\n", cls, pr,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    rm(fit); gc(verbose = FALSE)
  }
}
cat("\n=== DONE ===\n")
