# =============================================================================
# analyses/fit_reliability_within_patient.R
#
# Splits the Fig 3 reliability into a cross-sectional and a within-patient
# part. Fig 3 inverts the variance of visit-level frequencies POOLED over all
# visits of all patients:
#     V_theta = (V_obs - mu(1 - mu) E[1/m]) / (1 - E[1/m])
# so V_theta carries the stable between-patient differences. "Reliability 0.8
# at four to five colonies" therefore answers "can one culture tell high- from
# low-resistance patients?", not "can one culture resolve a change around the
# patient's own level?" - which is what the dynamics argument needs.
#
# Same moment model, same visits (m >= 2), same noise term E = E[theta(1-theta)]
# and same design effect as Fig 3, but the signal is split:
#   V_tot  pooled (reproduces Fig 3)                          -> R_cs(m)
#   V_w    visit deviations from the patient's own mean:
#            E[SS_i] = (n_i - 1)(V_w + mean_j s_ij),  s_ij = noise at visit j
#          so V_w = (SS_w - sum_i (n_i - 1) mean_j s_ij) / sum_i (n_i - 1)
#                                                             -> R_wp(m)
#   V_b    V_tot - V_w
#   V_dth  variance of the TRUE change between consecutive (m >= 2) visits:
#            E[dq^2] = E[dtheta^2] + s_j + s_j+1              -> R_chg(m)
# R_wp(m) = V_w / (V_w + E / (m m_eff)) is the reliability of one visit's
# deviation from its patient's mean; R_chg(m) = V_dth / (V_dth + 2E/(m m_eff))
# is the reliability of an observed visit-to-visit change.
#
# Two noise estimators for s_ij:
#   pooledE  E / m_ij with E = mu(1-mu) - V_tot, exactly the Fig 3 assumption
#   pervisit q(1-q)/(m-1), unbiased for theta(1-theta)/m under binomial picking
#            and free of the pooled E. One-per-morphotype picking (README issue
#            1) inflates within-culture disagreement, so this one biases V_w low.
#
# Two frames, as in Fig 3: drug-test frame (class_deff, the published Fig 3)
# and colony frame. In the drug-test frame correlated agents on one colony are
# treated as binomial draws when V is estimated, then credited m_eff per colony
# in R(m) - inherited from Fig 3, flagged rather than fixed here.
#
# Also re-derives the reliability used to deattenuate the 14.2 autocorrelation
# for the 14.3 power simulation (mean pooled R = 0.74). The autocorrelation is
# computed on patient-centred values, so it is attenuated by R_wp, not R_cs.
#
# Patient-level bootstrap (B = 2000); m_eff held at its full-data value.
#
# Run from the repo root:  Rscript analyses/fit_reliability_within_patient.R
# Output: results/reliability_split/   figures/reliability_split/
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2)
})

out_res <- "results/reliability_split"
out_fig <- "figures/reliability_split"
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
B <- as.integer(Sys.getenv("RELSPLIT_B", "2000"))
M_GRID <- 1:50

CLASSES   <- c("Aminoglycoside", "Carbapenem", "Fluoroquinolone", "Cephalosporin")
CLASS_KEY <- c(Aminoglycoside = "ag", Carbapenem = "carb",
               Fluoroquinolone = "fq", Cephalosporin = "ceph")
CLASS_COL <- c(Aminoglycoside = "#0072B2", Carbapenem = "#D55E00",
               Fluoroquinolone = "#009E73", Cephalosporin = "#CC79A7")
CLASS_DRUGS <- list(
  Aminoglycoside  = c("amikacin", "tobramycin", "gentamicin"),
  Carbapenem      = c("meropenem", "imipenem"),
  Fluoroquinolone = c("ciprofloxacin", "levofloxacin"),
  Cephalosporin   = c("cefepime", "ceftazidime"))
SIR <- c("Susceptible", "Intermediate", "Resistant")
NS  <- c("Resistant", "Intermediate")

iso <- readRDS("data/pa_isolates.rds")
vis <- readRDS("data/pa_visit_full.rds")
deff <- read_csv("results/manuscript/fig3_design_effect.csv", show_col_types = FALSE)
m_eff <- setNames(deff$m_eff_per_colony, deff$class)


# #############################################################################
# INPUTS: one row per culture per class, in each frame
# #############################################################################

frame_dt <- map_dfr(CLASSES, function(cl) {
  key <- CLASS_KEY[[cl]]
  tibble(id = vis$individual_id, date = as.Date(vis$date_collected_from_subject),
         q = vis[[paste0(key, "_prop_resistant_conservative")]],
         m = vis[[paste0(key, "_n_tested")]], class = cl)
}) %>% filter(!is.na(q), !is.na(m), m >= 1)

frame_col <- imap_dfr(CLASS_DRUGS, function(drugs, cl) {
  tests <- rowSums(sapply(drugs, function(x) iso[[x]] %in% SIR))
  ns    <- rowSums(sapply(drugs, function(x) iso[[x]] %in% NS)) > 0
  tibble(id = iso$individual_id, date = as.Date(iso$date_collected_from_subject),
         tests = tests, ns = ns) %>%
    filter(tests > 0) %>%
    group_by(id, date) %>%
    summarise(m = n(), k = sum(ns), .groups = "drop") %>%
    transmute(id, date, q = k / m, m, class = cl)
})

frames <- bind_rows(drugtest = frame_dt, colony = frame_col, .id = "frame")


# #############################################################################
# ESTIMATOR
# #############################################################################

# d: one class x frame, visits with m >= 2 (Fig 3's sample). Returns the
# variance components; works on a bootstrap draw where `id` is the draw's id.
components <- function(d) {
  d <- d %>% filter(m >= 2)
  mu <- mean(d$q); Einv <- mean(1 / d$m)
  V_tot <- (var(d$q) - mu * (1 - mu) * Einv) / (1 - Einv)
  E <- mu * (1 - mu) - V_tot
  d <- d %>% mutate(s_pooled = E / m, s_visit = q * (1 - q) / (m - 1))

  w <- d %>% group_by(id) %>% filter(n() >= 2) %>%
    summarise(n = n(), ss = sum((q - mean(q))^2),
              s_pooled = mean(s_pooled), s_visit = mean(s_visit), .groups = "drop")
  df_w <- sum(w$n - 1)
  V_w_pooled <- (sum(w$ss) - sum((w$n - 1) * w$s_pooled)) / df_w
  V_w_visit  <- (sum(w$ss) - sum((w$n - 1) * w$s_visit)) / df_w

  ch <- d %>% arrange(id, date) %>% group_by(id) %>%
    mutate(dq = q - lag(q), s2p = s_pooled + lag(s_pooled), s2v = s_visit + lag(s_visit),
           gap = as.numeric(date - lag(date))) %>%
    ungroup() %>% filter(!is.na(dq))

  tibble(visits = nrow(d), patients = n_distinct(d$id),
         patients_w = nrow(w), df_w = df_w, pairs = nrow(ch),
         median_gap = median(ch$gap),
         mu = mu, E = E, V_tot = V_tot,
         V_w_pooledE = V_w_pooled, V_w_pervisit = V_w_visit,
         V_b_pooledE = V_tot - V_w_pooled,
         V_dth_pooledE = mean(ch$dq^2) - mean(ch$s2p),
         V_dth_pervisit = mean(ch$dq^2) - mean(ch$s2v))
}

# V <= 0 (a moment estimate below zero) means no recoverable signal: R = 0
rel     <- function(V, E, m, me) if (V > 0) V / (V + E / (m * me)) else 0 * m
rel_chg <- function(V, E, m, me) if (V > 0) V / (V + 2 * E / (m * me)) else 0 * m
need <- function(f, target = 0.8) {
  r <- f(M_GRID); if (any(r >= target)) min(M_GRID[r >= target]) else Inf
}

# reliability summaries from one set of components
summarise_rel <- function(cp, me) {
  with(cp, tibble(
    share_within = V_w_pooledE / V_tot,
    R_cs_1 = rel(V_tot, E, 1, me),        R_cs_2 = rel(V_tot, E, 2, me),
    R_wp_1 = rel(V_w_pooledE, E, 1, me),  R_wp_2 = rel(V_w_pooledE, E, 2, me),
    R_wp_4 = rel(V_w_pooledE, E, 4, me),
    R_wpv_1 = rel(V_w_pervisit, E, 1, me), R_wpv_2 = rel(V_w_pervisit, E, 2, me),
    R_chg_1 = rel_chg(V_dth_pooledE, E, 1, me), R_chg_4 = rel_chg(V_dth_pooledE, E, 4, me),
    need_cs  = need(function(m) rel(V_tot, E, m, me)),
    need_wp  = need(function(m) rel(V_w_pooledE, E, m, me)),
    need_wpv = need(function(m) rel(V_w_pervisit, E, m, me)),
    need_chg = need(function(m) rel_chg(V_dth_pooledE, E, m, me))))
}

me_for <- function(frame, cl) if (frame == "colony") 1 else m_eff[[cl]]


# #############################################################################
# POINT ESTIMATES + PATIENT BOOTSTRAP
# #############################################################################

cells <- expand_grid(frame = c("drugtest", "colony"), class = CLASSES)

res <- pmap_dfr(cells, function(frame, class) {
  d  <- frames %>% filter(.data$frame == !!frame, .data$class == !!class)
  me <- me_for(frame, class)
  cp <- components(d)
  pt <- bind_cols(cp, summarise_rel(cp, me))
  # bootstrap: index once by patient for speed
  by_id <- split(d, d$id); ids <- names(by_id)
  bt <- map_dfr(seq_len(B), function(b) {
    pick <- sample(ids, length(ids), replace = TRUE)
    db <- bind_rows(map2(by_id[pick], seq_along(pick), ~ mutate(.x, id = .y)))
    cpb <- components(db)
    bind_cols(cpb %>% select(V_tot, V_w_pooledE, V_w_pervisit, V_b_pooledE, V_dth_pooledE),
              summarise_rel(cpb, me))
  })
  ci <- bt %>% summarise(across(everything(),
           list(lo = ~ quantile(.x, 0.025, na.rm = TRUE, type = 1),
                hi = ~ quantile(.x, 0.975, na.rm = TRUE, type = 1))))
  cat(sprintf("  %-9s %-16s done\n", frame, class))
  bind_cols(tibble(frame = frame, class = class, m_eff = me), pt, ci)
})

write_csv(res, file.path(out_res, "reliability_split.csv"))

# Check: R_cs must reproduce Fig 3 (class_deff and colony rows)
f3 <- read_csv("results/manuscript/fig3_reliability_at_1to4_colonies.csv", show_col_types = FALSE) %>%
  filter(method %in% c("class_deff", "colony")) %>%
  mutate(frame = ifelse(method == "colony", "colony", "drugtest"))
chk <- inner_join(res %>% select(frame, class, R_cs_1, R_cs_2), f3, by = c("frame", "class"))
stopifnot(nrow(chk) == 8, max(abs(chk$R_cs_1 - chk$R_1)) < 1e-8,
          max(abs(chk$R_cs_2 - chk$R_2)) < 1e-8)
cat("  R_cs reproduces Fig 3 (class_deff and colony) exactly\n")


# #############################################################################
# CURVES for the figure
# #############################################################################

curves <- pmap_dfr(cells, function(frame, class) {
  r <- res %>% filter(.data$frame == !!frame, .data$class == !!class)
  tibble(frame = frame, class = class, colonies = 1:12) %>%
    mutate(`Cross-sectional (pooled, Fig 3)` = rel(r$V_tot, r$E, colonies, r$m_eff),
           `Within-patient deviation` = rel(r$V_w_pooledE, r$E, colonies, r$m_eff),
           `Within-patient, lower CI` = rel(r$V_w_pooledE_lo, r$E, colonies, r$m_eff),
           `Within-patient, upper CI` = rel(r$V_w_pooledE_hi, r$E, colonies, r$m_eff),
           `Visit-to-visit change` = rel_chg(r$V_dth_pooledE, r$E, colonies, r$m_eff))
})
write_csv(curves, file.path(out_res, "reliability_split_curves.csv"))

# #############################################################################
# FIGURE: drawn by analyses/figures_new.R from the CSVs written here, at print
# size (analyses/fig_style.R). This script no longer draws.
# #############################################################################

# #############################################################################
# DEATTENUATION used by the 14.3 power simulation
#
# Pipeline: rho_latent(dt) = r_obs(dt) / mean(reliability_now), mean = 0.74
# (as_run pooled reliability at 1/mean(m), averaged over classes). The 14.2
# pairs are patient-centred and include m = 1 visits, so each pair is
# attenuated by sqrt(R_j R_j+1) with R the WITHIN-patient reliability at that
# visit's own depth. Drug-test frame, binomial in drug-tests (the assumption
# under which V was estimated).
# #############################################################################

fit4 <- read_csv("results/manuscript/fig4_autocorrelation_fit.csv", show_col_types = FALSE)
cat("\nFig 4 fit:\n"); print(as.data.frame(fit4))

att <- map_dfr(CLASSES, function(cl) {
  r <- res %>% filter(frame == "drugtest", class == cl)
  frame_dt %>% filter(class == cl) %>% arrange(id, date) %>% group_by(id) %>%
    mutate(R_cs = rel(r$V_tot, r$E, m, 1), R_wp = rel(r$V_w_pooledE, r$E, m, 1),
           att_cs = sqrt(R_cs * lag(R_cs)), att_wp = sqrt(R_wp * lag(R_wp))) %>%
    ungroup() %>% filter(!is.na(att_cs)) %>%
    summarise(class = cl, pairs = n(), att_cs = mean(att_cs), att_wp = mean(att_wp))
})
att_all <- att %>% summarise(class = "Pooled (as in 14.2)",
                             att_cs = weighted.mean(att_cs, pairs),
                             att_wp = weighted.mean(att_wp, pairs), pairs = sum(pairs))
att <- bind_rows(att, att_all) %>% mutate(pipeline_Rbar = 0.74)
write_csv(att, file.path(out_res, "deattenuation_factor.csv"))


# #############################################################################
# CONSISTENCY: can V_w carry the observed short-gap autocorrelation?
#
# Under the moment model noise is independent between visits, so the
# patient-centred correlation of consecutive visits is at most
# mean sqrt(R_wp,j R_wp,j+1) (latent rho = 1). Patient centring biases the
# observed r downwards (~ -1/(n_i - 1)), so r_obs above that ceiling cannot be
# a centring artefact. Gaps <= 42 days (Fig 4: r ~ 0.25-0.32 there).
# Each frame on all its visits (m >= 1, as in 14.2) and on m >= 2 pairs only
# (the visits V_w was estimated from).
# #############################################################################

short_gap_ac <- function(d) {
  cp <- components(d)
  map_dfr(c(1, 2), function(mmin) {
    p <- d %>% filter(m >= mmin) %>% arrange(id, date) %>% group_by(id) %>%
      mutate(xc = q - mean(q), lxc = lag(xc), gap = as.numeric(date - lag(date)),
             Rv = rel(cp$V_w_pooledE, cp$E, m, 1), ceil = sqrt(Rv * lag(Rv))) %>%
      ungroup() %>% filter(!is.na(lxc), gap <= 42)
    tibble(visits = ifelse(mmin == 1, "all", "m>=2"), pairs = nrow(p),
           r_obs = cor(p$xc, p$lxc), ceiling_wp = mean(p$ceil),
           excess = r_obs - ceiling_wp)
  })
}

consist <- pmap_dfr(cells, function(frame, class) {
  d <- frames %>% filter(.data$frame == !!frame, .data$class == !!class)
  by_id <- split(d, d$id); ids <- names(by_id)
  bt <- map_dfr(seq_len(B), function(b) {
    pick <- sample(ids, length(ids), replace = TRUE)
    short_gap_ac(bind_rows(map2(by_id[pick], seq_along(pick), ~ mutate(.x, id = .y))))
  })
  ci <- bt %>% group_by(visits) %>%
    summarise(across(c(r_obs, ceiling_wp, excess),
                     list(lo = ~ quantile(.x, 0.025, na.rm = TRUE),
                          hi = ~ quantile(.x, 0.975, na.rm = TRUE))),
              p_excess_le0 = mean(excess <= 0, na.rm = TRUE), .groups = "drop")
  short_gap_ac(d) %>% left_join(ci, by = "visits") %>%
    mutate(frame = frame, class = class, .before = 1)
})
write_csv(consist, file.path(out_res, "consistency_short_gap_ac.csv"))


# #############################################################################
# REPORT
# #############################################################################

f2 <- function(x) sprintf("%.2f", x)
ci <- function(v) sprintf("%s (%s–%s)", f2(res[[v]]), f2(res[[paste0(v, "_lo")]]),
                          f2(res[[paste0(v, "_hi")]]))
nd <- function(x) ifelse(is.infinite(x), ">50", as.character(x))
ndci <- function(v) sprintf("%s (%s–%s)", nd(res[[v]]), nd(res[[paste0(v, "_lo")]]),
                            nd(res[[paste0(v, "_hi")]]))

tab <- tibble(
  frame = res$frame, class = res$class,
  `visits / patients (m≥2)` = sprintf("%d / %d", res$visits, res$patients),
  V_tot = f2(res$V_tot), `V_within (pooled E)` = ci("V_w_pooledE"),
  `V_within (per-visit)` = ci("V_w_pervisit"), `within share` = ci("share_within"),
  `R_cs(1)` = f2(res$R_cs_1), `R_wp(1)` = ci("R_wp_1"), `R_wp(2)` = ci("R_wp_2"),
  `R_chg(1)` = ci("R_chg_1"),
  `colonies for 0.8: cs` = nd(res$need_cs), wp = ndci("need_wp"),
  `wp per-visit` = ndci("need_wpv"), change = ndci("need_chg"))

md_tab <- function(df) {
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
writeLines(c("# Reliability split: cross-sectional vs within-patient",
             "", sprintf("Patient bootstrap B = %d; 95%% percentile intervals.", B), "",
             md_tab(tab), "", "## Deattenuation factor for the 14.2 autocorrelation", "",
             md_tab(att %>% mutate(across(where(is.numeric) & !pairs, f2))), "",
             "## Consistency: short-gap (<= 42 d) within-patient autocorrelation vs the V_w ceiling", "",
             md_tab(consist %>% transmute(frame, class, visits, pairs,
               r_obs = sprintf("%s (%s–%s)", f2(r_obs), f2(r_obs_lo), f2(r_obs_hi)),
               ceiling = sprintf("%s (%s–%s)", f2(ceiling_wp), f2(ceiling_wp_lo), f2(ceiling_wp_hi)),
               excess = sprintf("%s (%s–%s)", f2(excess), f2(excess_lo), f2(excess_hi)),
               `boot P(excess <= 0)` = f2(p_excess_le0)))),
           file.path(out_res, "reliability_split.md"))

cat("\n"); print(as.data.frame(tab), row.names = FALSE)
cat("\n"); print(as.data.frame(att), digits = 3, row.names = FALSE)
cat("\n"); print(as.data.frame(consist), digits = 3, row.names = FALSE)
cat("\nwrote", out_res, "and", out_fig, "\n")
