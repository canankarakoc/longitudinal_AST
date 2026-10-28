# =============================================================================
# analyses/fit_pseudo_trajectories.R
#
# Sampling-induced pseudo-dynamics in the cohort: the same infection, read
# through routine one-colony sampling, can produce different apparent
# resistance histories depending only on which colony was picked.
#
# For every patient x class, each culture's tested colonies give a
# non-susceptible fraction p_j = k_j / m_j (colony frame: a colony is NS if any
# agent in the class reads R or I, as in analysis 1). A routine record picks
# one colony, so the culture reads NS with probability p_j, independently
# across cultures. Over R = 2000 picks per patient x class:
#   - apparent switches (S <-> NS between consecutive cultures)
#   - apparent history type: stable S | stable NS | acquired (one S -> NS) |
#     lost (one NS -> S) | fluctuating (>= 2 switches)
# Expected switches are also exact: sum_j p_j(1 - p_j+1) + (1 - p_j)p_j+1,
# split into switches fixed by the data (both cultures homogeneous) and
# switches that depend on the pick (either culture heterogeneous).
#
# Reference trajectory: all tested colonies, culture NS if ANY colony is NS
# (the most-resistant reading of the full record).
#
# THIS IS A LOWER BOUND. Heterogeneity is only visible at cultures with >= 2
# colonies; single-colony cultures enter as fixed. The "multi-colony only"
# variant keeps just the cultures whose heterogeneity is observed (longer gaps,
# fewer patients). Colony picking is purposive (README issue 1), so p_j is the
# NS share among the colonies the lab chose, not in the sputum.
#
# Patient bootstrap (B = 2000) on the cohort summaries.
#
# Run from the repo root:  Rscript analyses/fit_pseudo_trajectories.R
# Output: results/pseudo_trajectories/   figures/pseudo_trajectories/
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2)
})

out_res <- "results/pseudo_trajectories"
out_fig <- "figures/pseudo_trajectories"
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
R_PICK <- as.integer(Sys.getenv("PSEUDO_R", "2000"))
B      <- as.integer(Sys.getenv("PSEUDO_B", "2000"))

CLASSES   <- c("Aminoglycoside", "Carbapenem", "Fluoroquinolone", "Cephalosporin")
CLASS_COL <- c(Aminoglycoside = "#0072B2", Carbapenem = "#D55E00",
               Fluoroquinolone = "#009E73", Cephalosporin = "#CC79A7")
CLASS_DRUGS <- list(
  Aminoglycoside  = c("amikacin", "tobramycin", "gentamicin"),
  Carbapenem      = c("meropenem", "imipenem"),
  Fluoroquinolone = c("ciprofloxacin", "levofloxacin"),
  Cephalosporin   = c("cefepime", "ceftazidime"))
SIR <- c("Susceptible", "Intermediate", "Resistant")
NS  <- c("Resistant", "Intermediate")
TYPES <- c("stable S", "stable NS", "acquired", "lost", "fluctuating")

iso <- readRDS("data/pa_isolates.rds")

cult <- imap_dfr(CLASS_DRUGS, function(drugs, cl) {
  tests <- rowSums(sapply(drugs, function(x) iso[[x]] %in% SIR))
  ns    <- rowSums(sapply(drugs, function(x) iso[[x]] %in% NS)) > 0
  tibble(id = as.character(iso$individual_id),
         date = as.Date(iso$date_collected_from_subject), tests, ns) %>%
    filter(tests > 0) %>%
    group_by(id, date) %>%
    summarise(m = n(), k = sum(ns), .groups = "drop") %>%
    mutate(class = cl, p = k / m)
}) %>% arrange(class, id, date)


# #############################################################################
# PER PATIENT x CLASS
# #############################################################################

history_type <- function(x) {
  # x: 0/1 matrix, rows = picks, cols = cultures in time order
  sw <- rowSums(abs(x[, -1, drop = FALSE] - x[, -ncol(x), drop = FALSE]))
  first <- x[, 1]
  out <- ifelse(sw == 0, ifelse(first == 1, "stable NS", "stable S"),
         ifelse(sw == 1, ifelse(first == 0, "acquired", "lost"), "fluctuating"))
  list(sw = sw, type = factor(out, TYPES))
}

one_series <- function(d) {
  p <- d$p; J <- length(p)
  het <- p > 0 & p < 1
  # exact expected switches, split by whether either culture is heterogeneous
  e_pair <- p[-J] * (1 - p[-1]) + (1 - p[-J]) * p[-1]
  pick_dep <- het[-J] | het[-1]
  x <- matrix(rbinom(R_PICK * J, 1, rep(p, each = R_PICK)), nrow = R_PICK)
  h <- history_type(x)
  ref <- history_type(matrix(as.integer(d$k > 0), nrow = 1))
  tp <- prop.table(table(h$type))
  tibble(cultures = J, multi = sum(d$m >= 2), het = sum(het),
         e_switch = sum(e_pair), e_switch_fixed = sum(e_pair[!pick_dep]),
         e_switch_pick = sum(e_pair[pick_dep]),
         sd_switch = sd(h$sw),
         # two independent one-colony records of the same cultures
         p_differ2 = mean(abs(h$sw[seq(1, R_PICK, 2)] - h$sw[seq(2, R_PICK, 2)]) >= 2),
         p_switch_differs = 1 - sum(prop.table(table(h$sw))^2),
         ref_switch = ref$sw, ref_type = as.character(ref$type),
         p_matches_ref_type = mean(h$type == ref$type),
         modal_type = names(tp)[which.max(tp)], p_modal = max(tp),
         n_types_10pct = sum(tp >= 0.10),
         !!!setNames(as.list(as.numeric(tp)), paste0("P_", gsub(" ", "_", TYPES))))
}

run_variant <- function(dat, label) {
  dat %>% group_by(class, id) %>% filter(n() >= 2) %>%
    group_modify(~ one_series(.x)) %>% ungroup() %>% mutate(variant = label)
}

per_pat <- bind_rows(
  run_variant(cult, "All cultures"),
  run_variant(cult %>% filter(m >= 2), "Multi-colony cultures only"))
write_csv(per_pat, file.path(out_res, "pseudo_trajectories_per_patient.csv"))


# #############################################################################
# COHORT SUMMARIES + PATIENT BOOTSTRAP
# #############################################################################

cohort_stats <- function(pp) {
  pp %>% group_by(variant, class) %>%
    summarise(patients = n(), cultures = sum(cultures),
              pct_pat_het = mean(het > 0),
              e_switch_per_pat = mean(e_switch),
              ref_switch_per_pat = mean(ref_switch),
              share_switch_pick = sum(e_switch_pick) / sum(e_switch),
              pct_pat_switch_varies = mean(sd_switch > 0),
              mean_sd_switch = mean(sd_switch),
              mean_p_differ2 = mean(p_differ2),
              pct_pat_2plus_types = mean(n_types_10pct >= 2),
              pct_pat_type_unstable = mean(p_modal < 0.9),
              mean_p_matches_ref = mean(p_matches_ref_type),
              .groups = "drop")
}
point <- cohort_stats(per_pat)

boot <- map_dfr(seq_len(B), function(b) {
  per_pat %>% group_by(variant, class) %>%
    slice_sample(prop = 1, replace = TRUE) %>% ungroup() %>% cohort_stats()
})
ci <- boot %>% group_by(variant, class) %>%
  summarise(across(-c(patients, cultures),
                   list(lo = ~ quantile(.x, 0.025, na.rm = TRUE),
                        hi = ~ quantile(.x, 0.975, na.rm = TRUE))), .groups = "drop")
cohort <- point %>% left_join(ci, by = c("variant", "class"))
write_csv(cohort, file.path(out_res, "pseudo_trajectories_cohort.csv"))


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
pc_ <- function(v, lo, hi) sprintf("%.0f%% (%.0f–%.0f)", 100 * v, 100 * lo, 100 * hi)
nm_ <- function(v, lo, hi) sprintf("%.2f (%.2f–%.2f)", v, lo, hi)
tab <- cohort %>% transmute(
  variant, class, patients, cultures,
  `patients with a heterogeneous culture` = pc_(pct_pat_het, pct_pat_het_lo, pct_pat_het_hi),
  `switches/patient, one colony` = nm_(e_switch_per_pat, e_switch_per_pat_lo, e_switch_per_pat_hi),
  `switches/patient, all colonies any-NS` = nm_(ref_switch_per_pat, ref_switch_per_pat_lo, ref_switch_per_pat_hi),
  `share of switches set by the pick` = pc_(share_switch_pick, share_switch_pick_lo, share_switch_pick_hi),
  `patients whose switch count varies` = pc_(pct_pat_switch_varies, pct_pat_switch_varies_lo, pct_pat_switch_varies_hi),
  `SD of switch count across picks` = nm_(mean_sd_switch, mean_sd_switch_lo, mean_sd_switch_hi),
  `P(two records differ by >= 2 switches)` = pc_(mean_p_differ2, mean_p_differ2_lo, mean_p_differ2_hi),
  `patients with >= 2 histories at >= 10%` = pc_(pct_pat_2plus_types, pct_pat_2plus_types_lo, pct_pat_2plus_types_hi),
  `P(one-colony history = all-colony history)` = pc_(mean_p_matches_ref, mean_p_matches_ref_lo, mean_p_matches_ref_hi))

writeLines(c("# Pseudo-trajectories: one colony per culture",
             "", sprintf("%d random picks per patient x class; patient bootstrap B = %d.", R_PICK, B),
             "Colony frame, class-level NS (any agent R or I). Patients with >= 2 cultures in the class.",
             "History types: stable S, stable NS, acquired (one S->NS), lost (one NS->S), fluctuating (>= 2 switches).",
             "Lower bound: single-colony cultures enter as fixed.", "",
             md_tab(tab)),
           file.path(out_res, "pseudo_trajectories.md"))
print(as.data.frame(tab), row.names = FALSE)
cat("\nwrote", out_res, "and", out_fig, "\n")
