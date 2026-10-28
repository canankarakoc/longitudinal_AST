# =============================================================================
# analyses/fit_design_by_estimand.R
#
# Which sampling design for which question? No design is best for every
# estimand, so designs are compared estimand by estimand:
#
#   detect     report NS when an NS subpopulation is present (sensitivity, and
#              the false-positive rate at theta = 0 once assay error is on)
#   frequency  estimate the NS fraction theta (mean absolute error)
#   dominant   classify the majority population (theta > 0.5) correctly
#   change     P(two visits with the SAME theta report different calls) - the
#              false-change rate - and P(calls differ) for a true change
#   clinical   response to treatment: NOT evaluable from these data
#
# Designs (cost in AST runs):
#   random-n        n colonies at random (n = 1, 4, 8); cost n
#   stratified-k    k colonies from each VISIBLE morphotype; cost k x visible
#                   frequency: naive (pooled share of tested colonies) or
#                   weighted by plate abundance of each morphotype
#   pooled-10       10 random colonies pooled into one AST; cost 1
#   whole sputum    idealised assay reading theta itself; cost 1
# A pooled or whole-sample readout is ONE result whose mapping from population
# composition is unknown. It is bracketed by two observation rules: NS if the
# NS share >= 1% (sensitive, like the proportion method) or >= 50% (majority).
# These are assumptions, not validated readouts.
#
# Report rules: colony designs call NS for detection if ANY tested colony is NS
# (most-resistant colony), and the dominant population by the share of NS
# colonies (random: > 0.5; stratified: weighted share > 0.5; ties at random).
#
# Population structures (two morphotypes, A and B; B has abundance pi_B):
#   aligned     resistance confined to one morphotype: theta_B = 1,
#               theta_A = 0, pi_B = theta. Stratified picking excels here.
#   unaligned   resistance spread evenly: pi_B = 0.3, theta_A = theta_B =
#               theta. Morphotype carries no information.
# A morphotype is visible on the plate with prob 1 - (1 - abundance)^50.
# The cohort's own alignment (between- vs within-morphotype discordance inside
# a culture, colony frame) is computed first to place routine practice between
# these brackets.
#
# Assay error: each colony test (and each pooled / whole readout) flips with
# probability eps = 0 or 0.03 (illustrative; disc replicates of one isolate
# disagree 0-27% by cell in the deep samples, mostly ~0).
#
# Also an empirical check on the five heterogeneous Vanderwoude cells
# (analysis 3): random-n vs pooled-10 vs whole-sample on real populations.
#
# Run from the repo root:  Rscript analyses/fit_design_by_estimand.R
# Output: results/design_estimand/   figures/design_estimand/
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2)
})

out_res <- "results/design_estimand"
out_fig <- "figures/design_estimand"
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
D     <- as.integer(Sys.getenv("DESIGN_D", "20000"))
B     <- as.integer(Sys.getenv("DESIGN_B", "1000"))
PLATE <- 50
THETA <- c(0, seq(0.025, 0.975, by = 0.05), 1)   # 0.5 excluded: dominant undefined

CLASSES <- c("Aminoglycoside", "Carbapenem", "Fluoroquinolone", "Cephalosporin")
CLASS_DRUGS <- list(
  Aminoglycoside  = c("amikacin", "tobramycin", "gentamicin"),
  Carbapenem      = c("meropenem", "imipenem"),
  Fluoroquinolone = c("ciprofloxacin", "levofloxacin"),
  Cephalosporin   = c("cefepime", "ceftazidime"))
SIR <- c("Susceptible", "Intermediate", "Resistant")
NS  <- c("Resistant", "Intermediate")


# #############################################################################
# 1  COHORT ANCHOR: is resistance aligned with morphotype inside a culture?
#
# Within each culture, pairs of colonies of DIFFERENT morphotype vs pairs of the
# SAME morphotype: P(NS calls differ). Aligned structure -> within ~ 0, between
# high; unaligned -> equal. alignment = 1 - within / between. Each culture
# weighted once per pair type; patient bootstrap.
# #############################################################################

iso <- readRDS("data/pa_isolates.rds")
col <- imap_dfr(CLASS_DRUGS, function(drugs, cl) {
  tibble(id = as.character(iso$individual_id), date = iso$date_collected_from_subject,
         morph = iso$pseudomonas_aeruginosa_phenotype,
         tests = rowSums(sapply(drugs, function(x) iso[[x]] %in% SIR)),
         ns = rowSums(sapply(drugs, function(x) iso[[x]] %in% NS)) > 0, class = cl) %>%
    filter(tests > 0, !is.na(morph)) %>% mutate(colony = row_number())
})

pairs <- col %>%
  inner_join(col %>% select(class, id, date, colony, morph, ns),
             by = c("class", "id", "date"), suffix = c("", ".j"),
             relationship = "many-to-many") %>%
  filter(colony < colony.j) %>%
  mutate(type = ifelse(morph == morph.j, "within", "between"), differ = ns != ns.j) %>%
  group_by(class, id, date, type) %>% summarise(d = mean(differ), .groups = "drop")

align_stats <- function(p) {
  p %>% group_by(class, type) %>% summarise(d = mean(d), .groups = "drop") %>%
    pivot_wider(names_from = type, values_from = d) %>%
    mutate(alignment = 1 - within / between)
}
align_pt <- align_stats(pairs) %>%
  left_join(pairs %>% count(class, type) %>% pivot_wider(names_from = type, values_from = n,
                                                         names_prefix = "cultures_"), by = "class")
pp <- split(pairs, pairs$id); ids <- names(pp)
align_bt <- map_dfr(seq_len(B), function(b) align_stats(bind_rows(pp[sample(ids, length(ids), TRUE)])))
align <- align_pt %>% left_join(
  align_bt %>% group_by(class) %>%
    summarise(across(c(between, within, alignment),
                     list(lo = ~ quantile(.x, 0.025, na.rm = TRUE),
                          hi = ~ quantile(.x, 0.975, na.rm = TRUE))), .groups = "drop"),
  by = "class")
write_csv(align, file.path(out_res, "cohort_morphotype_alignment.csv"))
cat("\nCohort morphotype alignment (colony frame):\n"); print(as.data.frame(align), digits = 2)


# #############################################################################
# 2  SIMULATION
# #############################################################################

flip <- function(k, n, eps) {          # observed NS count after per-test error
  if (eps == 0) return(k)
  rbinom(length(k), k, 1 - eps) + rbinom(length(k), n - k, eps)
}
flip1 <- function(x, eps) if (eps == 0) x else xor(x, runif(length(x)) < eps)

structure_pars <- function(st, theta) {
  if (st == "aligned") list(piB = theta, thA = 0, thB = 1)
  else list(piB = 0.3, thA = theta, thB = theta)
}

# One visit per row; returns, per design, the detection call, the frequency
# estimate and the dominant call.
simulate_visits <- function(st, theta, eps, nd = D) {
  p <- structure_pars(st, theta)
  out <- list()
  draw_random <- function(n) {
    nb <- rbinom(nd, n, p$piB)
    k <- rbinom(nd, nb, p$thB) + rbinom(nd, n - nb, p$thA)
    flip(k, n, eps)
  }
  for (n in c(1, 4, 8)) {
    k <- draw_random(n)
    dom <- ifelse(k / n == 0.5, runif(nd) < 0.5, k / n > 0.5)
    out[[paste0("random-", n)]] <- tibble(cost = n, detect = k > 0, freq = k / n, dominant = dom)
  }
  visA <- runif(nd) < 1 - p$piB^PLATE           # A has abundance 1 - piB
  visB <- runif(nd) < 1 - (1 - p$piB)^PLATE
  none <- !visA & !visB; visA[none] <- p$piB < 0.5; visB[none] <- p$piB >= 0.5
  for (k_per in c(1, 2)) {
    kA <- k_per * visA; kB <- k_per * visB
    nsA <- flip(rbinom(nd, kA, p$thA), kA, eps)
    nsB <- flip(rbinom(nd, kB, p$thB), kB, eps)
    wA <- ifelse(visA, 1 - p$piB, 0); wB <- ifelse(visB, p$piB, 0)
    fw <- (wA * ifelse(kA > 0, nsA / pmax(kA, 1), 0) + wB * ifelse(kB > 0, nsB / pmax(kB, 1), 0)) /
      (wA + wB)
    fn <- (nsA + nsB) / (kA + kB)
    nm <- paste0("stratified-", k_per)
    domw <- ifelse(abs(fw - 0.5) < 1e-9, runif(nd) < 0.5, fw > 0.5)   # ties at random, as random-n
    out[[nm]] <- tibble(cost = kA + kB, detect = nsA + nsB > 0, freq = fw, dominant = domw,
                        freq_naive = fn)
  }
  kp <- rbinom(nd, 10, theta)
  for (tau in c(0.01, 0.5)) {
    lab <- if (tau < 0.1) "1%" else "50%"
    r <- flip1(kp / 10 >= tau, eps)
    out[[paste0("pooled-10, NS if >= ", lab)]] <- tibble(cost = 1, detect = r, freq = NA_real_, dominant = r)
    w <- flip1(rep(theta >= tau, nd), eps)
    out[[paste0("whole sputum, NS if >= ", lab)]] <- tibble(cost = 1, detect = w, freq = NA_real_, dominant = w)
  }
  bind_rows(out, .id = "design")
}

grid <- expand_grid(structure = c("aligned", "unaligned"), theta = THETA, eps = c(0, 0.03))
sim <- pmap_dfr(grid, function(structure, theta, eps) {
  a <- simulate_visits(structure, theta, eps)
  b <- simulate_visits(structure, theta, eps)        # an independent second visit
  a %>% mutate(detect2 = b$detect, dominant2 = b$dominant) %>%
    group_by(design) %>%
    summarise(cost = mean(cost),
              p_report_NS = mean(detect),
              freq_mae = mean(abs(freq - theta)),
              freq_bias = mean(freq - theta),
              freq_naive_mae = if (all(is.na(freq_naive))) NA_real_ else mean(abs(freq_naive - theta), na.rm = TRUE),
              dominant_correct = mean(dominant == (theta > 0.5)),
              false_change_detect = mean(detect != detect2),
              false_change_dominant = mean(dominant != dominant2), .groups = "drop") %>%
    mutate(structure = structure, theta = theta, eps = eps, .before = 1)
})
write_csv(sim, file.path(out_res, "design_by_estimand_grid.csv"))

# Power for a true change theta1 -> theta2, against the false-change rate at theta1
CHANGES <- tribble(~theta1, ~theta2, 0, 0.225, 0.075, 0.375, 0.225, 0.725, 0.725, 0.975)
change <- pmap_dfr(expand_grid(structure = c("aligned", "unaligned"), eps = c(0, 0.03),
                               CHANGES), function(structure, eps, theta1, theta2) {
  a <- simulate_visits(structure, theta1, eps); b <- simulate_visits(structure, theta2, eps)
  a0 <- simulate_visits(structure, theta1, eps)
  a %>% mutate(d2 = b$detect, m2 = b$dominant, d0 = a0$detect, m0 = a0$dominant,
               f2 = b$freq, f0 = a0$freq) %>%
    group_by(design) %>%
    summarise(power_detect = mean(detect != d2), false_detect = mean(detect != d0),
              power_dominant = mean(dominant != m2), false_dominant = mean(dominant != m0),
              # frequency designs: call a change when the estimate moves >= 0.25
              power_freq = mean(abs(f2 - freq) >= 0.25), false_freq = mean(abs(f0 - freq) >= 0.25),
              .groups = "drop") %>%
    mutate(structure = structure, eps = eps, theta1 = theta1, theta2 = theta2, .before = 1)
})
write_csv(change, file.path(out_res, "design_change_power.csv"))


# #############################################################################
# 3  EMPIRICAL CHECK: the heterogeneous Vanderwoude cells (analysis 3)
# Each isolate's NS call from its mean zone; random-n vs pooled-10 vs whole.
# #############################################################################

deep_pop <- read_csv("results/deep_downsampling/deep_reference_populations.csv", show_col_types = FALSE) %>%
  filter(heterogeneous_SIR)
deep <- pmap_dfr(deep_pop %>% select(sample, drug, N, theta_NS), function(sample, drug, N, theta_NS) {
  K <- round(theta_NS * N); pop <- c(rep(1, K), rep(0, N - K))
  dr <- function(n) matrix(pop[as.vector(replicate(2 * D / 10, sample.int(N, n)))], ncol = n, byrow = TRUE)
  res <- list()
  for (n in c(1, 4, 8)) {
    m <- dr(n); k <- rowSums(m); h <- length(k) / 2
    dom <- ifelse(k / n == 0.5, runif(length(k)) < 0.5, k / n > 0.5)
    res[[paste0("random-", n)]] <- tibble(detect = k > 0, dominant = dom, freq = k / n, half = rep(1:2, each = h))
  }
  m <- dr(10); k <- rowSums(m); h <- length(k) / 2
  for (tau in c(0.01, 0.5)) {
    lab <- if (tau < 0.1) "1%" else "50%"
    res[[paste0("pooled-10, NS if >= ", lab)]] <- tibble(detect = k / 10 >= tau, dominant = k / 10 >= tau,
                                                        freq = NA_real_, half = rep(1:2, each = h))
    res[[paste0("whole sputum, NS if >= ", lab)]] <- tibble(detect = rep(theta_NS >= tau, length(k)),
                                                           dominant = rep(theta_NS >= tau, length(k)),
                                                           freq = NA_real_, half = rep(1:2, each = h))
  }
  bind_rows(res, .id = "design") %>% group_by(design) %>%
    summarise(p_report_NS = mean(detect),
              freq_mae = mean(abs(freq - theta_NS)),
              dominant_correct = mean(dominant == (theta_NS > 0.5)),
              false_change_detect = mean(detect[half == 1] != detect[half == 2]),
              false_change_dominant = mean(dominant[half == 1] != dominant[half == 2]),
              .groups = "drop") %>%
    mutate(sample = sample, drug = drug, theta_NS = theta_NS, .before = 1)
})
write_csv(deep, file.path(out_res, "deep_empirical_designs.csv"))


# #############################################################################
# FIGURE: drawn by analyses/figures_new.R from the CSVs written here, at print
# size (analyses/fig_style.R). This script no longer draws.
# #############################################################################

DESIGNS <- c("random-1", "random-4", "random-8", "stratified-1", "stratified-2",
             "pooled-10, NS if >= 1%", "pooled-10, NS if >= 50%",
             "whole sputum, NS if >= 1%", "whole sputum, NS if >= 50%")   # report row order

# #############################################################################
# REPORT
# #############################################################################

md_tab <- function(df) {
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
pc_ <- function(x) ifelse(is.na(x), "–", sprintf("%.0f%%", 100 * x))
f2_ <- function(x) ifelse(is.na(x), "–", sprintf("%.2f", x))

snap <- function(th, e) sim %>% filter(abs(theta - th) < 1e-9, eps == e) %>%
  mutate(design = factor(design, DESIGNS)) %>% arrange(structure, design) %>%
  transmute(structure, design, cost = sprintf("%.1f", cost),
            `P(report NS)` = pc_(p_report_NS), `freq MAE` = f2_(freq_mae),
            `naive freq MAE` = f2_(freq_naive_mae), `dominant correct` = pc_(dominant_correct),
            `false change (detect)` = pc_(false_change_detect),
            `false change (dominant)` = pc_(false_change_dominant))

writeLines(c(
  "# Sampling design by estimand",
  "", sprintf("D = %d simulated visits per cell; plate of %d colonies for morphotype visibility.", D, PLATE),
  "Pooled / whole-sputum readouts use ASSUMED observation rules (NS if NS share >= 1% or >= 50%).",
  "", "## Cohort anchor: morphotype alignment inside cultures (colony frame, patient bootstrap)", "",
  md_tab(align %>% transmute(class, cultures_between, cultures_within,
                             `P(differ), different morphotype` = sprintf("%.2f (%.2f–%.2f)", between, between_lo, between_hi),
                             `P(differ), same morphotype` = sprintf("%.2f (%.2f–%.2f)", within, within_lo, within_hi),
                             alignment = sprintf("%.2f (%.2f–%.2f)", alignment, alignment_lo, alignment_hi))),
  "", "alignment = 1 - within/between: 1 = resistance follows morphotype exactly, 0 = morphotype uninformative.",
  "", "## Snapshot: theta = 0.025 (rare minority), no assay error", "", md_tab(snap(0.025, 0)),
  "", "## Snapshot: theta = 0.225, no assay error", "", md_tab(snap(0.225, 0)),
  "", "## Snapshot: theta = 0.225, 3% per-test error", "", md_tab(snap(0.225, 0.03)),
  "", "## Snapshot: theta = 0 (false positives), 3% per-test error", "", md_tab(snap(0, 0.03)),
  "", "## True change vs false change (no assay error)", "",
  md_tab(change %>% filter(eps == 0) %>% mutate(design = factor(design, DESIGNS)) %>%
           arrange(structure, theta1, design) %>%
           transmute(structure, change = sprintf("%.3f -> %.3f", theta1, theta2), design,
                     `P(differ), detect call` = pc_(power_detect), `false, detect` = pc_(false_detect),
                     `P(differ), dominant call` = pc_(power_dominant), `false, dominant` = pc_(false_dominant),
                     `P(|d freq| >= 0.25)` = pc_(power_freq), `false, freq` = pc_(false_freq))),
  "", "## Empirical: heterogeneous Vanderwoude cells", "",
  md_tab(deep %>% mutate(design = factor(design, DESIGNS)) %>% arrange(sample, drug, design) %>%
           transmute(sample, drug, NS = pc_(theta_NS), design, `P(report NS)` = pc_(p_report_NS),
                     `freq MAE` = f2_(freq_mae), `dominant correct` = pc_(dominant_correct),
                     `false change (detect)` = pc_(false_change_detect),
                     `false change (dominant)` = pc_(false_change_dominant)))),
  file.path(out_res, "design_by_estimand.md"))

cat("\nwrote", out_res, "and", out_fig, "\n")
