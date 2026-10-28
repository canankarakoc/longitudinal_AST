# =============================================================================
# analyses/fit_deep_downsampling.R
#
# Downsampling the deep-sampled sputa of Vanderwoude et al. 2024 (mBio; CC BY
# 4.0; supplement mbio.03109-23-s0008.xlsx, sheets S1-S4): four sputum
# samples, 75 isolates each, disc-diffusion zones in biological triplicate for
# six agents. Their P1 and P3 are our patients 7 and 38; P2 and P4 are not in
# the cohort but are the same kind of sample.
#
# Each isolate's call is CLSI S/I/R on its MEAN zone (the pipeline's thresholds,
# AMR_pipeline.R Section 7). The 75 calls are the reference population. Draw n
# colonies at random without replacement (D = 5000 draws per n) and ask, per
# estimand:
#   detect     P(at least one NS colony)              - minority detection
#   freq       mean |k/n - theta_NS|                   - frequency estimation
#   dominant   P(plurality call of n = plurality call of 75)
#   switch_*   P(two independent n-colony draws from the SAME sputum report
#              different calls) - an apparent change with no true change.
#              Under two report rules: plurality, and most-resistant colony.
# Two observation modes:
#   mean       the isolate's call from its mean zone (colony sampling only)
#   single     the call from ONE random biological replicate (colony sampling
#              plus disc-to-disc assay noise, as a routine single test reads)
# The triplicates also give assay reproducibility directly: P(two replicates
# of the same isolate get different calls), against P(two different isolates
# get different calls).
#
# Random picking is the design here: the reference has no morphotype labels,
# so morphotype-stratified picking cannot be simulated from these sheets.
#
# Run from the repo root:  Rscript analyses/fit_deep_downsampling.R
# Output: results/deep_downsampling/   figures/deep_downsampling/
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2)
  library(readxl)
})

out_res <- "results/deep_downsampling"
out_fig <- "figures/deep_downsampling"
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
D <- as.integer(Sys.getenv("DEEP_D", "5000"))
N_GRID <- c(1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30)

supp_path <- "mbio.03109-23-s0008.xlsx"
SAMPLES <- c(S1 = "P1 (cohort id 7)", S2 = "P2 (not in cohort)",
             S3 = "P3 (cohort id 38)", S4 = "P4 (not in cohort)")
DRUGS <- c(AK = "amikacin", TOB = "tobramycin", MEM = "meropenem",
           CAZ = "ceftazidime", CIP = "ciprofloxacin", TZP = "piperacillin-tazobactam")
# CLSI disc diffusion, ZOI mm: R if <= R, S if >= S, else I (as AMR_pipeline.R 7)
THR <- tibble(code = c("AK", "TOB", "MEM", "CAZ", "CIP", "TZP"),
              R = c(14, 12, 15, 14, 18, 14), S = c(17, 15, 19, 18, 25, 21))
SIR <- c("S", "I", "R")
call_zoi <- function(z, R, S) factor(ifelse(z <= R, "R", ifelse(z >= S, "S", "I")), SIR)


# #############################################################################
# DATA: isolate x replicate x drug
# #############################################################################

# The per-replicate columns are exactly half the sheet's AVG columns (checked
# below), i.e. radii; AVG, which AMR_pipeline.R Section 7 reads, is the
# diameter. Replicates are therefore doubled to diameters before applying the
# CLSI diameter breakpoints. Whether AVG follows the CLSI convention (diameter
# including the 6 mm disc) is a question for the source's Methods; it is
# inherited unchanged from Section 7.
read_sheet <- function(sh) {
  read_excel(supp_path, sheet = sh, .name_repair = "minimal") %>%
    select(Isolate, `Biological Replicate`, all_of(names(DRUGS)),
           all_of(paste(names(DRUGS), "AVG"))) %>%
    tidyr::fill(Isolate)
}
avg_chk <- map_dfr(names(SAMPLES), function(sh) {
  d <- read_sheet(sh)
  map_dfr(names(DRUGS), function(cd) {
    d %>% group_by(Isolate) %>%
      summarise(avg = first(na.omit(.data[[paste(cd, "AVG")]])),
                twice_mean = 2 * mean(.data[[cd]], na.rm = TRUE), .groups = "drop")
  })
})
stopifnot(max(abs(avg_chk$avg - avg_chk$twice_mean), na.rm = TRUE) < 1e-6)
cat("  sheet AVG = 2 x mean of replicate columns, every isolate and drug\n")

reps <- imap_dfr(SAMPLES, function(lab, sh) {
  read_sheet(sh) %>%
    select(Isolate, `Biological Replicate`, all_of(names(DRUGS))) %>%
    pivot_longer(all_of(names(DRUGS)), names_to = "code", values_to = "radius") %>%
    transmute(sample = lab, isolate = Isolate, rep = `Biological Replicate`, code,
              zoi = 2 * radius)
}) %>% filter(!is.na(zoi)) %>%
  left_join(THR, by = "code") %>%
  mutate(drug = DRUGS[code], call_rep = call_zoi(zoi, R, S))

iso_calls <- reps %>% group_by(sample, isolate, drug, code, R, S) %>%
  summarise(zoi_mean = mean(zoi), n_rep = n(), .groups = "drop") %>%
  mutate(call = call_zoi(zoi_mean, R, S))

# Check against the pipeline's split (mean of triplicates = the sheet's AVG)
old <- read_csv("results/deep/deep_sample_RIS_split.csv", show_col_types = FALSE)
chk <- iso_calls %>% filter(sample %in% SAMPLES[c("S1", "S3")]) %>%
  group_by(sample, drug) %>%
  summarise(pctR = round(100 * mean(call == "R")), pctI = round(100 * mean(call == "I")),
            .groups = "drop")
cat("isolates per sample x drug:", paste(range(count(iso_calls, sample, drug)$n), collapse = "-"), "\n")
cat("check vs results/deep (P1, P3): pctR", paste(sort(chk$pctR), collapse = " "), "| old",
    paste(sort(old$pctR), collapse = " "), "\n")
stopifnot(identical(sort(chk$pctR), sort(old$pctR)), identical(sort(chk$pctI), sort(old$pctI)))

truth <- iso_calls %>% group_by(sample, drug) %>%
  summarise(N = n(), pS = mean(call == "S"), pI = mean(call == "I"), pR = mean(call == "R"),
            theta_NS = pI + pR, .groups = "drop") %>%
  mutate(plurality = SIR[max.col(cbind(pS, pI, pR), ties.method = "last")],
         heterogeneous_SIR = pmax(pS, pI, pR) < 0.97,
         heterogeneous_NS  = theta_NS > 0.03 & theta_NS < 0.97)
write_csv(truth, file.path(out_res, "deep_reference_populations.csv"))


# #############################################################################
# ASSAY vs COLONY reproducibility
# #############################################################################

pair_disc <- function(x) {          # P(two distinct draws differ), x = factor calls
  n <- length(x); if (n < 2) return(NA_real_)
  t <- table(x); 1 - sum(t * (t - 1)) / (n * (n - 1))
}
assay <- reps %>% group_by(sample, drug, isolate) %>%
  summarise(d = pair_disc(call_rep), .groups = "drop") %>%
  group_by(sample, drug) %>% summarise(assay_rep_disc = mean(d, na.rm = TRUE), .groups = "drop")
colony <- iso_calls %>% group_by(sample, drug) %>%
  summarise(colony_disc_mean = pair_disc(call), .groups = "drop")
colony_single <- reps %>% group_by(sample, drug) %>%
  summarise(colony_disc_single = {
    # two different isolates, one random replicate each
    iso_ids <- unique(isolate)
    a <- sample(iso_ids, 20000, TRUE); b <- sample(iso_ids, 20000, TRUE); keep <- a != b
    ca <- call_rep[match(paste(a, sample(1:3, 20000, TRUE)), paste(isolate, rep))]
    cb <- call_rep[match(paste(b, sample(1:3, 20000, TRUE)), paste(isolate, rep))]
    mean((ca != cb)[keep], na.rm = TRUE)
  }, .groups = "drop")
noise <- truth %>% select(sample, drug, pS, pI, pR) %>%
  left_join(assay, by = c("sample", "drug")) %>%
  left_join(colony, by = c("sample", "drug")) %>%
  left_join(colony_single, by = c("sample", "drug"))
write_csv(noise, file.path(out_res, "assay_vs_colony_discordance.csv"))


# #############################################################################
# DOWNSAMPLING
# #############################################################################

# integer codes: S = 1, I = 2, R = 3
plurality_call <- function(m) {           # rows = draws, cols = picked colonies
  cnt <- cbind(rowSums(m == 1), rowSums(m == 2), rowSums(m == 3))
  # ties broken toward the more resistant call (the cautious report)
  max.col(cnt + matrix(c(0, 1e-6, 2e-6), nrow(cnt), 3, byrow = TRUE))
}

one_cell <- function(calls_mean, calls_rep_mat, truth_row) {
  N <- length(calls_mean)
  pop_mean <- as.integer(calls_mean)
  plur75 <- match(truth_row$plurality, SIR)
  map_dfr(N_GRID, function(n) {
    map_dfr(c("mean", "single"), function(mode) {
      pick <- function() {
        idx <- t(replicate(D, sample.int(N, n)))
        if (mode == "mean") matrix(pop_mean[idx], nrow = D)
        else matrix(calls_rep_mat[cbind(as.vector(idx), sample.int(3, length(idx), TRUE))], nrow = D)
      }
      A <- pick(); Bm <- pick()
      nsA <- rowSums(A >= 2); nsB <- rowSums(Bm >= 2)
      plA <- plurality_call(A); plB <- plurality_call(Bm)
      wA <- apply(A, 1, max); wB <- apply(Bm, 1, max)
      tibble(n = n, mode = mode,
             detect = mean(nsA > 0),
             freq_abs_err = mean(abs(nsA / n - truth_row$theta_NS)),
             dominant_correct = mean(plA == plur75),
             switch_plurality = mean(plA != plB),
             switch_worst = mean(wA != wB),
             switch_NS_any = mean((nsA > 0) != (nsB > 0)))
    })
  })
}

cells <- truth %>% select(sample, drug)
down <- pmap_dfr(cells, function(sample, drug) {
  ic <- iso_calls %>% filter(.data$sample == !!sample, .data$drug == !!drug) %>% arrange(isolate)
  rm <- reps %>% filter(.data$sample == !!sample, .data$drug == !!drug) %>%
    select(isolate, rep, call_rep) %>%
    pivot_wider(names_from = rep, values_from = call_rep) %>%
    arrange(isolate)
  stopifnot(identical(rm$isolate, ic$isolate))
  rep_mat <- sapply(rm[, -1], as.integer)
  # an isolate with a missing replicate falls back to its mean call there
  rep_mat[is.na(rep_mat)] <- as.integer(ic$call)[row(rep_mat)[is.na(rep_mat)]]
  tr <- truth %>% filter(.data$sample == !!sample, .data$drug == !!drug)
  one_cell(ic$call, rep_mat, tr) %>% mutate(sample = sample, drug = drug, .before = 1)
}) %>% left_join(truth, by = c("sample", "drug"))
write_csv(down, file.path(out_res, "deep_downsampling.csv"))

# exact check: hypergeometric detection, mean mode
down_chk <- down %>% filter(mode == "mean") %>%
  mutate(K = round(theta_NS * N),
         exact = 1 - choose(N - K, n) / choose(N, n))
cat(sprintf("max |simulated - exact| detection probability: %.3f\n",
            max(abs(down_chk$detect - down_chk$exact))))

# colonies needed, per heterogeneous cell
needed <- down %>% filter(heterogeneous_SIR) %>%
  group_by(sample, drug, mode) %>%
  summarise(n_detect_95 = if (any(theta_NS > 0 & detect >= 0.95)) min(n[detect >= 0.95]) else NA,
            n_dominant_90 = if (any(dominant_correct >= 0.90)) min(n[dominant_correct >= 0.90]) else NA,
            n_switch_plur_10 = if (any(switch_plurality <= 0.10)) min(n[switch_plurality <= 0.10]) else NA,
            n_switch_worst_10 = if (any(switch_worst <= 0.10)) min(n[switch_worst <= 0.10]) else NA,
            .groups = "drop")
write_csv(needed, file.path(out_res, "deep_colonies_needed.csv"))


# #############################################################################
# FIGURE: drawn by analyses/figures_new.R from the CSVs written here, at print
# size (analyses/fig_style.R). This script no longer draws.
# #############################################################################

het <- down %>% filter(heterogeneous_SIR)            # the report's heterogeneous cells

# #############################################################################
# REPORT
# #############################################################################

md_tab <- function(df) {
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
pc_ <- function(x) sprintf("%.0f%%", 100 * x)

writeLines(c(
  "# Downsampling the Vanderwoude et al. 2024 deep samples",
  "", sprintf("75 isolates per sputum, triplicate discs; D = %d draws per n. CLSI calls on mean zone.", D),
  "", "## Reference populations", "",
  md_tab(truth %>% transmute(sample, drug, S = pc_(pS), I = pc_(pI), R = pc_(pR),
                             NS = pc_(theta_NS), plurality,
                             heterogeneous = ifelse(heterogeneous_SIR, "yes", ""))),
  "", "## Assay vs colony discordance (every cell)", "",
  md_tab(noise %>% transmute(sample, drug,
                             `same isolate, two discs` = pc_(assay_rep_disc),
                             `two isolates, mean zone` = pc_(colony_disc_mean),
                             `two isolates, one disc each` = pc_(colony_disc_single))),
  "", "## Heterogeneous cells: estimands by colonies tested (colony sampling only / + single-disc noise)", "",
  md_tab(het %>% filter(n %in% c(1, 2, 4, 8, 15)) %>%
           select(sample, drug, n, mode, detect, freq_abs_err, dominant_correct,
                  switch_plurality, switch_worst) %>%
           pivot_wider(names_from = mode,
                       values_from = c(detect, freq_abs_err, dominant_correct,
                                       switch_plurality, switch_worst)) %>%
           transmute(sample, drug, n,
                     detect = pc_(detect_mean),
                     `|freq err|` = sprintf("%.2f", freq_abs_err_mean),
                     `dominant correct` = sprintf("%s / %s", pc_(dominant_correct_mean), pc_(dominant_correct_single)),
                     `switch, plurality` = sprintf("%s / %s", pc_(switch_plurality_mean), pc_(switch_plurality_single)),
                     `switch, most-resistant` = sprintf("%s / %s", pc_(switch_worst_mean), pc_(switch_worst_single)))),
  "", "## Colonies needed (heterogeneous cells)", "",
  "detect >= 95%; dominant call right >= 90%; apparent-switch probability <= 10%. NA = not reached by 30.", "",
  md_tab(needed %>% mutate(across(everything(), ~ ifelse(is.na(.x), "NA", as.character(.x)))))),
  file.path(out_res, "deep_downsampling.md"))

cat("\nwrote", out_res, "and", out_fig, "\n")
