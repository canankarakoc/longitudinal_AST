# =============================================================================
# analyses/fit_sir_disagreement_matrix.R
#
# The direct empirical result: what would a different colony pick from the SAME
# culture have done to the reported S / I / R call?
#
# Builds on AMR_sir_ordinal.R Section 2 (same nine drugs, same visit_key, same
# ordered colony pairs), which reported only the diagonal complement - "a
# second colony disagrees 29% after S, 34% after R, 74% after I" - pooled over
# drugs and weighted by ordered pairs. Here:
#
#   1. the full 3 x 3 matrix P(partner call | index call)
#   2. two weightings
#        pairs    every ordered colony pair counts once (reproduces 29/34/74;
#                 a culture with n colonies counts n(n-1) times, so deep
#                 cultures dominate, and depth tracks clinical concern)
#        culture  every culture x drug counts once: the matrix is the chance
#                 that two colonies picked at random from that culture give
#                 (i, j). This is the "different pick" estimand and is primary.
#   3. overall discordance (unordered), its two-step part (S <-> R), and the
#      binary flip the pipeline uses (S vs R + I)
#   4. by class and by drug
#   5. same vs different morphotype (mucoid / non-mucoid). The lab picks one
#      colony per morphotype (README issue 1), so different-morphotype pairs
#      are over-represented relative to random picking; the split shows how
#      much of the disagreement that selection carries.
#   6. comparator: a colony at this culture vs a colony at the patient's NEXT
#      culture (same drug), weighted per consecutive-culture pair. If
#      within-culture discordance is close to next-culture discordance, a
#      single-colony "change" between visits is mostly a different pick.
#      Reported for all index cultures and, like-for-like, for index cultures
#      with >= 2 colonies tested (the cultures within-culture pairs come from).
#
# Patient-cluster bootstrap (B = 1000) on pre-aggregated weights.
#
# Run from the repo root:  Rscript analyses/fit_sir_disagreement_matrix.R
# Output: results/sir_matrix/   figures/sir_matrix/
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2)
})

out_res <- "results/sir_matrix"
out_fig <- "figures/sir_matrix"
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
B <- as.integer(Sys.getenv("SIRMAT_B", "1000"))

SIR <- c("Susceptible", "Intermediate", "Resistant")
SIR_AB <- c(Susceptible = "S", Intermediate = "I", Resistant = "R")
CLASS_ORDER <- c("Aminoglycoside", "Carbapenem", "Fluoroquinolone", "Cephalosporin",
                 "Penicillin", "Monobactam")
CLASS_COL <- c(Aminoglycoside = "#0072B2", Carbapenem = "#D55E00",
               Fluoroquinolone = "#009E73", Cephalosporin = "#CC79A7",
               Penicillin = "#E69F00", Monobactam = "#56B4E9")

iso <- readRDS("data/pa_isolates.rds")


# #############################################################################
# LONG FRAME - identical drug screen and keys to AMR_sir_ordinal.R Section 0
# #############################################################################

pa_panel <- c("amikacin", "gentamicin", "tobramycin", "imipenem", "meropenem",
              "cefepime", "ceftazidime", "ciprofloxacin", "levofloxacin",
              "piperacillin_tazobactam", "aztreonam", "colistin")
drugs_use <- keep(pa_panel, function(d) {
  v <- iso[[d]]; v <- v[v %in% SIR]; length(v) >= 500 && sum(v == "Intermediate") >= 50
})
drug_class <- c(
  amikacin = "Aminoglycoside", gentamicin = "Aminoglycoside", tobramycin = "Aminoglycoside",
  meropenem = "Carbapenem", imipenem = "Carbapenem",
  cefepime = "Cephalosporin", ceftazidime = "Cephalosporin",
  levofloxacin = "Fluoroquinolone", ciprofloxacin = "Fluoroquinolone",
  piperacillin_tazobactam = "Penicillin", aztreonam = "Monobactam", colistin = "Polymyxin")
cat("Drugs:", paste(drugs_use, collapse = ", "), "\n")

long <- iso %>%
  mutate(isolate_id = sprintf("iso%04d", row_number()),
         visit_key  = paste(individual_id, visit_number, sep = "_"),
         phenotype  = pseudomonas_aeruginosa_phenotype) %>%
  select(individual_id, visit_key, visit_number, isolate_id, phenotype, all_of(drugs_use)) %>%
  pivot_longer(all_of(drugs_use), names_to = "drug", values_to = "sir") %>%
  filter(sir %in% SIR) %>%
  mutate(class = unname(drug_class[drug]), individual_id = as.character(individual_id))


# #############################################################################
# PAIRS
# #############################################################################

# within culture: ordered pairs of distinct colonies, same visit and drug
within <- long %>%
  inner_join(long %>% select(visit_key, drug, isolate_id, sir, phenotype),
             by = c("visit_key", "drug"), suffix = c(".i", ".j"),
             relationship = "many-to-many") %>%
  filter(isolate_id.i != isolate_id.j) %>%
  group_by(visit_key, drug) %>% mutate(w_culture = 1 / n()) %>% ungroup() %>%
  transmute(individual_id, unit = paste(visit_key, drug), drug, class,
            index = sir.i, partner = sir.j, w_pairs = 1, w_culture,
            morph = case_when(is.na(phenotype.i) | is.na(phenotype.j) ~ "unknown",
                              phenotype.i == phenotype.j ~ "same morphotype",
                              TRUE ~ "different morphotype"),
            comparison = "Same culture")

# next culture: every colony at visit v x every colony at the patient's next
# visit with a call for the drug; each consecutive pair of cultures counts once
vis_seq <- long %>% distinct(individual_id, drug, visit_number) %>%
  arrange(individual_id, drug, visit_number) %>%
  group_by(individual_id, drug) %>% mutate(next_visit = lead(visit_number)) %>%
  ungroup() %>% filter(!is.na(next_visit))
nextv <- long %>%
  inner_join(vis_seq, by = c("individual_id", "drug", "visit_number")) %>%
  inner_join(long %>% select(individual_id, drug, visit_number, sir2 = sir),
             by = c("individual_id", "drug", "next_visit" = "visit_number"),
             relationship = "many-to-many") %>%
  group_by(individual_id, drug, visit_number) %>% mutate(w_culture = 1 / n()) %>% ungroup() %>%
  transmute(individual_id, unit = paste(individual_id, drug, visit_number), drug, class,
            index = sir, partner = sir2, w_pairs = 1, w_culture,
            morph = "n/a", comparison = "Next culture")

# Like-for-like: within-culture pairs exist only where >= 2 colonies were
# tested, and the lab tests more colonies when it sees more morphotypes. Restrict
# the comparator to index cultures that ALSO had >= 2 colonies for the drug.
multi_idx <- long %>% count(individual_id, drug, visit_number, name = "n_col") %>%
  filter(n_col >= 2) %>% transmute(unit = paste(individual_id, drug, visit_number))
nextv_m <- nextv %>% semi_join(multi_idx, by = "unit") %>%
  mutate(comparison = "Next culture, multi-colony index")

pairs <- bind_rows(within, nextv, nextv_m) %>%
  mutate(index = factor(index, SIR), partner = factor(partner, SIR))
cat(sprintf("within-culture ordered pairs %d (%d cultures x drug, %d patients); next-culture pairs %d\n",
            nrow(within), n_distinct(within$unit), n_distinct(within$individual_id), nrow(nextv)))

# pre-aggregate: patient x stratum x index x partner, both weightings. Every
# statistic below is a ratio of sums of these, so the bootstrap only resamples
# patients and re-sums.
agg <- function(p, ...) {
  p %>% group_by(individual_id, ..., index, partner) %>%
    summarise(w_pairs = sum(w_pairs), w_culture = sum(w_culture), .groups = "drop")
}
strata <- list(
  pooled       = agg(pairs, comparison) %>% mutate(stratum = "All drugs"),
  class        = agg(pairs, comparison, class) %>% rename(stratum = class),
  drug         = agg(pairs, comparison, drug) %>% rename(stratum = drug),
  morph        = agg(pairs %>% filter(comparison == "Same culture", morph != "unknown"),
                     comparison, morph) %>% rename(stratum = morph),
  morph_class  = agg(pairs %>% filter(comparison == "Same culture", morph != "unknown"),
                     comparison, morph, class) %>%
                 mutate(stratum = paste(class, morph, sep = " / ")) %>% select(-class, -morph))
strata <- imap_dfr(strata, ~ mutate(.x, level = .y))


# #############################################################################
# STATISTICS
# #############################################################################

stats_from <- function(a) {
  long_w <- a %>% pivot_longer(c(w_pairs, w_culture), names_to = "weighting", values_to = "w") %>%
    mutate(weighting = sub("w_", "", weighting))
  cells <- long_w %>% group_by(level, comparison, stratum, weighting, index, partner) %>%
    summarise(w = sum(w), .groups = "drop") %>%
    complete(nesting(level, comparison, stratum, weighting), index, partner, fill = list(w = 0))
  rows <- cells %>% group_by(level, comparison, stratum, weighting, index) %>%
    mutate(p = w / sum(w)) %>% ungroup() %>%
    transmute(level, comparison, stratum, weighting,
              stat = paste0("P(", SIR_AB[as.character(partner)], "|", SIR_AB[as.character(index)], ")"),
              value = p)
  overall <- cells %>% group_by(level, comparison, stratum, weighting) %>%
    summarise(discord   = sum(w[index != partner]) / sum(w),
              two_step  = sum(w[(index == "Susceptible" & partner == "Resistant") |
                                (index == "Resistant" & partner == "Susceptible")]) / sum(w),
              flip_S_vs_NS = sum(w[(index == "Susceptible") != (partner == "Susceptible")]) / sum(w),
              flip_R_vs_SI = sum(w[(index == "Resistant") != (partner == "Resistant")]) / sum(w),
              .groups = "drop") %>%
    pivot_longer(c(discord, two_step, flip_S_vs_NS, flip_R_vs_SI), names_to = "stat")
  bind_rows(rows, overall)
}

point <- stats_from(strata)
n_tab <- strata %>% group_by(level, comparison, stratum) %>%
  summarise(patients = n_distinct(individual_id), ordered_pairs = sum(w_pairs),
            cultures = sum(w_culture), .groups = "drop")

by_pat <- split(strata, strata$individual_id); pats <- names(by_pat)
boot <- map_dfr(seq_len(B), function(b) {
  stats_from(bind_rows(by_pat[sample(pats, length(pats), replace = TRUE)])) %>% mutate(rep = b)
})
ci <- boot %>% group_by(level, comparison, stratum, weighting, stat) %>%
  summarise(lo = quantile(value, 0.025, na.rm = TRUE),
            hi = quantile(value, 0.975, na.rm = TRUE), .groups = "drop")
res <- point %>% left_join(ci, by = c("level", "comparison", "stratum", "weighting", "stat")) %>%
  left_join(n_tab, by = c("level", "comparison", "stratum"))
write_csv(res, file.path(out_res, "sir_disagreement_all.csv"))

# Check: pair-weighted pooled within-culture rows reproduce AMR_sir_ordinal.R 2.2
old <- read_csv("results/sir/reproducibility_by_call.csv", show_col_types = FALSE) %>%
  filter(grepl("^B", null)) %>% transmute(index, old = obs_disagree)
new <- res %>% filter(level == "pooled", comparison == "Same culture", weighting == "pairs",
                      stat %in% c("P(S|S)", "P(I|I)", "P(R|R)")) %>%
  transmute(index = c(`P(S|S)` = "Susceptible", `P(I|I)` = "Intermediate",
                      `P(R|R)` = "Resistant")[stat], new = 1 - value)
chk <- inner_join(old, new, by = "index"); print(chk)
stopifnot(nrow(chk) == 3, max(abs(chk$old - chk$new)) < 1e-10)
cat("  pair-weighted rows reproduce AMR_sir_ordinal.R 2.2 (29 / 74 / 34%)\n")

# Is the morphotype split a real difference? bootstrap of the contrast
morph_diff <- boot %>%
  filter(level == "morph", stat %in% c("discord", "flip_S_vs_NS")) %>%
  select(rep, weighting, stat, stratum, value) %>%
  pivot_wider(names_from = stratum, values_from = value) %>%
  mutate(diff = `different morphotype` - `same morphotype`) %>%
  group_by(weighting, stat) %>%
  summarise(lo = quantile(diff, 0.025), hi = quantile(diff, 0.975),
            p_le0 = mean(diff <= 0), .groups = "drop") %>%
  left_join(point %>% filter(level == "morph", stat %in% c("discord", "flip_S_vs_NS")) %>%
              select(weighting, stat, stratum, value) %>%
              pivot_wider(names_from = stratum, values_from = value) %>%
              mutate(diff = `different morphotype` - `same morphotype`),
            by = c("weighting", "stat"))
write_csv(morph_diff, file.path(out_res, "morphotype_contrast.csv"))

# Same culture vs next culture, per class
ratio_tab <- function(d) d %>%
  pivot_wider(names_from = comparison, values_from = value) %>%
  mutate(ratio_all = `Same culture` / `Next culture`,
         ratio_matched = `Same culture` / `Next culture, multi-colony index`)
next_diff <- boot %>%
  filter(level %in% c("pooled", "class"), stat %in% c("discord", "flip_S_vs_NS")) %>%
  select(rep, weighting, stat, stratum, comparison, value) %>% ratio_tab() %>%
  group_by(weighting, stat, stratum) %>%
  summarise(across(c(ratio_all, ratio_matched),
                   list(lo = ~ quantile(.x, 0.025, na.rm = TRUE),
                        hi = ~ quantile(.x, 0.975, na.rm = TRUE))), .groups = "drop") %>%
  left_join(point %>% filter(level %in% c("pooled", "class"), stat %in% c("discord", "flip_S_vs_NS")) %>%
              select(weighting, stat, stratum, comparison, value) %>% ratio_tab(),
            by = c("weighting", "stat", "stratum"))
write_csv(next_diff, file.path(out_res, "same_vs_next_culture.csv"))


# #############################################################################
# FIGURE: drawn by analyses/figures_new.R from the CSVs written here, at print
# size (analyses/fig_style.R). This script no longer draws.
# #############################################################################

# #############################################################################
# REPORT
# #############################################################################

pc <- function(v, lo, hi) sprintf("%.0f%% (%.0f–%.0f)", 100 * v, 100 * lo, 100 * hi)
md_tab <- function(df) {
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
matrix_md <- function(lv, cmp, st, wt) {
  d <- res %>% filter(level == lv, comparison == cmp, stratum == st, weighting == wt,
                      grepl("^P\\(", stat)) %>%
    mutate(partner = substr(stat, 3, 3), index = substr(stat, 5, 5),
           cell = pc(value, lo, hi)) %>%
    select(index, partner, cell) %>%
    pivot_wider(names_from = partner, values_from = cell) %>%
    mutate(index = factor(index, c("S", "I", "R"))) %>% arrange(index) %>%
    select(`first \\ other` = index, S, I, R)
  md_tab(d)
}
summary_tab <- function(lv, wt) {
  res %>% filter(level == lv, weighting == wt,
                 stat %in% c("discord", "two_step", "flip_S_vs_NS", "flip_R_vs_SI")) %>%
    mutate(cell = pc(value, lo, hi)) %>%
    select(comparison, stratum, patients, cultures, stat, cell) %>%
    mutate(cultures = round(cultures)) %>%
    pivot_wider(names_from = stat, values_from = cell) %>%
    rename(`any S/I/R difference` = discord, `S<->R` = two_step,
           `S vs R+I flip` = flip_S_vs_NS, `R vs S+I flip` = flip_R_vs_SI)
}

writeLines(c(
  "# Within-culture S/I/R disagreement",
  "", sprintf("Nine drugs (%s). Patient-cluster bootstrap B = %d; 95%% percentile intervals.",
              paste(drugs_use, collapse = ", "), B),
  "'culture' weighting: each culture x drug counts once (primary). 'pairs': each ordered colony pair counts once (as in AMR_sir_ordinal.R).",
  "", "## Matrix, all drugs, same culture, culture-weighted", "",
  matrix_md("pooled", "Same culture", "All drugs", "culture"),
  "", "## Matrix, all drugs, same culture, pair-weighted (reproduces 29 / 74 / 34%)", "",
  matrix_md("pooled", "Same culture", "All drugs", "pairs"),
  "", "## Matrix, same culture, same morphotype (culture-weighted)", "",
  matrix_md("morph", "Same culture", "same morphotype", "culture"),
  "", "## Matrix, same culture, different morphotype (culture-weighted)", "",
  matrix_md("morph", "Same culture", "different morphotype", "culture"),
  "", "## Matrix, all drugs, next culture, multi-colony index cultures (culture-weighted)", "",
  matrix_md("pooled", "Next culture, multi-colony index", "All drugs", "culture"),
  "", "## Matrix, all drugs, next culture, all index cultures (culture-weighted)", "",
  matrix_md("pooled", "Next culture", "All drugs", "culture"),
  "", "## Overall discordance, all drugs and by class (culture-weighted)", "",
  md_tab(bind_rows(summary_tab("pooled", "culture"), summary_tab("class", "culture"))),
  "", "## By drug (culture-weighted)", "",
  md_tab(summary_tab("drug", "culture")),
  "", "## By morphotype pairing (culture-weighted)", "",
  md_tab(bind_rows(summary_tab("morph", "culture"), summary_tab("morph_class", "culture"))),
  "", "## Different minus same morphotype", "",
  md_tab(morph_diff %>% transmute(weighting, stat,
                                  same = sprintf("%.2f", `same morphotype`),
                                  different = sprintf("%.2f", `different morphotype`),
                                  diff = sprintf("%.2f (%.2f–%.2f)", diff, lo, hi),
                                  `P(diff <= 0)` = sprintf("%.3f", p_le0))),
  "", "## Same culture / next culture discordance ratio", "",
  md_tab(next_diff %>% filter(weighting == "culture") %>%
           transmute(stat, stratum, same = sprintf("%.2f", `Same culture`),
                     `next, all` = sprintf("%.2f", `Next culture`),
                     `ratio, all` = sprintf("%.2f (%.2f–%.2f)", ratio_all, ratio_all_lo, ratio_all_hi),
                     `next, multi-colony index` = sprintf("%.2f", `Next culture, multi-colony index`),
                     `ratio, matched` = sprintf("%.2f (%.2f–%.2f)", ratio_matched,
                                                ratio_matched_lo, ratio_matched_hi)))),
  file.path(out_res, "sir_disagreement.md"))

cat("\nwrote", out_res, "and", out_fig, "\n")
