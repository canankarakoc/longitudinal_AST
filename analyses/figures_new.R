# =============================================================================
# analyses/figures_new.R  -  the figures for analyses 1-6, drawn at print size
#
# Reads saved results only (fits nothing), so figures can be redrawn in
# seconds. The analysis scripts no longer draw; they write results/.
#
#   New 1  figures/reliability_split/reliability_split        analysis 1
#   New 2  figures/sir_matrix/sir_disagreement_matrix         analysis 2
#   New 3  figures/pseudo_trajectories/pseudo_trajectories    analysis 3
#   New 4  figures/deep_downsampling/deep_downsampling        analysis 3
#   New 5  figures/design_estimand/design_by_estimand         analysis 4
#   New 6  figures/power_recalibrated/power_frontiers         analysis 5
#   New 7  figures/calendar/calendar_effects                  analysis 6
#
# Style: analyses/fig_style.R (174 mm wide, 8 pt ticks, 9-10 pt text).
#
# Run from the repo root:  Rscript analyses/figures_new.R [1 2 ... 7]
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr); library(ggplot2); library(patchwork)
})
source("analyses/fig_style.R")
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
which_figs <- if (length(args)) as.integer(args) else 1:7

CLASSES   <- c("Aminoglycoside", "Carbapenem", "Fluoroquinolone", "Cephalosporin")
CLASS_COL <- c(Aminoglycoside = "#0072B2", Carbapenem = "#D55E00",
               Fluoroquinolone = "#009E73", Cephalosporin = "#CC79A7",
               Penicillin = "#E69F00", Monobactam = "#56B4E9")
rd <- function(f) read_csv(f, show_col_types = FALSE)
tag_theme <- theme(plot.tag = element_text(size = 11, face = "bold"))


# #############################################################################
# NEW 1  Reliability: cross-sectional vs within-patient
# #############################################################################
if (1 %in% which_figs) {
  cv <- rd("results/reliability_split/reliability_split_curves.csv") %>%
    mutate(class = factor(class, CLASSES),
           frame = factor(ifelse(frame == "colony", "Colony frame", "Drug-test frame"),
                          c("Drug-test frame", "Colony frame")))
  ln <- cv %>%
    pivot_longer(c(`Cross-sectional (pooled, Fig 3)`, `Within-patient deviation`, `Visit-to-visit change`),
                 names_to = "quantity", values_to = "R") %>%
    mutate(quantity = factor(recode(quantity,
             `Cross-sectional (pooled, Fig 3)` = "Cross-sectional (between patients)",
             `Within-patient deviation` = "Within patient",
             `Visit-to-visit change` = "Visit-to-visit change"),
             c("Cross-sectional (between patients)", "Within patient", "Visit-to-visit change")))
  p <- ggplot(ln, aes(colonies, R)) +
    geom_ribbon(data = cv, aes(x = colonies, ymin = `Within-patient, lower CI`,
                               ymax = `Within-patient, upper CI`),
                inherit.aes = FALSE, fill = "grey85") +
    geom_hline(yintercept = 0.8, linetype = 2, linewidth = 0.35, colour = "grey40") +
    geom_line(aes(linetype = quantity, colour = class), linewidth = 0.8) +
    facet_grid(frame ~ class) +
    scale_colour_manual(values = CLASS_COL, guide = "none") +
    scale_linetype_manual(values = c(1, 2, 3)) +
    scale_x_continuous(breaks = c(1, 4, 8, 12)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.4, 0.8)) +
    labs(x = "Colonies per culture", y = "Reliability", linetype = NULL,
         caption = "Grey band: patient-bootstrap 95% interval for the within-patient curve. Dashed line: 0.8.") +
    theme_fig() + theme(legend.position = "top", strip.text.y = element_text(angle = -90, hjust = 0.5))
  save_fig(p, "figures/reliability_split/reliability_split", FIG_W2, 120)
  cat("  new 1 done\n")
}


# #############################################################################
# NEW 2  Within-culture S/I/R disagreement
# #############################################################################
if (2 %in% which_figs) {
  res <- rd("results/sir_matrix/sir_disagreement_all.csv")
  NEXT <- "Next culture, multi-colony index"
  mat <- bind_rows(
    res %>% filter(level == "morph", weighting == "culture"),
    res %>% filter(level == "pooled", comparison == NEXT, weighting == "culture") %>%
      mutate(stratum = "next culture")) %>%
    filter(grepl("^P\\(", stat)) %>%
    mutate(partner = factor(substr(stat, 3, 3), c("S", "I", "R")),
           index   = factor(substr(stat, 5, 5), c("R", "I", "S")),
           panel = factor(c(`same morphotype` = "Same culture,\nsame morphotype",
                            `different morphotype` = "Same culture,\ndifferent morphotype",
                            `next culture` = "Next culture\n(same index cultures)")[stratum],
                          c("Same culture,\nsame morphotype", "Same culture,\ndifferent morphotype",
                            "Next culture\n(same index cultures)")))
  pa <- ggplot(mat, aes(partner, index, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.0f", 100 * value), colour = value > 0.5), size = pt(9)) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = "#F2F2F2", high = "#08306B", limits = c(0, 1), guide = "none") +
    facet_wrap(~panel, nrow = 1) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Call of the other colony", y = "Call of the first colony",
         title = "Chance of the other colony's call, given the first (%)") +
    theme_fig() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                        aspect.ratio = 1, panel.spacing = unit(14, "pt"))

  cls <- bind_rows(
    res %>% filter(level == "morph_class", weighting == "culture", stat == "discord") %>%
      separate(stratum, c("class", "what"), sep = " / "),
    res %>% filter(level == "class", comparison == NEXT, weighting == "culture", stat == "discord") %>%
      mutate(class = stratum, what = "next culture")) %>%
    mutate(class = factor(class, rev(c(CLASSES, "Penicillin", "Monobactam"))),
           what = factor(what, c("same morphotype", "different morphotype", "next culture"),
                         c("Same culture, same morphotype", "Same culture, different morphotype",
                           "Next culture")))
  pb <- ggplot(cls, aes(value, class, colour = what)) +
    geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.65),
                    size = 0.35, linewidth = 0.6) +
    scale_colour_manual(values = c("#4D4D4D", "#D55E00", "#0072B2")) +
    scale_x_continuous(limits = c(0, 0.6), labels = scales::percent) +
    labs(x = "Two colonies give different S/I/R calls", y = NULL, colour = NULL,
         title = "Discordance by class (each culture counted once; 95% CI)") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_fig() + theme(legend.position = "top")
  p <- (pa / pb) + plot_layout(heights = c(1, 1.5)) + plot_annotation(tag_levels = "a") & tag_theme
  save_fig(p, "figures/sir_matrix/sir_disagreement_matrix", FIG_W2, 200)
  cat("  new 2 done\n")
}


# #############################################################################
# NEW 3  Pseudo-trajectories
# #############################################################################
if (3 %in% which_figs) {
  per_pat <- rd("results/pseudo_trajectories/pseudo_trajectories_per_patient.csv") %>%
    mutate(id = as.character(id))
  cohort  <- rd("results/pseudo_trajectories/pseudo_trajectories_cohort.csv")
  iso <- readRDS("data/pa_isolates.rds")
  AG <- c("amikacin", "tobramycin", "gentamicin"); SIR <- c("Susceptible", "Intermediate", "Resistant")
  cult <- tibble(id = as.character(iso$individual_id), date = as.Date(iso$date_collected_from_subject),
                 tests = rowSums(sapply(AG, function(x) iso[[x]] %in% SIR)),
                 ns = rowSums(sapply(AG, function(x) iso[[x]] %in% c("Resistant", "Intermediate"))) > 0) %>%
    filter(tests > 0) %>% group_by(id, date) %>% summarise(m = n(), k = sum(ns), .groups = "drop") %>%
    mutate(p = k / m) %>% arrange(id, date)
  ex_ids <- per_pat %>%
    filter(variant == "All cultures", class == "Aminoglycoside", cultures >= 8) %>%
    arrange(desc(e_switch_pick / pmax(e_switch, 1e-9)), desc(het)) %>% slice_head(n = 3) %>% pull(id)
  NPICK <- 30                                         # matches the Fig 6 caption
  ex <- map_dfr(ex_ids, function(i) {
    d <- cult %>% filter(id == i)
    x <- matrix(rbinom(NPICK * nrow(d), 1, rep(d$p, each = NPICK)), nrow = NPICK)
    bind_rows(tibble(row = 0L, j = seq_len(nrow(d)), ns = as.integer(d$k > 0)),
              map_dfr(seq_len(NPICK), ~ tibble(row = .x, j = seq_len(nrow(d)), ns = x[.x, ]))) %>%
      mutate(id = paste("Patient", i))
  }) %>% mutate(y = ifelse(row == 0, NPICK + 2, NPICK + 1 - row),
                state = factor(ifelse(ns == 1, "Non-susceptible", "Susceptible"),
                               c("Susceptible", "Non-susceptible")))
  pa <- ggplot(ex, aes(j, y, fill = state)) +
    geom_tile(colour = "white", linewidth = 0.25, height = 0.95) +
    facet_wrap(~id, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c(Susceptible = "#D9E3EA", `Non-susceptible` = "#B2182B")) +
    scale_y_continuous(breaks = c(NPICK + 2, NPICK / 2 + 0.5),
                       labels = c("All colonies", sprintf("%d one-colony\nrecords", NPICK)),
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = function(l) unique(pmax(1, round(pretty(l, 4)))), expand = c(0, 0)) +
    labs(x = "Culture, in time order", y = NULL, fill = NULL,
         title = "Aminoglycoside: the same cultures, one colony picked at random") +
    theme_fig() + theme(legend.position = "top", axis.line = element_blank(), axis.ticks.y = element_blank())

  pb_dat <- cohort %>% filter(variant == "All cultures") %>%
    transmute(class = factor(class, rev(CLASSES)),
              `Fixed by the data` = e_switch_per_pat * (1 - share_switch_pick),
              `Set by which colony was picked` = e_switch_per_pat * share_switch_pick) %>%
    pivot_longer(-class, names_to = "part") %>%
    mutate(part = factor(part, c("Set by which colony was picked", "Fixed by the data")))
  pb <- ggplot(pb_dat, aes(value, class, fill = part)) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = c(`Fixed by the data` = "grey70", `Set by which colony was picked` = "#B2182B")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Expected apparent switches between S and NS, per patient", y = NULL, fill = NULL,
         title = "Where apparent switches come from") +
    theme_fig() + theme(legend.position = "top")

  TYPES <- c("stable S", "stable NS", "acquired", "lost", "fluctuating")
  type_cols <- c(`stable S` = "#D9E3EA", `stable NS` = "#B2182B", acquired = "#F4A582",
                 lost = "#92C5DE", fluctuating = "#5E3C99")
  pc_dat <- per_pat %>% filter(variant == "All cultures") %>%
    select(class, id, p_modal, starts_with("P_", ignore.case = FALSE)) %>%
    pivot_longer(starts_with("P_", ignore.case = FALSE), names_to = "type", values_to = "prob") %>%
    mutate(type = factor(gsub("_", " ", sub("^P_", "", type)), TYPES), class = factor(class, CLASSES)) %>%
    group_by(class) %>% mutate(idc = reorder(paste(class, id), -p_modal)) %>% ungroup()
  pc <- ggplot(pc_dat, aes(idc, prob, fill = type)) +
    geom_col(width = 1) +
    facet_wrap(~class, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = type_cols) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = "Patients, sorted by the probability of their most likely history",
         y = "One-colony records", fill = NULL, title = "Apparent history depends on the pick") +
    theme_fig() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                        legend.position = "top")
  p <- (pa / pb / pc) + plot_layout(heights = c(1.5, 0.75, 1)) +
    plot_annotation(tag_levels = "a") & tag_theme
  save_fig(p, "figures/pseudo_trajectories/pseudo_trajectories", FIG_W2, 225)
  cat("  new 3 done\n")
}


# #############################################################################
# NEW 4  Deep-sample downsampling
# #############################################################################
if (4 %in% which_figs) {
  down <- rd("results/deep_downsampling/deep_downsampling.csv")
  noise <- rd("results/deep_downsampling/assay_vs_colony_discordance.csv")
  short_sample <- function(s) sub(" \\(cohort id (\\d+)\\)", " (pt \\1)", sub(" \\(not in cohort\\)", "", s))
  short_drug <- c(amikacin = "amikacin", ciprofloxacin = "ciprofloxacin", ceftazidime = "ceftazidime",
                  `piperacillin-tazobactam` = "pip-tazo", meropenem = "meropenem", tobramycin = "tobramycin")
  het <- down %>% filter(heterogeneous_SIR) %>%
    mutate(cell = sprintf("%s %s (S/I/R %.0f/%.0f/%.0f%%)", short_sample(sample), short_drug[drug],
                          100 * pS, 100 * pI, 100 * pR))
  cell_cols <- setNames(c("#0072B2", "#56B4E9", "#D55E00", "#009E73", "#CC79A7"), sort(unique(het$cell)))
  xs <- scale_x_continuous(trans = "log2", breaks = c(1, 2, 4, 8, 15, 30))
  pa <- het %>% pivot_longer(c(switch_plurality, switch_worst), names_to = "rule") %>%
    mutate(rule = recode(rule, switch_plurality = "Report the most common call",
                         switch_worst = "Report the most resistant colony")) %>%
    ggplot(aes(n, value, colour = cell, linetype = mode)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1) +
    facet_wrap(~rule) + xs +
    scale_colour_manual(values = cell_cols) +
    scale_linetype_manual(values = c(mean = 1, single = 2),
                          labels = c(mean = "colony sampling only", single = "+ single-disc assay noise")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(x = "Colonies tested", y = "Two tests of the same\nsputum disagree", colour = NULL, linetype = NULL,
         title = "Apparent change with no true change") +
    guides(colour = guide_legend(ncol = 2, order = 1), linetype = guide_legend(order = 2)) +
    theme_fig() + theme(legend.position = "top", legend.box = "vertical", legend.box.just = "left")
  pb <- het %>% filter(mode == "mean", theta_NS > 0, theta_NS < 1) %>%
    pivot_longer(c(detect, freq_abs_err), names_to = "what") %>%
    mutate(what = recode(what, detect = "Detect any NS colony",
                         freq_abs_err = "Mean absolute error of the NS fraction")) %>%
    ggplot(aes(n, value, colour = cell)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1) +
    facet_wrap(~what, scales = "free_y") + xs +
    scale_colour_manual(values = cell_cols) +
    labs(x = "Colonies tested", y = NULL, colour = NULL,
         title = "Minority detection and frequency estimation") +
    guides(colour = "none") + theme_fig()
  pc_dat <- noise %>% filter(pmax(assay_rep_disc, colony_disc_mean, colony_disc_single) > 0.02) %>%
    mutate(cell = sprintf("%s %s", short_sample(sample), short_drug[drug])) %>%
    pivot_longer(c(assay_rep_disc, colony_disc_mean, colony_disc_single)) %>%
    mutate(name = factor(recode(name, assay_rep_disc = "Same isolate, two discs",
                                colony_disc_mean = "Two isolates, mean zone",
                                colony_disc_single = "Two isolates, one disc each"),
                         c("Same isolate, two discs", "Two isolates, mean zone", "Two isolates, one disc each")))
  pc <- ggplot(pc_dat, aes(value, reorder(cell, value), colour = name)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = c("#E69F00", "#0072B2", "#56B4E9")) +
    scale_x_continuous(labels = scales::percent, limits = c(0, 0.7)) +
    labs(x = "Different S/I/R call", y = NULL, colour = NULL,
         title = "Assay noise vs colony-to-colony differences (cells with any disagreement)") +
    guides(colour = guide_legend(nrow = 1)) + theme_fig() + theme(legend.position = "top")
  pc <- if (exists("free", asNamespace("patchwork"))) patchwork::free(pc) else pc
  p <- (pa / pb / pc) + plot_layout(heights = c(1.2, 0.8, 1.05)) +
    plot_annotation(tag_levels = "a") & tag_theme
  save_fig(p, "figures/deep_downsampling/deep_downsampling", FIG_W2, 250)
  cat("  new 4 done\n")
}


# #############################################################################
# NEW 5  Sampling design by estimand
# #############################################################################
if (5 %in% which_figs) {
  sim <- rd("results/design_estimand/design_by_estimand_grid.csv")
  DES <- c(`random-1` = "Random, 1 colony", `random-4` = "Random, 4 colonies", `random-8` = "Random, 8 colonies",
           `stratified-1` = "1 per morphotype", `stratified-2` = "2 per morphotype",
           `pooled-10, NS if >= 1%` = "Pool of 10, reads NS at ≥1%",
           `pooled-10, NS if >= 50%` = "Pool of 10, reads NS at ≥50%",
           `whole sputum, NS if >= 1%` = "Whole sputum, ≥1%",
           `whole sputum, NS if >= 50%` = "Whole sputum, ≥50%")
  DCOL <- setNames(c("#BDBDBD", "#737373", "#000000", "#6BAED6", "#08519C",
                     "#FDAE6B", "#D94801", "#BCBDDC", "#6A51A3"), DES)
  EST <- c(p_report_NS = "Detect\nP(report NS)", freq_mae = "Frequency\nmean abs. error",
           dominant_correct = "Dominant\nP(correct)", false_change_detect = "Change\nP(false change)")
  fd <- sim %>% filter(eps == 0) %>%
    select(structure, theta, design, all_of(names(EST))) %>%
    pivot_longer(all_of(names(EST)), names_to = "estimand") %>% filter(!is.na(value)) %>%
    mutate(estimand = factor(EST[estimand], EST),
           structure = factor(ifelse(structure == "aligned", "Resistance in one morphotype",
                                     "Resistance spread across morphotypes"),
                              c("Resistance in one morphotype", "Resistance spread across morphotypes")),
           design = factor(DES[design], DES))
  p <- ggplot(fd, aes(theta, value, colour = design)) +
    geom_line(linewidth = 0.75) +
    facet_grid(estimand ~ structure, scales = "free_y", switch = "y") +
    scale_colour_manual(values = DCOL) +
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    labs(x = "True non-susceptible fraction in the sputum", y = NULL, colour = NULL,
         caption = paste("Pool and whole-sputum readouts follow ASSUMED rules (NS when ≥1% or ≥50% is NS).",
                         "Change = two visits with the same fraction report different detection calls.",
                         "No assay error here.", sep = "\n")) +
    guides(colour = guide_legend(ncol = 3, byrow = TRUE)) +
    theme_fig() +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90, hjust = 0.5),
          legend.position = "bottom", legend.key.width = unit(16, "pt"))
  save_fig(p, "figures/design_estimand/design_by_estimand", FIG_W2, 215)
  cat("  new 5 done\n")
}


# #############################################################################
# NEW 6  Power by lever, recalibrated
# #############################################################################
if (6 %in% which_figs) {
  g <- rd("results/power_recalibrated/power_frontiers.csv") %>% filter(params == "Aminoglycoside (fitted tau)")
  eff <- rd("results/power_recalibrated/power_effect_sizes.csv")
  LEV <- c(`interval (L fixed)` = "Interval (days)\n12 visits",
           `visits (follow-up fixed)` = "Visits per patient\nsame follow-up",
           `follow-up (interval fixed)` = "Visits per patient\nsame 84-d interval",
           depth = "Colonies\nper visit")
  BUD <- c(`N fixed at 90` = "Same 90 patients", `fixed cost: base` = "Fixed budget",
           `fixed cost: visits dear` = "Fixed budget, cultures dear")
  fd <- g %>% pivot_longer(c(power_coupling, power_step), names_to = "estimand", values_to = "power") %>%
    mutate(estimand = factor(ifelse(estimand == "power_coupling",
                                    sprintf("Coupling\n(β = %.1f at 84 d)", eff$beta_at_84d),
                                    sprintf("Step change\n(%g logit)", eff$delta_logit))),
           lever = factor(LEV[lever], LEV), budget = factor(BUD[budget], BUD))
  base_pts <- fd %>% filter((lever == LEV[1] & x == 84) | (lever %in% LEV[2:3] & x == 12) |
                              (lever == LEV[4] & x == 1.5), budget == BUD[1])
  p <- ggplot(fd, aes(x, power, colour = budget)) +
    geom_hline(yintercept = 0.8, linetype = 2, linewidth = 0.35, colour = "grey50") +
    geom_line(linewidth = 0.75) + geom_point(size = 1.3) +
    geom_point(data = base_pts, shape = 21, size = 3.4, colour = "black", fill = NA, stroke = 0.6) +
    facet_grid(estimand ~ lever, scales = "free_x") +
    scale_x_continuous(trans = "log2", breaks = function(lim) {
      if (lim[2] > 100) c(14, 28, 84, 168) else if (lim[2] > 10) c(4, 12, 24, 48) else c(1, 2, 4, 8)
    }) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.4, 0.8)) +
    scale_colour_manual(values = c("grey55", "#0072B2", "#D55E00")) +
    labs(x = NULL, y = "Power", colour = NULL,
         caption = paste("Circle: the cohort's own design (90 patients, 12 visits, 84 days, 1.5 colonies).",
                         "Aminoglycoside calibration. Budget weights: colony 1, culture 2 (dear: 10), patient 20.",
                         sep = "\n")) +
    theme_fig() + theme(legend.position = "top", strip.text.y = element_text(angle = -90, hjust = 0.5))
  save_fig(p, "figures/power_recalibrated/power_frontiers", FIG_W2, 130)
  cat("  new 6 done\n")
}


# #############################################################################
# NEW 7  Calendar effects
# #############################################################################
if (7 %in% which_figs) {
  half_date <- function(h) as.Date(paste0(substr(h, 1, 4), ifelse(grepl("H1", h), "-04-01", "-10-01")))
  CLSI <- as.Date("2019-01-01"); COVID <- as.Date("2020-03-01")
  best <- rd("results/calendar/changepoint_best.csv")
  c_star <- best$change[best$drug == "levofloxacin"]
  xs <- scale_x_date(breaks = as.Date(paste0(2017:2022, "-01-01")), date_labels = "%Y",
                     limits = as.Date(c("2017-01-01", "2023-01-01")), expand = c(0, 0))
  marks <- list(geom_vline(xintercept = CLSI, linetype = 2, linewidth = 0.4, colour = "grey30"),
                geom_vline(xintercept = COVID, linetype = 3, linewidth = 0.5, colour = "grey30"))

  des <- rd("results/calendar/design_by_half_year.csv") %>% mutate(date = half_date(half)) %>%
    select(date, cultures, colonies_tested_mean, pct_single_tested, median_gap_days) %>%
    pivot_longer(-date) %>%
    mutate(name = factor(recode(name, cultures = "Cultures per half-year",
                                colonies_tested_mean = "Colonies tested per culture",
                                pct_single_tested = "Single-colony cultures (%)",
                                median_gap_days = "Median days between cultures"),
                         c("Cultures per half-year", "Colonies tested per culture",
                           "Single-colony cultures (%)", "Median days between cultures")))
  pa <- ggplot(des, aes(date, value)) + marks +
    geom_line(linewidth = 0.6) + geom_point(size = 1.3) +
    facet_wrap(~name, ncol = 2, scales = "free_y") + xs +
    labs(x = NULL, y = NULL, title = "Sampling design") + theme_fig()

  comp <- rd("results/calendar/sir_composition_by_half_year.csv") %>% filter(drug == "levofloxacin") %>%
    mutate(date = half_date(half)) %>% pivot_longer(c(S, I, R), names_to = "call") %>%
    mutate(call = factor(call, c("R", "I", "S")))
  pb <- ggplot(comp, aes(date, value, fill = call)) +
    geom_col(width = 160) + marks +
    geom_vline(xintercept = c_star, colour = "#D55E00", linewidth = 0.7) +
    annotate("text", x = CLSI + 20, y = 1.04, label = "CLSI revision", hjust = 0, vjust = 0, size = pt(8)) +
    scale_fill_manual(values = c(S = "#D9E3EA", I = "#F4A582", R = "#B2182B")) +
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1), limits = c(0, 1.1)) + xs +
    labs(x = NULL, y = "Colonies", fill = NULL,
         title = sprintf("Levofloxacin calls by half-year (orange: best step, %s)", format(c_star, "%b %Y"))) +
    theme_fig() + theme(legend.position = "right")

  DRUG_LAB <- c(levofloxacin = "Levofloxacin", tobramycin = "Tobramycin", amikacin = "Amikacin",
                gentamicin = "Gentamicin", meropenem = "Meropenem", cefepime = "Cefepime",
                ceftazidime = "Ceftazidime", piperacillin_tazobactam = "Pip-tazo", aztreonam = "Aztreonam")
  scan <- rd("results/calendar/changepoint_scan.csv") %>%
    mutate(drug = factor(DRUG_LAB[drug], DRUG_LAB))
  pc <- ggplot(scan, aes(change, lr)) + marks[1] +
    geom_line(linewidth = 0.6) +
    facet_wrap(~drug, ncol = 3) +
    scale_x_date(breaks = as.Date(paste0(c(2018, 2020, 2022), "-01-01")), date_labels = "%Y") +
    scale_y_continuous(breaks = c(0, 20, 40)) +
    labs(x = "Candidate date of a step", y = "Evidence for a step (LR)",
         title = "Change-point scan, every drug",
         caption = paste("Dashed: CLSI M100 29th ed., Jan 2019. Dotted: March 2020.",
                         "Tested against a null that keeps within-patient dynamics, only fluoroquinolones step at Jan 2019",
                         "(p = 0.005); carbapenems rise in late 2017 (p = 0.015).", sep = "\n")) +
    theme_fig()
  p <- (pa / pb / pc) + plot_layout(heights = c(1.15, 0.6, 1.3)) +
    plot_annotation(tag_levels = "a") & tag_theme
  save_fig(p, "figures/calendar/calendar_effects", FIG_W2, 240)
  cat("  new 7 done\n")
}
