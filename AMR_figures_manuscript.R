# =============================================================================
# AMR_figures_manuscript.R  -  Figures 1-5 for the AST measurement manuscript
#
# Author: Canan Karakoc
#
# "Sampling depth, not sampling frequency, sets what susceptibility testing can
#  resolve in chronic Pseudomonas aeruginosa infection"
#
#   Figure 1  Sampling depth in routine practice (colonies and drug-tests per
#             culture, by class)
#   Figure 2  Variance components across four drug classes (Table 1)
#   Figure 3  Reliability against colonies per culture
#   Figure 4  Within-patient autocorrelation against visit interval
#   Figure 5  Design frontier: patients against visits per patient at 80% power
#
# Reads data/ and results/ only; fits no brms model. Figure 4 recomputes the
# binned autocorrelation from pa_visit_full (and checks it against
# results/design/interval_autocorrelation.csv) so it can add a patient-level
# bootstrap. Everything else is read from the pipeline's result CSVs.
#
# Deliberately does NOT source AMR_common.R: that file turns on showtext with
# Lato and a transparent talk theme, which is exactly what these figures need
# to get away from. White background, Helvetica, journal column widths.
#
# Output: figures/manuscript/  PDF (vector, Helvetica embedded), PNG and TIFF
#         (600 dpi). results/manuscript/  every number drawn, for checking text.
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
  library(ggplot2); library(patchwork); library(scales)
})

out_fig <- "figures/manuscript"
out_res <- "results/manuscript"
dir.create(file.path(out_fig, "variants"), recursive = TRUE, showWarnings = FALSE)
dir.create(out_res, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

# Which reliability calculation Figure 3 draws. See FIGURE 3 for what each does.
#   "class_deff" - each class's own design effect, E[1/m] moment inversion
#   "as_run"     - AMR_pipeline.R 14.1 exactly (AG design effect for all classes)
#   "colony"     - colony frame: k of n colonies non-susceptible, no design effect
# The two not chosen are written to figures/manuscript/variants/.
REL_METHOD <- "class_deff"

N_BOOT <- 1000   # patient-level bootstrap replicates for Figure 4


# #############################################################################
# STYLE
# #############################################################################

FONT <- "Helvetica"

theme_pub <- function(base_size = 9.5) {   # 8 -> 9.5: text was too small to read at print size
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      axis.line         = element_line(colour = "black", linewidth = 0.3),
      axis.ticks        = element_line(colour = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.text         = element_text(colour = "black", size = base_size - 1),
      axis.title        = element_text(colour = "black", size = base_size),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_size - 1),
      legend.key.size   = unit(8, "pt"),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, hjust = 0,
                                       margin = margin(b = 3)),
      plot.tag          = element_text(size = base_size + 2, face = "bold"),
      plot.margin       = margin(4, 6, 4, 4)
    )
}

# Okabe-Ito: colour-blind safe, and the same class always gets the same colour.
# Order follows Table 1.
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
NS  <- c("Resistant", "Intermediate")          # conservative coding, R + I

# PDF through the Quartz device: it embeds Helvetica and handles Unicode
# (rho, tau, >=). The base pdf() device does neither, and cairo_pdf needs
# XQuartz. Off macOS, fall back to pdf() with a warning.
pdf_dev <- function(filename, width, height, ...) {
  if (capabilities("aqua")) {
    grDevices::quartz(file = filename, type = "pdf", width = width,
                      height = height, family = FONT, bg = "white")
  } else {
    warning("Quartz unavailable: PDF fonts will not be embedded.")
    grDevices::pdf(filename, width = width, height = height, family = FONT)
  }
}

save_fig <- function(p, name, width_mm, height_mm, dir = out_fig) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = pdf_dev)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(file.path(dir, paste0(name, ".tiff")), p, width = w, height = h,
         dpi = 600, device = ragg::agg_tiff, bg = "white", compression = "lzw")
  cat("  wrote", file.path(dir, name), sprintf("(%.0f x %.0f mm)\n", width_mm, height_mm))
}

# Journal widths (mm). Most CF / microbiology journals: 85 single, 174-180 double.
W1 <- 85; W15 <- 120; W2 <- 174


# #############################################################################
# DATA
# #############################################################################

iso <- readRDS("data/pa_isolates.rds")
vis <- readRDS("data/pa_visit_full.rds")
cat(sprintf("\nLoaded %d isolates, %d visits, %d patients\n",
            nrow(iso), nrow(vis), n_distinct(vis$individual_id)))

tested_matrix <- function(drugs, levels) {
  sapply(drugs, function(x) iso[[x]] %in% levels)
}

# One row per culture (patient x collection date) per class, counting colonies
# with at least one result in the class and the drug-test results they carry.
depth <- imap_dfr(CLASS_DRUGS, function(drugs, cl) {
  n_tests <- rowSums(tested_matrix(drugs, SIR))
  n_ns    <- rowSums(tested_matrix(drugs, NS))
  tibble(individual_id = iso$individual_id,
         visit_date    = as.Date(iso$date_collected_from_subject),
         tests = n_tests, colony_ns = n_ns > 0) %>%
    filter(tests > 0) %>%
    group_by(individual_id, visit_date) %>%
    summarise(colonies = n(), colonies_ns = sum(colony_ns), tests = sum(tests),
              .groups = "drop") %>%
    mutate(class = cl)
}) %>% mutate(class = factor(class, levels = CLASSES))

depth_summary <- depth %>%
  group_by(class) %>%
  summarise(cultures        = n(),
            patients        = n_distinct(individual_id),
            colonies_mean   = mean(colonies),
            pct_one_colony  = 100 * mean(colonies == 1),
            n_two_plus      = sum(colonies >= 2),
            tests_mean      = mean(tests),
            pct_one_test    = 100 * mean(tests == 1),
            tests_per_colony = sum(tests) / sum(colonies),
            .groups = "drop")
cat("\nSampling depth by class (cultures with >= 1 colony tested in the class)\n")
print(as.data.frame(depth_summary), digits = 3, row.names = FALSE)
write_csv(depth_summary, file.path(out_res, "fig1_depth_summary.csv"))

all_cultures <- iso %>% count(individual_id, date_collected_from_subject, name = "n")
cat(sprintf(paste0("All PA cultures, counting every isolate whether or not it was ",
                   "tested: %d cultures, mean %.2f colonies, %.1f%% single\n"),
            nrow(all_cultures), mean(all_cultures$n), 100 * mean(all_cultures$n == 1)))


# #############################################################################
# FIGURE 1 - sampling depth in routine practice
#
# (a) colonies carried to AST per culture, by class. Near-identical across
#     classes: a colony is tested against the whole panel.
# (b) the same cultures counted in drug-test results. AG looks three times
#     deeper than it is because one colony carries three AG results; carbapenem
#     and FQ carry one (imipenem and ciprofloxacin are almost never run).
# #############################################################################

col_cap  <- 5
test_cap <- 9

f1a_dat <- depth %>%
  mutate(x = pmin(colonies, col_cap)) %>%
  count(class, x) %>%
  group_by(class) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()

f1b_dat <- depth %>%
  mutate(x = pmin(tests, test_cap)) %>%
  count(class, x) %>%
  group_by(class) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>%
  complete(class, x = 1:test_cap, fill = list(n = 0, pct = 0))

lab_cap <- function(cap) function(b) ifelse(b == cap, paste0("≥", cap), b)

f1a <- ggplot(f1a_dat, aes(x, pct, fill = class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_x_continuous(breaks = 1:col_cap, labels = lab_cap(col_cap)) +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0, 60, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = CLASS_COL) +
  labs(x = "Colonies tested per culture", y = "Cultures (%)") +
  theme_pub() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))

f1b_means <- depth_summary %>%
  transmute(class, lab = sprintf("mean %.1f", tests_mean))

f1b <- ggplot(f1b_dat, aes(x, pct, fill = class)) +
  geom_col(width = 0.75) +
  geom_text(data = f1b_means, aes(x = test_cap, y = 60, label = lab),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = 8 / .pt,
            family = FONT) +
  facet_wrap(~class, nrow = 1) +
  scale_x_continuous(breaks = c(1, 3, 6, test_cap),
                     labels = lab_cap(test_cap)(c(1, 3, 6, test_cap))) +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0, 60, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = CLASS_COL, guide = "none") +
  labs(x = "Drug-test results per culture", y = "Cultures (%)") +
  theme_pub()

fig1 <- f1a / f1b + plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a")
save_fig(fig1, "Fig1_sampling_depth", W2, 110)


# #############################################################################
# FIGURE 2 - variance components across four drug classes (Table 1)
#
# Source: results/noise/variance_decomposition_4class.csv (pipeline 3.2).
# Isolate-level Bernoulli, patient and visit-within-patient random effects;
# shares on the latent scale with the residual fixed at pi^2/3, so the
# between-colony share is partly the link-function convention (legend says so).
# #############################################################################

vd <- read_csv("results/noise/variance_decomposition_4class.csv", show_col_types = FALSE)

comp_lab <- c(between_patient              = "Between patients",
              between_isolate_within_visit = "Between colonies\n(within culture)",
              between_visit_within_patient = "Between visits\n(within patient)")

f2_dat <- vd %>%
  mutate(class = factor(class, levels = CLASSES),
         comp  = factor(comp_lab[component], levels = comp_lab),
         # manual vertical dodge: classes stacked within each component
         y = as.numeric(comp) + (as.numeric(class) - 2.5) * 0.17)

write_csv(f2_dat %>% select(class, component, median, lo, hi),
          file.path(out_res, "fig2_variance_components.csv"))

# (b) the aminoglycoside drug x isolate model (pipeline 10.3/10.3b), where the
# between-colony term is ESTIMATED and the fixed residual is drug-to-drug
# variation on one colony. Shares are computed per posterior draw here, so they
# carry intervals; the pipeline CSV reports shares of posterior-median
# variances. Draws are read once from the two cached fits and kept in a CSV.
f2b_cache <- file.path(out_res, "fig2b_estimated_shares.csv")
if (!file.exists(f2b_cache)) {
  fits <- c(`All visits (1,061)` = "models/sensitivity/m_var_isolate.rds",
            `Visits with ≥ 2 colonies (443)` = "models/sensitivity/m_var_isolate_multi.rds")
  f2b_raw <- imap_dfr(fits, function(f, lab) {
    dr <- brms::as_draws_df(readRDS(f))
    v <- tibble(Patient = dr$sd_individual_id__Intercept^2,
                Visit   = dr$sd_visit_key__Intercept^2,
                Colony  = dr$sd_isolate_id__Intercept^2,
                Residual = pi^2 / 3)
    tot <- rowSums(v)
    sh <- map_dfr(names(v), function(nm) {
      s <- v[[nm]] / tot
      tibble(component = nm, median = median(s),
             lo = quantile(s, 0.025), hi = quantile(s, 0.975),
             share_of_median_variances = median(v[[nm]]) /
               (median(v$Patient) + median(v$Visit) + median(v$Colony) + pi^2 / 3))
    })
    sdv <- dr$sd_visit_key__Intercept
    sh %>% mutate(subset = lab, sd_visit_median = median(sdv),
                  sd_visit_lo = quantile(sdv, 0.025), sd_visit_hi = quantile(sdv, 0.975))
  })
  write_csv(f2b_raw, f2b_cache)
}
f2b <- read_csv(f2b_cache, show_col_types = FALSE)
cat("\nFigure 2b: estimated four-component shares (AG drug x isolate)\n")
print(as.data.frame(f2b %>% select(subset, component, median, lo, hi,
                                   share_of_median_variances)), digits = 3, row.names = FALSE)

comp_lab_b <- c(Patient  = "Between patients",
                Colony   = "Between colonies\n(within culture)",
                Visit    = "Between visits\n(within patient)",
                Residual = "Between agents on one\ncolony (fixed, π²/3)")
f2b_dat <- f2b %>%
  mutate(comp = factor(comp_lab_b[component], levels = comp_lab_b),
         subset = factor(subset, levels = unique(subset)),
         y = as.numeric(comp) + (as.numeric(subset) - 1.5) * 0.3)

fig2b <- ggplot(f2b_dat, aes(median, y, shape = subset)) +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = CLASS_COL[["Aminoglycoside"]],
                 linewidth = 0.55) +
  geom_point(size = 1.6, colour = CLASS_COL[["Aminoglycoside"]], fill = "white", stroke = 0.6) +
  scale_shape_manual(values = c(16, 21)) +
  scale_y_reverse(breaks = seq_along(comp_lab_b), labels = comp_lab_b,
                  expand = expansion(add = 0.45)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Share of variance (latent scale)", y = NULL) +
  guides(shape = guide_legend(nrow = 2)) +
  theme_pub() +
  theme(legend.position = "top", axis.ticks.y = element_blank(),
        legend.margin = margin(0, 0, 0, 0))

fig2 <- ggplot(f2_dat, aes(median, y, colour = class)) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.55) +
  geom_point(size = 1.5) +
  scale_y_reverse(breaks = seq_along(comp_lab), labels = comp_lab,
                  expand = expansion(add = 0.45)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = expansion(mult = c(0.01, 0.02))) +
  scale_colour_manual(values = CLASS_COL) +
  labs(x = "Share of variance (latent scale)", y = NULL) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_pub() +
  theme(legend.position = "top", axis.ticks.y = element_blank(),
        legend.margin = margin(0, 0, 0, 0))

fig2 <- (fig2 | fig2b) + plot_annotation(tag_levels = "a")
save_fig(fig2, "Fig2_variance_components", W2, 82)


# #############################################################################
# FIGURE 3 - reliability against colonies per culture
#
# A visit's resistance frequency is k/m from m tests. Its variance is signal
# plus binomial noise:  V_obs = V_theta + E[theta(1 - theta)] * E[1/m]
# with E[theta(1 - theta)] = mu(1 - mu) - V_theta, which inverts to
#     V_theta = (V_obs - mu(1 - mu) E[1/m]) / (1 - E[1/m])
# on visits with m >= 2 (at m = 1 signal and noise are not separable).
# Reliability of a visit measured with m_eff independent tests:
#     R(m_eff) = V_theta / (V_theta + E[theta(1 - theta)] / m_eff)
# A colony tested against k correlated agents is worth
#     m_eff = k / (1 + (k - 1) rho)
# independent tests, where rho is the mean within-colony correlation between
# the agents' non-susceptible calls.
#
# Three versions:
#   as_run     - pipeline 14.1: 1/mean(m) in place of E[1/m], and the AG design
#                effect (k = 3, 1.40 tests per colony) applied to EVERY class.
#                Carbapenem and FQ colonies carry ONE result (imipenem and
#                ciprofloxacin are almost never run), so this over-credits them.
#   class_deff - E[1/m], and each class's own k and rho.
#   colony     - no design effect at all: the visit frequency is computed on
#                colonies (k of n colonies non-susceptible to any agent in the
#                class). This is the quantity the Discussion proposes reporting,
#                but it is a different estimand from the drug-test frequency.
# #############################################################################

# within-colony correlation of non-susceptible calls, per class
within_rho <- function(drugs) {
  dd <- iso %>% filter(if_all(all_of(drugs), ~ .x %in% SIR)) %>%
    transmute(across(all_of(drugs), ~ .x %in% NS))
  if (nrow(dd) < 30 || length(drugs) < 2) return(NA_real_)
  cmb <- combn(drugs, 2)
  mean(apply(cmb, 2, function(p) cor(dd[[p[1]]], dd[[p[2]]])))
}

deff_tab <- imap_dfr(CLASS_DRUGS, function(drugs, cl) {
  # agents routinely run per colony: the rounded mean number of results
  k   <- round(depth_summary$tests_per_colony[depth_summary$class == cl])
  # rho over the agents actually run; for k = 1 there is nothing to correct
  run <- names(sort(colSums(tested_matrix(drugs, SIR)), decreasing = TRUE))[seq_len(k)]
  rho <- if (k >= 2) within_rho(run) else NA_real_
  tibble(class = cl, k = k, agents = paste(run, collapse = "+"), rho = rho,
         m_eff_per_colony = if (k >= 2) k / (1 + (k - 1) * rho) else 1)
})
cat("\nDesign effect per class\n"); print(as.data.frame(deff_tab), digits = 3, row.names = FALSE)

rel_components <- map_dfr(CLASSES, function(cl) {
  key <- CLASS_KEY[[cl]]
  p <- vis[[paste0(key, "_prop_resistant_conservative")]]
  m <- vis[[paste0(key, "_n_tested")]]
  ok <- !is.na(p) & !is.na(m) & m >= 2
  pp <- p[ok]; mm <- m[ok]
  mu <- mean(pp); Vobs <- var(pp); mbar <- mean(mm); Einv <- mean(1 / mm)
  # drug-test frame, as run (1/mbar, clamped at zero as in the pipeline)
  V_run <- max((Vobs - mu * (1 - mu) / mbar) / (1 - 1 / mbar), 0)
  # drug-test frame, E[1/m]
  V_dt  <- (Vobs - mu * (1 - mu) * Einv) / (1 - Einv)
  # colony frame, E[1/m]
  cf <- depth %>% filter(class == cl, colonies >= 2) %>%
    mutate(q = colonies_ns / colonies)
  mu_c <- mean(cf$q); Einv_c <- mean(1 / cf$colonies)
  V_col <- (var(cf$q) - mu_c * (1 - mu_c) * Einv_c) / (1 - Einv_c)
  tibble(class = cl, visits_m2 = sum(ok), visits_colony_m2 = nrow(cf),
         mu = mu, V_obs = Vobs, mean_m = mbar, E_inv_m = Einv,
         V_as_run = V_run, E_as_run = mu * (1 - mu) - V_run,
         V_dt = V_dt, E_dt = mu * (1 - mu) - V_dt,
         V_col = V_col, E_col = mu_c * (1 - mu_c) - V_col)
})

# The as-run components must reproduce the pipeline's CSV exactly, or this
# script is not drawing what the pipeline computed.
chk <- read_csv("results/design/depth_reliability.csv", show_col_types = FALSE) %>%
  inner_join(rel_components, by = "class")
stopifnot(nrow(chk) == 4,
          max(abs(chk$V_theta - chk$V_as_run)) < 1e-10,
          max(abs(chk$E_theta1theta - chk$E_as_run)) < 1e-10)
cat("  as-run components reproduce results/design/depth_reliability.csv\n")

m_eff_ag <- deff_tab$m_eff_per_colony[deff_tab$class == "Aminoglycoside"]

rel_setup <- function(method) {
  rel_components %>%
    left_join(deff_tab %>% select(class, m_eff_per_colony), by = "class") %>%
    transmute(class,
              V = switch(method, as_run = V_as_run, class_deff = V_dt, colony = V_col),
              E = switch(method, as_run = E_as_run, class_deff = E_dt, colony = E_col),
              m_eff = switch(method, as_run = m_eff_ag,
                             class_deff = m_eff_per_colony, colony = 1))
}
rel_curve <- function(method, colonies = 1:10) {
  rel_setup(method) %>%
    tidyr::expand_grid(colonies = colonies) %>%
    mutate(reliability = V / (V + E / (colonies * m_eff)), method = method)
}
colonies_needed <- function(method, target = 0.8) {
  rel_curve(method, 1:40) %>% filter(reliability >= target) %>%
    group_by(class) %>% summarise(colonies_needed = min(colonies), .groups = "drop") %>%
    mutate(method = method)
}

methods_all <- c("as_run", "class_deff", "colony")
need_tab <- map_dfr(methods_all, colonies_needed) %>%
  pivot_wider(names_from = method, values_from = colonies_needed) %>%
  mutate(class = factor(class, levels = CLASSES)) %>% arrange(class)
rel_at <- map_dfr(methods_all, rel_curve, colonies = 1:4) %>%
  select(method, class, colonies, reliability) %>%
  pivot_wider(names_from = colonies, values_from = reliability, names_prefix = "R_")

cat("\nColonies per culture needed for reliability >= 0.8, by method\n")
print(as.data.frame(need_tab), row.names = FALSE)
if (!identical(as.integer(need_tab$as_run),
               as.integer(read_csv("results/design/colonies_for_reliability.csv",
                                   show_col_types = FALSE) %>%
                          mutate(class = factor(class, levels = CLASSES)) %>%
                          arrange(class) %>% pull(colonies_needed))))
  warning("as_run colonies-needed does not match results/design/colonies_for_reliability.csv")
cat("\nReliability at 1-4 colonies\n"); print(as.data.frame(rel_at), digits = 2, row.names = FALSE)

write_csv(rel_components, file.path(out_res, "fig3_reliability_components.csv"))
write_csv(deff_tab,       file.path(out_res, "fig3_design_effect.csv"))
write_csv(need_tab,       file.path(out_res, "fig3_colonies_needed_by_method.csv"))
write_csv(rel_at,         file.path(out_res, "fig3_reliability_at_1to4_colonies.csv"))

cur_colonies <- mean(depth_summary$colonies_mean)

# FQ and cephalosporin curves coincide under class_deff (R = 0.52, 0.68, 0.76,
# 0.81 at 1-4 colonies, to two decimals), so classes get their own point shape
# and a small horizontal offset; otherwise one line hides the other.
CLASS_SHAPE <- c(Aminoglycoside = 16, Carbapenem = 17, Fluoroquinolone = 15, Cephalosporin = 18)

plot_reliability <- function(method) {
  cur  <- rel_curve(method)
  need <- colonies_needed(method)
  labs_cls <- setNames(sprintf("%s (%d)", need$class, need$colonies_needed), need$class)
  cur <- cur %>% mutate(class = factor(class, levels = CLASSES),
                        x = colonies + (as.numeric(class) - 2.5) * 0.05)
  ggplot(cur, aes(x, reliability, colour = class, shape = class)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_vline(xintercept = cur_colonies, linetype = "dotted", colour = "grey30", linewidth = 0.4) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.3) +
    annotate("text", x = cur_colonies + 0.15, y = 0.315, hjust = 0, vjust = 0,
             size = 8 / .pt, family = FONT, colour = "grey20",
             label = sprintf("current mean\n%.1f colonies", cur_colonies)) +
    annotate("text", x = 10, y = 0.81, hjust = 1, vjust = 0,
             size = 8 / .pt, family = FONT, colour = "grey20", label = "target 0.8") +
    scale_x_continuous(breaks = 1:10, expand = expansion(add = 0.3)) +
    scale_y_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1),
                       expand = expansion(mult = c(0, 0.01))) +
    scale_colour_manual(values = CLASS_COL, labels = labs_cls,
                        name = "Colonies for reliability ≥ 0.8") +
    scale_shape_manual(values = CLASS_SHAPE, labels = labs_cls,
                       name = "Colonies for reliability ≥ 0.8") +
    labs(x = "Colonies tested per culture",
         y = "Reliability of the visit resistance frequency") +
    theme_pub() +
    theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
          legend.title = element_text(size = 8.5))
}

# (a) what depth buys clinically: the chance that a culture contains at least
# one resistant colony when a fraction f of the population is resistant,
#     P = 1 - (1 - f)^n,  n = colonies (no design effect: one colony is one draw,
#     however many agents it is tested against).
# Current practice is averaged over the OBSERVED colony distribution, not
# evaluated at the mean depth (Jensen: plugging in 1.5 overstates it).
F_RES <- c(0.1, 0.2, 0.3)
need_rng <- range(need_tab[[REL_METHOD]])
catch <- tidyr::expand_grid(colonies = 1:10, f = F_RES) %>%
  mutate(p = 1 - (1 - f)^colonies)
depth_obs <- depth %>% filter(class == "Aminoglycoside") %>% pull(colonies)
catch_now <- tibble(f = F_RES) %>%
  mutate(p_observed_mix = map_dbl(f, ~ mean(1 - (1 - .x)^depth_obs)),
         p_at_mean_depth = 1 - (1 - f)^cur_colonies,
         p_at_4 = 1 - (1 - f)^4, p_at_5 = 1 - (1 - f)^5)
cat("\nFigure 3a: chance of catching a resistant subpopulation\n")
print(as.data.frame(catch_now), digits = 2, row.names = FALSE)
write_csv(catch_now, file.path(out_res, "fig3a_catch_resistant_now.csv"))

F_COL <- setNames(c("grey70", "grey45", "black"), sprintf("%.0f%% resistant", 100 * F_RES))
catch <- catch %>% mutate(lab = factor(sprintf("%.0f%% resistant", 100 * f), levels = names(F_COL)))
catch_now_pt <- catch_now %>%
  mutate(lab = factor(sprintf("%.0f%% resistant", 100 * f), levels = names(F_COL)))

shade_need <- annotate("rect", xmin = need_rng[1] - 0.3, xmax = need_rng[2] + 0.3,
                       ymin = -Inf, ymax = Inf, fill = "grey92")

fig3a <- ggplot(catch, aes(colonies, p, colour = lab)) +
  shade_need +
  geom_vline(xintercept = cur_colonies, linetype = "dotted", colour = "grey30", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  geom_point(data = catch_now_pt, aes(x = cur_colonies, y = p_observed_mix),
             shape = 21, fill = "white", size = 1.8, stroke = 0.6) +
  scale_colour_manual(values = F_COL, name = "Population\nresistant") +
  scale_x_continuous(breaks = 1:10, expand = expansion(add = 0.3)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     labels = label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Colonies tested per culture",
       y = "Chance a culture contains a resistant colony") +
  theme_pub() +
  theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
        legend.title = element_text(size = 8.5))

reliability_panel <- function(method) {
  plot_reliability(method) +
    annotate("rect", xmin = range(need_tab[[method]])[1] - 0.3,
             xmax = range(need_tab[[method]])[2] + 0.3,
             ymin = -Inf, ymax = Inf, fill = "grey92") -> p
  # move the shading under the data layers
  p$layers <- c(p$layers[length(p$layers)], p$layers[-length(p$layers)])
  p
}
fig3 <- (fig3a | reliability_panel(REL_METHOD)) + plot_annotation(tag_levels = "a")
save_fig(fig3, "Fig3_depth", W2, 78)
for (mth in setdiff(methods_all, REL_METHOD))
  save_fig(reliability_panel(mth), paste0("Fig3b_reliability_", mth), W1, 75,
           dir = file.path(out_fig, "variants"))


# #############################################################################
# FIGURE 4 - within-patient autocorrelation against visit interval
#
# Recomputes AMR_pipeline.R 14.2 exactly: per class, visit-level R + I
# frequency centred on the patient's own mean; consecutive-visit pairs pooled
# over the four classes; Pearson r within interval bins holding >= 60 pairs;
# r(t) = a exp(-t / tau) by nls weighted by pairs. Adds a patient-level
# bootstrap (patients resampled with replacement, all their pairs in all classes
# kept together) for the bin intervals and for a.
# #############################################################################

ac_pairs <- map_dfr(CLASSES, function(cl) {
  col <- paste0(CLASS_KEY[[cl]], "_prop_resistant_conservative")
  vis %>%
    arrange(individual_id, date_collected_from_subject) %>%
    filter(!is.na(.data[[col]])) %>%
    group_by(individual_id) %>%
    mutate(gap = as.numeric(as.Date(date_collected_from_subject) -
                              lag(as.Date(date_collected_from_subject))),
           xc = .data[[col]] - mean(.data[[col]]), lxc = lag(xc)) %>%
    ungroup() %>%
    filter(!is.na(lxc), !is.na(gap)) %>%
    transmute(individual_id, class = cl, gap, xc, lxc)
})
BREAKS <- c(0, 14, 28, 42, 60, 90, 150, 240, 400, Inf)
bin_ac <- function(d) {
  d %>% mutate(bin = cut(gap, BREAKS)) %>%
    group_by(bin) %>% filter(n() >= 60) %>%
    summarise(n = n(), gap_mid = median(gap), r = cor(xc, lxc), .groups = "drop")
}
fit_ac <- function(b) {
  f <- try(nls(r ~ a * exp(-gap_mid / tau), data = b,
               start = list(a = 0.4, tau = 60), weights = b$n), silent = TRUE)
  if (inherits(f, "try-error")) return(c(a = NA_real_, tau = NA_real_))
  coef(f)
}

ac_bins <- bin_ac(ac_pairs)
ac_coef <- fit_ac(ac_bins)
a_hat <- unname(ac_coef["a"]); tau_hat <- unname(ac_coef["tau"])

old <- read_csv("results/design/interval_autocorrelation.csv", show_col_types = FALSE) %>%
  filter(src == "observed")
stopifnot(nrow(old) == nrow(ac_bins), max(abs(old$r - ac_bins$r)) < 1e-10)
cat(sprintf("\nFigure 4: %d pairs, %d bins; reproduces interval_autocorrelation.csv\n",
            nrow(ac_pairs), nrow(ac_bins)))
cat(sprintf("  fit: r(t) = %.3f exp(-t / %.1f d)\n", a_hat, tau_hat))

ids <- unique(ac_pairs$individual_id)
pairs_by_id <- split(ac_pairs, ac_pairs$individual_id)
boot <- map_dfr(seq_len(N_BOOT), function(b) {
  d <- bind_rows(pairs_by_id[as.character(sample(ids, replace = TRUE))])
  bb <- d %>% mutate(bin = cut(gap, BREAKS)) %>%
    filter(bin %in% ac_bins$bin) %>%
    group_by(bin) %>%
    summarise(n = n(), gap_mid = median(gap), r = cor(xc, lxc), .groups = "drop")
  cf <- fit_ac(bb)
  bb %>% mutate(rep = b, a = cf[["a"]], tau = cf[["tau"]])
})
bin_ci <- boot %>% group_by(bin) %>%
  summarise(r_lo = quantile(r, 0.025, na.rm = TRUE),
            r_hi = quantile(r, 0.975, na.rm = TRUE), .groups = "drop")
par_ci <- boot %>% distinct(rep, a, tau) %>%
  summarise(fits_ok = sum(!is.na(a)),
            a_lo = quantile(a, 0.025, na.rm = TRUE), a_hi = quantile(a, 0.975, na.rm = TRUE),
            tau_lo = quantile(tau, 0.025, na.rm = TRUE), tau_hi = quantile(tau, 0.975, na.rm = TRUE))
ac_bins <- ac_bins %>% left_join(bin_ci, by = "bin")
cat(sprintf("  bootstrap (%d reps, %d fits converged): a 95%% CI [%.3f, %.3f], tau [%.0f, %.0f] d\n",
            N_BOOT, par_ci$fits_ok, par_ci$a_lo, par_ci$a_hi, par_ci$tau_lo, par_ci$tau_hi))

gap_median <- median(ac_pairs$gap)
fitted_at <- tibble(days = c(14, 28, 56, 84), r = a_hat * exp(-days / tau_hat))
write_csv(ac_bins, file.path(out_res, "fig4_autocorrelation_bins.csv"))
write_csv(tibble(a = a_hat, tau_days = tau_hat, par_ci, pairs = nrow(ac_pairs),
                 patients = n_distinct(ac_pairs$individual_id),
                 median_interval_days = gap_median),
          file.path(out_res, "fig4_autocorrelation_fit.csv"))
write_csv(fitted_at, file.path(out_res, "fig4_fitted_at_intervals.csv"))

ac_curve <- tibble(gap = seq(0, 620, 2)) %>% mutate(r = a_hat * exp(-gap / tau_hat))

fig4 <- ggplot(ac_bins, aes(gap_mid, r)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = a_hat, linetype = "dashed", colour = "grey35", linewidth = 0.35) +
  geom_vline(xintercept = gap_median, linetype = "dotted", colour = "grey45", linewidth = 0.4) +
  geom_line(data = ac_curve, aes(gap, r), colour = "black", linewidth = 0.55) +
  geom_linerange(aes(ymin = r_lo, ymax = r_hi), colour = "grey40", linewidth = 0.4) +
  geom_point(aes(size = n), shape = 21, fill = "white", colour = "black", stroke = 0.5) +
  annotate("text", x = 620, y = a_hat + 0.012, hjust = 1, vjust = 0,
           size = 8 / .pt, family = FONT, colour = "grey20",
           label = sprintf("ceiling a = %.2f", a_hat)) +
  annotate("text", x = gap_median + 8, y = -0.33, hjust = 0, vjust = 0,
           size = 8 / .pt, family = FONT, colour = "grey20",
           label = sprintf("median interval\n%.0f days", gap_median)) +
  scale_size_area(max_size = 3.2, breaks = c(200, 500, 900), name = "Visit pairs") +
  scale_x_continuous(breaks = seq(0, 600, 100)) +
  scale_y_continuous(breaks = round(seq(-0.3, 0.5, 0.1), 1),
                     labels = label_number(accuracy = 0.1)) +
  # coord limits, not scale limits: scale limits DROP any interval that runs
  # past them (the first and last bins did), coord limits only clip.
  coord_cartesian(xlim = c(0, 620), ylim = c(-0.35, 0.56), expand = FALSE) +
  labs(x = "Interval between consecutive visits (days)",
       y = "Within-patient lag-1 correlation") +
  theme_pub() +
  theme(legend.position = c(0.99, 0.995), legend.justification = c(1, 1),
        legend.direction = "horizontal", legend.title = element_text(size = 8.5),
        plot.margin = margin(6, 8, 4, 4))

save_fig(fig4, "Fig4_autocorrelation_interval", W1, 70)


# #############################################################################
# FIGURE 5 - design frontier
#
# Source: results/design/power_grid_dt.csv (pipeline 14.3, 300 reps per cell).
# Power of the cross-lagged mixed model to detect the reference coupling
# beta = 0.2 at a 28-day interval, on the simulated 3 x 3 grid of patients x
# visits per patient, at 3 and 10 tests per visit. The 80% frontier is the
# 0.8 contour of power interpolated bilinearly in log(patients) x log(visits);
# the grid values themselves are printed on the points so the interpolation
# is visible for what it is.
# #############################################################################

pg <- read_csv("results/design/power_grid_dt.csv", show_col_types = FALSE)
REF_BETA <- 0.2; REF_DT <- 28

f5_pts <- pg %>% filter(beta == REF_BETA, dt == REF_DT) %>%
  mutate(depth = factor(sprintf("%d tests per visit", m),
                        levels = sprintf("%d tests per visit", sort(unique(m)))),
         total_visits = n_pat * L)
write_csv(f5_pts %>% select(m, n_pat, L, total_visits, power, power_se),
          file.path(out_res, "fig5_power_grid_beta0.2_dt28.csv"))

bilinear <- function(xs, ys, Z, x, y) {
  i  <- findInterval(x, xs, rightmost.closed = TRUE, all.inside = TRUE)
  j  <- findInterval(y, ys, rightmost.closed = TRUE, all.inside = TRUE)
  tx <- (x - xs[i]) / (xs[i + 1] - xs[i]); ty <- (y - ys[j]) / (ys[j + 1] - ys[j])
  (1 - tx) * (1 - ty) * Z[cbind(i, j)] + tx * (1 - ty) * Z[cbind(i + 1, j)] +
    (1 - tx) * ty * Z[cbind(i, j + 1)] + tx * ty * Z[cbind(i + 1, j + 1)]
}
f5_surface <- f5_pts %>% group_split(depth) %>% map_dfr(function(d) {
  xs <- log(sort(unique(d$n_pat))); ys <- log(sort(unique(d$L)))
  Z  <- d %>% arrange(n_pat, L) %>% pull(power) %>%
    matrix(nrow = length(xs), byrow = TRUE)
  g <- tidyr::expand_grid(lx = seq(min(xs), max(xs), length.out = 150),
                          ly = seq(min(ys), max(ys), length.out = 150))
  g %>% mutate(power = bilinear(xs, ys, Z, lx, ly),
               n_pat = exp(lx), L = exp(ly), depth = d$depth[1])
})

# Lines of constant total visits (patients x visits per patient), the same in
# both panels. Labelled once, in the right panel, in the gaps between grid rows
# and columns: in the left panel the frontier crosses every gap a label fits in.
iso_visits <- tidyr::expand_grid(V = c(1000, 2000, 4000),
                                 n_pat = exp(seq(log(20), log(300), length.out = 80))) %>%
  mutate(L = V / n_pat)
iso_lab <- tibble(V = c(1000, 2000, 4000), n_pat = c(45, 52, 134)) %>%
  mutate(L = V / n_pat, lab = paste(comma(V), "visits"),
         depth = factor("10 tests per visit", levels = levels(f5_pts$depth)))

current <- tibble(n_pat = n_distinct(vis$individual_id), L = 12,
                  depth = factor("3 tests per visit", levels = levels(f5_pts$depth)))

x_lim <- c(24, 250); y_lim <- c(10, 56)
fig5 <- ggplot(f5_pts, aes(n_pat, L)) +
  geom_line(data = iso_visits, aes(group = V), colour = "grey78", linewidth = 0.3) +
  geom_label(data = iso_lab, aes(label = lab), size = 8 / .pt, family = FONT,
             colour = "grey45", fill = "white", border.colour = NA,
             label.padding = unit(0.5, "pt")) +
  geom_contour(data = f5_surface, aes(z = power, linetype = "80% power frontier"),
               breaks = 0.8, colour = "black", linewidth = 0.6) +
  geom_point(aes(fill = power), shape = 21, size = 8.6, colour = "black", stroke = 0.3) +
  geom_text(aes(label = sub("^0", "", sprintf("%.2f", power)),
                colour = power > 0.6), size = 8 / .pt, family = FONT) +
  geom_point(data = current, shape = 23, size = 2.6, fill = "#D55E00",
             colour = "black", stroke = 0.3, position = position_nudge(x = 0.13)) +
  geom_text(data = current, label = "current\ncohort", hjust = 0.5, vjust = 0, lineheight = 0.9,
            size = 8 / .pt, family = FONT, colour = "grey20",
            position = position_nudge(x = 0.13, y = 0.035)) +
  facet_wrap(~depth, nrow = 1) +
  # rectangles, not the default raster: the Quartz PDF device draws raster
  # images upside down, which inverted the colour bar against its labels.
  scale_fill_gradient(low = "white", high = "#08306B", limits = c(0, 1),
                      breaks = c(0, 0.4, 0.8), name = "Power",
                      guide = guide_colourbar(display = "rectangles", nbin = 60, order = 1)) +
  scale_linetype_manual(values = "solid", name = NULL, guide = guide_legend(order = 2)) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_x_log10(breaks = c(30, 90, 200)) +
  scale_y_log10(breaks = c(12, 24, 48),
                sec.axis = sec_axis(~ . * REF_DT / 365.25, breaks = c(1, 2, 3),
                                    name = sprintf("Follow-up at %d-day interval (years)", REF_DT))) +
  coord_cartesian(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  labs(x = "Patients", y = "Visits per patient") +
  theme_pub() +
  theme(legend.position = "right", legend.title = element_text(size = 8.5),
        legend.key.height = unit(10, "pt"), panel.spacing.x = unit(10, "pt"))

save_fig(fig5, "Fig5_design_frontier", W2, 88)


# #############################################################################
# SUPPLEMENTARY FIGURES
#
#   S1  Cohort structure: patient spans, patients seen per half-year, visits per
#       patient, intervals by the clinical status of the later visit
#   S2  Agent coverage: results per agent across isolates
#   S3  Every patient's trajectory, by class
#   S4  Identifiability of the between-visit SD       (needs analyses/fit_...)
#   S5  Prior sensitivity of the variance SDs          (needs analyses/fit_...)
#   S6  Colonies needed for reliability 0.8, by estimator
#   S7  Autocorrelation against interval, by exacerbation structure
#   S8  Multispatial CCM: skill against library size, and against correlation
#   S9  Morphotype composition at two-colony visits against binomial picking
#   S10 Class resistance frequency by year of follow-up
# Numbered in order of first citation in the manuscript; code order differs.
#
# Time since a patient's first visit is recomputed from the filtered dates:
# days_since_first_visit is wrong for 26 patients (README open issue 5).
# #############################################################################

out_sup <- file.path(out_fig, "supplement")
dir.create(out_sup, showWarnings = FALSE)
save_sup <- function(p, name, w, h) save_fig(p, name, w, h, dir = out_sup)

vis_t <- vis %>%
  mutate(date = as.Date(date_collected_from_subject)) %>%
  group_by(individual_id) %>%
  mutate(first_date = min(date), years = as.numeric(date - first_date) / 365.25) %>%
  ungroup()

# ---- S1 cohort structure ----------------------------------------------------
spans <- vis_t %>% group_by(individual_id) %>%
  summarise(first = min(date), last = max(date), n_visits = n(), .groups = "drop") %>%
  arrange(first, last) %>% mutate(rank = row_number())
s1a <- ggplot(spans) +
  geom_segment(aes(x = first, xend = last, y = rank, yend = rank),
               colour = "grey60", linewidth = 0.25) +
  geom_point(data = vis_t %>% left_join(spans %>% select(individual_id, rank), by = "individual_id"),
             aes(date, rank), size = 0.25, colour = "black") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = 0.01)) +
  scale_y_reverse(expand = expansion(add = 1)) +
  labs(x = NULL, y = "Patients (by first visit)") +
  theme_pub() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

half_year <- vis_t %>%
  mutate(window = as.Date(paste0(format(date, "%Y"), ifelse(as.integer(format(date, "%m")) <= 6, "-01-01", "-07-01")))) %>%
  group_by(window) %>% summarise(patients = n_distinct(individual_id), visits = n(), .groups = "drop")
s1b <- ggplot(half_year, aes(window + 91, patients)) +
  geom_col(width = 160, fill = "grey55") +
  geom_text(aes(label = patients), vjust = -0.4, size = 8 / .pt, family = FONT) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = 0.01),
               limits = range(c(spans$first, spans$last)) + c(-30, 30)) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 20), expand = expansion(mult = c(0, 0))) +
  labs(x = "Collection date", y = "Patients seen\nper half-year") +
  theme_pub()

s1c <- ggplot(spans, aes(n_visits)) +
  geom_histogram(binwidth = 2, boundary = 0.5, fill = "grey55", colour = "white", linewidth = 0.2) +
  labs(x = "Cultures per patient", y = "Patients") + theme_pub()

gaps <- vis_t %>% arrange(individual_id, date) %>% group_by(individual_id) %>%
  mutate(gap = as.numeric(date - lag(date))) %>% ungroup() %>% filter(!is.na(gap)) %>%
  mutate(status = factor(coalesce(as.character(clinical_status), "Not recorded"),
                         levels = c("APE IV Antibiotics", "APE Outpatient",
                                    "Clinically Stable", "Not recorded"),
                         labels = c("Exacerbation, IV", "Exacerbation, outpatient",
                                    "Clinically stable", "Not recorded")))
gap_breaks <- c(0, 14, 28, 42, 60, 90, 150, 240, 400, Inf)
gap_tab <- gaps %>% mutate(bin = cut(gap, gap_breaks, labels = c("1–14", "15–28", "29–42",
                                                                  "43–60", "61–90", "91–150",
                                                                  "151–240", "241–400", ">400"))) %>%
  count(bin, status)
write_csv(gap_tab, file.path(out_res, "figS1_interval_by_status.csv"))
s1d <- ggplot(gap_tab, aes(bin, n, fill = status)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("#D55E00", "#E69F00", "#56B4E9", "grey80")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Days since the patient's previous culture", y = "Cultures",
       fill = "Status at the later culture") +
  theme_pub() +
  guides(fill = guide_legend(nrow = 2, title.position = "top")) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top",
        legend.title = element_text(size = 8.5), legend.margin = margin(0, 0, 0, 0))

figS1 <- (s1a / s1b / (s1c | s1d)) +
  plot_layout(heights = c(2.2, 0.8, 1.3)) + plot_annotation(tag_levels = "a")
save_sup(figS1, "FigS1_cohort_structure", W2, 175)
cat(sprintf("  S1: median gap %.0f d; %.0f%% of gaps > 28 d\n", median(gaps$gap), 100 * mean(gaps$gap > 28)))

# ---- S2 agent coverage ------------------------------------------------------
AGENTS <- tribble(
  ~agent,                     ~label,                      ~cls,
  "amikacin",                 "Amikacin",                  "Aminoglycoside",
  "tobramycin",               "Tobramycin",                "Aminoglycoside",
  "gentamicin",               "Gentamicin",                "Aminoglycoside",
  "meropenem",                "Meropenem",                 "Carbapenem",
  "imipenem",                 "Imipenem",                  "Carbapenem",
  "levofloxacin",             "Levofloxacin",              "Fluoroquinolone",
  "ciprofloxacin",            "Ciprofloxacin",             "Fluoroquinolone",
  "cefepime",                 "Cefepime",                  "Cephalosporin",
  "ceftazidime",              "Ceftazidime",               "Cephalosporin",
  "piperacillin_tazobactam",  "Piperacillin–tazobactam", "Other",
  "aztreonam",                "Aztreonam",                 "Other",
  "colistin",                 "Colistin",                  "Other")
cov_tab <- AGENTS %>%
  mutate(n = map_int(agent, ~ sum(iso[[.x]] %in% SIR)), pct = 100 * n / nrow(iso),
         cls = factor(cls, levels = c(CLASSES, "Other")),
         label = factor(label, levels = rev(label)))
write_csv(cov_tab, file.path(out_res, "figS2_agent_coverage.csv"))
figS2 <- ggplot(cov_tab, aes(pct, label, fill = cls)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = comma(n)), hjust = -0.15, size = 8 / .pt, family = FONT) +
  facet_grid(cls ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c(CLASS_COL, Other = "grey70"), guide = "none") +
  scale_x_continuous(limits = c(0, 118), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = sprintf("Isolates with a result (%% of %s)", comma(nrow(iso))), y = NULL) +
  theme_pub() +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5),
        panel.spacing.y = unit(3, "pt"))
save_sup(figS2, "FigS2_agent_coverage", W15, 95)

# ---- S3 morphotype composition at two-colony visits -------------------------
phe <- iso %>% filter(!is.na(pseudomonas_aeruginosa_phenotype)) %>%
  mutate(visit_date = as.Date(date_collected_from_subject),
         mucoid = pseudomonas_aeruginosa_phenotype == "mucoid")
q_patient <- phe %>% group_by(individual_id) %>% summarise(q_pat = mean(mucoid), .groups = "drop")
two <- phe %>% group_by(individual_id, visit_date) %>%
  filter(n() == 2) %>% summarise(k = sum(mucoid), .groups = "drop") %>%
  left_join(q_patient, by = "individual_id")
q_pool <- mean(two$k) / 2
obs <- as.numeric(table(factor(two$k, levels = 0:2)))
exp_pool <- nrow(two) * dbinom(0:2, 2, q_pool)
exp_pat  <- c(sum((1 - two$q_pat)^2), sum(2 * two$q_pat * (1 - two$q_pat)), sum(two$q_pat^2))
chi_pool <- sum((obs - exp_pool)^2 / exp_pool)
p_pool   <- pchisq(chi_pool, df = 1, lower.tail = FALSE)     # 3 cells - 1 - 1 estimated q
# patient-specific null by simulation: is the number of mixed visits this high?
sim_mixed <- replicate(10000, sum(rbinom(nrow(two), 2, two$q_pat) == 1))
p_pat <- (1 + sum(sim_mixed >= obs[2])) / (1 + length(sim_mixed))
# fraction of colonies that are mucoid, by colonies phenotyped at the visit
depth_muc <- phe %>% group_by(individual_id, visit_date) %>%
  summarise(nc = n(), k = sum(mucoid), .groups = "drop") %>%
  mutate(depth = factor(pmin(nc, 4), levels = 1:4, labels = c("1", "2", "3", "≥4"))) %>%
  group_by(depth) %>%
  summarise(visits = n(), colonies = sum(nc), mucoid = sum(k), .groups = "drop") %>%
  mutate(frac = mucoid / colonies,
         lo = qbeta(0.025, mucoid + 0.5, colonies - mucoid + 0.5),
         hi = qbeta(0.975, mucoid + 0.5, colonies - mucoid + 0.5))
morph <- tibble(composition = factor(c("Both non-mucoid", "One of each", "Both mucoid"),
                                     levels = c("Both non-mucoid", "One of each", "Both mucoid")),
                Observed = obs, `Random picking, pooled` = exp_pool,
                `Random picking, patient-specific` = exp_pat)
write_csv(morph %>% mutate(q_pooled = q_pool, chisq = chi_pool, p_pooled = p_pool,
                           p_patient_specific_sim = p_pat, visits = nrow(two)),
          file.path(out_res, "figS9_morphotype_two_colony.csv"))
write_csv(depth_muc, file.path(out_res, "figS9_mucoid_fraction_by_depth.csv"))
cat(sprintf("  S9: %d two-colony visits; mixed observed %d vs %.0f pooled / %.0f patient-specific; chi2 p = %.2g; sim p = %.2g\n",
            nrow(two), obs[2], exp_pool[2], exp_pat[2], p_pool, p_pat))

s3a <- morph %>% pivot_longer(-composition, names_to = "series", values_to = "visits") %>%
  mutate(series = factor(series, levels = c("Observed", "Random picking, pooled",
                                            "Random picking, patient-specific"))) %>%
  ggplot(aes(composition, visits, fill = series)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("black", "grey60", "grey85"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Morphotypes of the two colonies", y = "Two-colony visits") +
  theme_pub() + theme(legend.position = "top", legend.direction = "vertical")
s3b <- ggplot(depth_muc, aes(depth, frac)) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.5) +
  geom_point(size = 1.6) +
  geom_text(aes(y = hi, label = comma(colonies)), vjust = -0.6, size = 8 / .pt, family = FONT) +
  scale_y_continuous(limits = c(0, 1), labels = label_percent(accuracy = 1)) +
  labs(x = "Colonies phenotyped at the visit", y = "Colonies that are mucoid") +
  theme_pub()
save_sup((s3a | s3b) + plot_layout(widths = c(1.5, 1)) + plot_annotation(tag_levels = "a"),
         "FigS9_morphotype_picking", W2, 70)

# ---- S4 identifiability of sd_visit ----------------------------------------
sim_file <- "results/sensitivity/sd_visit_recovery.csv"
obs_sdv <- f2b %>% filter(grepl("^All visits", subset)) %>% slice(1)
if (file.exists(sim_file)) {
  simr <- read_csv(sim_file, show_col_types = FALSE)
  figS4 <- ggplot(simr, aes(sd_visit_true + (rep - 3) * 0.018, est)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = obs_sdv$sd_visit_lo, ymax = obs_sdv$sd_visit_hi,
             fill = "#CC79A7", alpha = 0.12) +
    geom_hline(yintercept = obs_sdv$sd_visit_median, colour = "#CC79A7", linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = "grey35", linewidth = 0.45) +
    geom_point(aes(shape = rep == 1), size = 1.5, fill = "white", stroke = 0.5) +
    annotate("text", x = 0.84, y = obs_sdv$sd_visit_median, hjust = 1, vjust = -0.5,
             size = 8 / .pt, family = FONT, colour = "#8B3A6B",
             label = sprintf("observed %.2f", obs_sdv$sd_visit_median)) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), breaks = c("TRUE", "FALSE"),
                       labels = c(`TRUE` = "Seed of the original analysis", `FALSE` = "Four further data sets"),
                       name = NULL) +
    scale_x_continuous(breaks = c(0, 0.3, 0.5, 0.8)) +
    coord_cartesian(xlim = c(-0.08, 0.88), ylim = c(0, 1.45), expand = FALSE) +
    labs(x = "Between-visit SD used to simulate (logit scale)",
         y = "Between-visit SD recovered") +
    theme_pub() + theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
  save_sup(figS4, "FigS4_sd_visit_identifiability", W1, 75)
  write_csv(simr %>% group_by(sd_visit_true) %>%
              summarise(reps = n(), median_of_medians = median(est), min_est = min(est),
                        max_est = max(est), mean_lo = mean(lo), mean_hi = mean(hi), .groups = "drop"),
            file.path(out_res, "figS4_sd_visit_recovery_summary.csv"))
} else cat("  S4 skipped: run analyses/fit_identifiability_prior_sensitivity.R first\n")

# ---- S5 prior sensitivity ---------------------------------------------------
prior_file <- "results/sensitivity/prior_sensitivity.csv"
PRIOR_LAB <- c(exp1 = "Exponential(1) (used)", halft = "Half-t(3, 0, 2.5)",
               exp02 = "Exponential(0.2)", exp5_vis = "Exponential(5), visit SD only")
PRIOR_COL <- c("black", "#0072B2", "#E69F00", "#CC79A7")
if (file.exists(prior_file)) {
  ps <- read_csv(prior_file, show_col_types = FALSE) %>%
    mutate(prior = factor(PRIOR_LAB[prior], levels = PRIOR_LAB))
  s5a_dat <- ps %>% filter(model == "AG drug x isolate", stat == "sd") %>%
    mutate(component = factor(component, levels = c("sd_patient", "sd_colony", "sd_visit"),
                              labels = c("Between patients", "Between colonies", "Between visits")))
  s5a <- ggplot(s5a_dat, aes(median, factor(prior, levels = rev(PRIOR_LAB)), colour = prior)) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.5) + geom_point(size = 1.4) +
    facet_wrap(~component, scales = "free_x", nrow = 1) +
    scale_colour_manual(values = PRIOR_COL, guide = "none") +
    expand_limits(x = 0) +
    labs(x = "Posterior SD (logit scale), aminoglycoside drug × isolate model", y = NULL) +
    theme_pub()
  s5b_dat <- ps %>% filter(model == "Table 1", component %in% c("share_patient", "share_visit")) %>%
    mutate(class = factor(class, levels = rev(CLASSES)),
           component = factor(component, levels = c("share_patient", "share_visit"),
                              labels = c("Between patients", "Between visits")),
           y = as.numeric(class) - (as.numeric(prior) - 2.5) * 0.17)   # same top-down order as panel a
  s5b <- ggplot(s5b_dat, aes(median, y, colour = prior)) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.45) + geom_point(size = 1.1) +
    facet_wrap(~component, scales = "free_x", nrow = 1) +
    scale_y_continuous(breaks = seq_along(CLASSES), labels = rev(CLASSES)) +
    scale_colour_manual(values = PRIOR_COL, name = "Prior on the SDs") +
    expand_limits(x = 0) +
    labs(x = "Proportion of total variance, Table 1 models", y = NULL) +
    theme_pub() + theme(legend.position = "bottom", axis.ticks.y = element_blank()) +
    guides(colour = guide_legend(nrow = 2))
  save_sup((s5a / s5b) + plot_layout(heights = c(1, 1.5)) + plot_annotation(tag_levels = "a"),
           "FigS5_prior_sensitivity", W2, 120)
} else cat("  S5 skipped: run analyses/fit_identifiability_prior_sensitivity.R first\n")

# ---- S6 colonies needed, by estimator ---------------------------------------
EST_LAB <- c(as_run = "As first run (1/m̄; aminoglycoside design effect for all)",
             class_deff = "Used (E[1/m]; each class's own design effect)",
             colony = "Colony frame (k of n colonies; no design effect)")
s6a_dat <- need_tab %>% pivot_longer(-class, names_to = "method", values_to = "colonies") %>%
  mutate(method = factor(EST_LAB[method], levels = EST_LAB),
         class = factor(class, levels = rev(CLASSES)),
         y = as.numeric(class) + (as.numeric(method) - 2) * 0.22)
s6a <- ggplot(s6a_dat, aes(colonies, y, shape = method)) +
  geom_point(size = 1.8, fill = "white", stroke = 0.6) +
  geom_text(aes(label = colonies), nudge_x = 0.35, size = 8 / .pt, family = FONT) +
  scale_shape_manual(values = c(21, 16, 24), name = NULL) +
  scale_y_continuous(breaks = seq_along(CLASSES), labels = rev(CLASSES)) +
  scale_x_continuous(limits = c(0, 10.5), breaks = 0:10) +
  labs(x = "Colonies per culture for reliability ≥ 0.8", y = NULL) +
  theme_pub() + theme(legend.position = "top", legend.direction = "vertical",
                      axis.ticks.y = element_blank())
s6b_dat <- rel_components %>%
  transmute(class, `Drug-test frame_Signal` = V_dt, `Drug-test frame_Noise` = E_dt,
            `Colony frame_Signal` = V_col, `Colony frame_Noise` = E_col) %>%
  pivot_longer(-class, names_to = c("frame", "part"), names_sep = "_") %>%
  mutate(class = factor(class, levels = rev(CLASSES)),
         frame = factor(frame, levels = c("Drug-test frame", "Colony frame")),
         part = factor(part, levels = c("Noise", "Signal"),
                       labels = c("Noise, E[θ(1 − θ)]", "Signal, Vθ")))
s6b <- ggplot(s6b_dat, aes(value, class, fill = part)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  facet_wrap(~frame, nrow = 1) +
  scale_fill_manual(values = c("grey70", "grey25"), name = NULL,
                    guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), breaks = c(0, 0.05, 0.1, 0.15),
                     labels = c("0", "", "0.1", "")) +
  labs(x = "Variance component", y = NULL) +
  theme_pub() + theme(legend.position = "top", axis.ticks.y = element_blank(),
                      panel.spacing.x = unit(14, "pt"))
save_sup((s6a | s6b) + plot_layout(widths = c(1, 1.3)) + plot_annotation(tag_levels = "a"),
         "FigS6_colonies_by_estimator", W2, 80)

# ---- S7 autocorrelation by exacerbation structure ---------------------------
is_ape <- function(x) !is.na(x) & x %in% c("APE Outpatient", "APE IV Antibiotics")
ac_cs <- map_dfr(CLASSES, function(cl) {
  col <- paste0(CLASS_KEY[[cl]], "_prop_resistant_conservative")
  vis %>% arrange(individual_id, date_collected_from_subject) %>%
    filter(!is.na(.data[[col]])) %>% group_by(individual_id) %>%
    mutate(gap = as.numeric(as.Date(date_collected_from_subject) -
                              lag(as.Date(date_collected_from_subject))),
           xc = .data[[col]] - mean(.data[[col]]), lxc = lag(xc),
           cs_now = as.character(clinical_status), cs_prev = lag(as.character(clinical_status))) %>%
    ungroup() %>% filter(!is.na(lxc), !is.na(gap)) %>%
    transmute(individual_id, gap, xc, lxc,
              same_episode = is_ape(cs_now) & is_ape(cs_prev) & gap <= 30,
              both_stable = !is.na(cs_now) & !is.na(cs_prev) &
                cs_now == "Clinically Stable" & cs_prev == "Clinically Stable")
})
bin_ac30 <- function(d) d %>% mutate(bin = cut(gap, BREAKS)) %>% group_by(bin) %>%
  filter(n() >= 30) %>% summarise(n = n(), gap_mid = median(gap), r = cor(xc, lxc), .groups = "drop")
SUBSETS <- list(`All pairs` = ac_cs,
                `Excluding same-episode pairs` = ac_cs %>% filter(!same_episode),
                `Both visits clinically stable` = ac_cs %>% filter(both_stable))
s7_bins <- imap_dfr(SUBSETS, ~ bin_ac30(.x) %>% mutate(subset = .y, pairs = nrow(.x)))
s7_fit  <- imap_dfr(SUBSETS, function(d, nm) {
  cf <- fit_ac(bin_ac30(d)); tibble(subset = nm, a = cf[["a"]], tau = cf[["tau"]], pairs = nrow(d))
})
conf_old <- read_csv("results/design/interval_autocorrelation_confound.csv", show_col_types = FALSE)
stopifnot(max(abs(sort(conf_old$a) - sort(s7_fit$a))) < 1e-6)
write_csv(s7_fit, file.path(out_res, "figS7_autocorrelation_by_subset.csv"))
s7_lab <- s7_fit %>% mutate(lab = sprintf("%s (%s), a = %.2f", subset, comma(pairs), a))
s7_curve <- s7_fit %>% tidyr::expand_grid(gap = seq(0, 620, 3)) %>% mutate(r = a * exp(-gap / tau))
SUB_COL <- setNames(c("black", "#E69F00", "#0072B2"), names(SUBSETS))
s7_bins <- s7_bins %>% mutate(subset = factor(subset, levels = names(SUBSETS)))
s7_curve <- s7_curve %>% mutate(subset = factor(subset, levels = names(SUBSETS)))
figS7 <- ggplot(s7_bins, aes(gap_mid, r, colour = subset)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_line(data = s7_curve, aes(gap, r), linewidth = 0.5) +
  geom_point(aes(size = n), shape = 21, fill = "white", stroke = 0.6) +
  scale_colour_manual(values = SUB_COL, labels = setNames(s7_lab$lab, s7_lab$subset), name = NULL) +
  scale_size_area(max_size = 3, guide = "none") +
  scale_y_continuous(breaks = round(seq(-0.3, 0.5, 0.1), 1), labels = label_number(accuracy = 0.1)) +
  coord_cartesian(xlim = c(0, 620), ylim = c(-0.3, 0.5), expand = FALSE) +
  labs(x = "Interval between consecutive visits (days)", y = "Within-patient lag-1 correlation") +
  theme_pub() + theme(legend.position = "top", legend.direction = "vertical",
                      plot.margin = margin(4, 8, 4, 4))
save_sup(figS7, "FigS7_autocorrelation_exacerbation", W1, 85)

# ---- S8 class resistance frequency by year of follow-up ---------------------
yr_long <- map_dfr(CLASSES, function(cl) {
  col <- paste0(CLASS_KEY[[cl]], "_prop_resistant_conservative")
  vis_t %>% filter(!is.na(.data[[col]])) %>%
    transmute(individual_id, year = floor(years), freq = .data[[col]], class = cl)
}) %>% filter(year <= 5)
pat_year <- yr_long %>% group_by(class, year, individual_id) %>%
  summarise(freq = mean(freq), .groups = "drop")
boot_year <- function(d, B = 1000) {
  ids <- unique(d$individual_id)
  est <- replicate(B, { s <- sample(ids, replace = TRUE)
    mean(d$freq[match(s, d$individual_id)]) })
  tibble(mean = mean(d$freq), lo = quantile(est, 0.025), hi = quantile(est, 0.975),
         patients = length(ids))
}
s8_dat <- pat_year %>% group_by(class, year) %>% group_modify(~ boot_year(.x)) %>% ungroup() %>%
  mutate(class = factor(class, levels = CLASSES))
write_csv(s8_dat, file.path(out_res, "figS10_class_frequency_by_year.csv"))
s8_n <- s8_dat %>% filter(class == "Aminoglycoside") %>% select(year, patients)
figS8 <- ggplot(s8_dat, aes(year + (as.numeric(class) - 2.5) * 0.08, mean, colour = class)) +
  geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.45) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.3) +
  geom_text(data = s8_n %>% mutate(lab = ifelse(year == min(year), paste("n =", patients), patients)),
            aes(x = year, y = 0.03, label = lab), inherit.aes = FALSE,
            size = 8 / .pt, family = FONT, colour = "grey30") +
  scale_colour_manual(values = CLASS_COL, name = NULL) +
  scale_x_continuous(breaks = 0:5, labels = c("0–1", "1–2", "2–3", "3–4", "4–5", "5–6")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(-0.4, 5.3), clip = "off") +
  labs(x = "Years since the patient's first culture",
       y = "Mean visit-level resistance frequency (R + I)") +
  theme_pub() + theme(legend.position = "top") + guides(colour = guide_legend(nrow = 2))
save_sup(figS8, "FigS10_resistance_by_year", W1, 80)

# ---- S9 every patient's trajectory ------------------------------------------
traj <- map_dfr(CLASSES, function(cl) {
  col <- paste0(CLASS_KEY[[cl]], "_prop_resistant_conservative")
  vis_t %>% filter(!is.na(.data[[col]])) %>%
    transmute(individual_id, years, freq = .data[[col]], class = cl)
}) %>% mutate(class = factor(class, levels = CLASSES))
# patient-weighted, as in S8: each patient's mean within the half-year, then
# the mean over patients, so frequent attenders do not dominate
traj_mean <- traj %>% mutate(bin = floor(years * 2) / 2 + 0.25) %>%
  group_by(class, bin, individual_id) %>% summarise(freq = mean(freq), .groups = "drop") %>%
  group_by(class, bin) %>% filter(n() >= 10) %>%
  summarise(freq = mean(freq), .groups = "drop")
figS9 <- ggplot(traj, aes(years, freq)) +
  geom_line(aes(group = individual_id), colour = "grey60", linewidth = 0.15, alpha = 0.6) +
  geom_line(data = traj_mean, aes(bin, freq, colour = class), linewidth = 0.8) +
  facet_wrap(~class, nrow = 1) +
  scale_colour_manual(values = CLASS_COL, guide = "none") +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = label_percent(accuracy = 1)) +
  scale_x_continuous(breaks = 0:6) +
  labs(x = "Years since the patient's first culture",
       y = "Visit-level resistance frequency (R + I)") +
  theme_pub() + theme(panel.spacing.x = unit(8, "pt"))
save_sup(figS9, "FigS3_trajectories", W2, 60)

# ---- S10 multispatial CCM ---------------------------------------------------
msr <- read_csv("results/msccm/msccm_results.csv", show_col_types = FALSE)
mss <- read_csv("results/msccm/msccm_summary.csv", show_col_types = FALSE)
flow_dir <- function(pair, direction) ifelse(sub(" ->.*$", "", direction) == sub(" <->.*$", "", pair),
                                             "First class drives second", "Second class drives first")
s10a_dat <- msr %>% filter(!is.na(rho)) %>% mutate(flow = flow_dir(pair, direction),
                                                  pair = gsub(" <-> ", " and\n", pair))
s10a <- ggplot(s10a_dat, aes(L, rho, colour = flow, fill = flow)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = rho - 1.96 * sdev, ymax = rho + 1.96 * sdev), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.45) +
  facet_wrap(~pair, ncol = 3) +
  scale_x_continuous(breaks = c(0, 400, 800)) +
  scale_colour_manual(values = c("#D55E00", "#0072B2"), name = "Cross-map direction") +
  scale_fill_manual(values = c("#D55E00", "#0072B2"), name = "Cross-map direction") +
  coord_cartesian(ylim = c(-0.4, 1)) +
  labs(x = "Library size (visits pooled across patients)", y = "Cross-map skill (ρ)") +
  theme_pub() + theme(legend.position = "top", strip.text = element_text(size = 8.5))
# short labels, defined in the legend
ABBR_CLS <- c(Aminoglycoside = "AG", Carbapenem = "CB", Cephalosporin = "CP", Fluoroquinolone = "FQ")
# static correlation per unordered pair, joined on a sorted key
pair_key <- function(a, b) map2_chr(a, b, ~ paste(sort(c(.x, .y)), collapse = "|"))
ccm_cor <- read_csv("results/msccm/ccm_vs_static_correlation.csv", show_col_types = FALSE) %>%
  transmute(key = pair_key(a, b), static_r)
s10b_dat <- mss %>%
  tidyr::separate(pair, into = c("a", "b"), sep = " <-> ", remove = FALSE) %>%
  mutate(key = pair_key(a, b), flow = flow_dir(pair, direction)) %>%
  left_join(ccm_cor, by = "key") %>%
  mutate(x = static_r + ifelse(flow == "First class drives second", -0.006, 0.006),
         short = paste(ABBR_CLS[a], ABBR_CLS[b], sep = "–"))
stopifnot(!anyNA(s10b_dat$static_r))
s10b <- ggplot(s10b_dat, aes(x, rho_end, colour = flow)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.4) +
  geom_point(size = 1.3) +
  scale_colour_manual(values = c("#D55E00", "#0072B2"), guide = "none") +
  ggrepel::geom_text_repel(data = s10b_dat %>% group_by(short) %>% slice_max(rho_end, n = 1),
                           aes(x = static_r, label = short), colour = "black", size = 8 / .pt,
                           family = FONT, min.segment.length = 0.3, seed = 1, segment.size = 0.2,
                           box.padding = 0.3) +
  coord_equal(xlim = c(0, 0.8), ylim = c(-0.15, 0.8)) +
  labs(x = "Contemporaneous correlation",
       y = "Cross-map skill at the largest library") +
  theme_pub()
save_sup((s10a | s10b) + plot_layout(widths = c(2.2, 1)) + plot_annotation(tag_levels = "a"),
         "FigS8_multispatial_ccm", W2, 95)


# #############################################################################
# NUMBERS FOR THE TEXT
# #############################################################################

cat("\n=== Numbers drawn, for checking against the manuscript ===\n")
cat(sprintf("Fig 1: colonies per culture (by class) mean %.2f, %.0f-%.0f%% single colony;\n",
            cur_colonies, min(depth_summary$pct_one_colony), max(depth_summary$pct_one_colony)))
cat(sprintf("       %d-%d cultures with >= 2 colonies; drug-tests per culture AG %.1f, Carb %.1f, FQ %.1f, Ceph %.1f\n",
            min(depth_summary$n_two_plus), max(depth_summary$n_two_plus),
            depth_summary$tests_mean[1], depth_summary$tests_mean[2],
            depth_summary$tests_mean[3], depth_summary$tests_mean[4]))
cat(sprintf("Fig 3 (%s): colonies for 0.8 = %s\n", REL_METHOD,
            paste(sprintf("%s %d", need_tab$class, need_tab[[REL_METHOD]]), collapse = ", ")))
cat(sprintf("Fig 4: a = %.3f [%.3f, %.3f], tau = %.0f d; fitted r at 14/28/84 d = %.2f/%.2f/%.2f\n",
            a_hat, par_ci$a_lo, par_ci$a_hi, tau_hat,
            fitted_at$r[1], fitted_at$r[2], fitted_at$r[4]))
cat("Figures in", out_fig, "| numbers in", out_res, "\n")
