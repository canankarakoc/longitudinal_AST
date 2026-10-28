# =============================================================================
# AMR_msccm_ecology.R  -  multispatial CCM: ecology -> resistance
#
# Author: Canan Karakoc
#
# WHY THIS FILE, AND WHY TREATMENT IS NOT IN IT
# -----------------------------------------------------------------------------
# Section 13 of AMR_pipeline.R runs msCCM among the four RESISTANCE classes, so
# every edge it can draw is cross-resistance. That is a real result but it is
# not the prediction the project is trying to make. This file asks instead
# whether the co-colonising community drives resistance.
#
# TREATMENT -> RESISTANCE IS DELIBERATELY EXCLUDED. Two reasons, both fatal:
#
#   1. Structural missingness. Every *_prescribed column is 47-48% NA, because
#      prescriptions are only recorded at exacerbation visits. msCCM needs
#      contiguous per-patient series with NA only at patient boundaries; a
#      series that is half missing cannot build lagged embedding vectors.
#
#   2. Coding the NAs to 0 does not rescue it. "Any antibiotic prescribed"
#      under NA -> 0 agrees with "is this an exacerbation visit" 95.5% of the
#      time (486/47/1/544). The treatment series IS the exacerbation series, so
#      any coupling CCM found would be exacerbation -> resistance wearing a
#      different label.
#
# The question "does prior prescription predict resistance" is already answered
# properly by M3-M7 in the pipeline, which adjust for the confounders and do
# not need the series to be dynamically reconstructable. Use those.
#
# Mucoid is also excluded: 58.5% NA, 64 patients, median series length 6.
# (Note the hand-rolled CCM in Section 7 does include mucoid - those edges rest
# on that same 58% missing series and should be read with that in mind.)
#
# What IS usable: has_staph, has_mrsa, has_coinfection, n_coinfections - all
# 0% missing, 92 patients, median 12-visit series.
#
# Runtime: ~4 drivers x 4 classes x 2 directions. Set MSCCM_ITER low to test.
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
source("AMR_common.R")
dir.create("results/msccm", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/msccm", recursive = TRUE, showWarnings = FALSE)
set.seed(42)

MSCCM_ITER <- 100     # bootstrap iterations inside CCM_boot

.dat <- load_amr_data(); pa_visit_full <- .dat$pa_visit_full; rm(.dat)

# -----------------------------------------------------------------------------
# 1  Frame: resistance + ecology on the same visits, one NA between patients
# -----------------------------------------------------------------------------

eco_src <- pa_visit_full %>%
  arrange(individual_id, date_collected_from_subject) %>%
  transmute(individual_id,
            ag_R   = ag_prop_resistant_conservative,
            carb_R = carb_prop_resistant_conservative,
            ceph_R = ceph_prop_resistant_conservative,
            fq_R   = fq_prop_resistant_conservative,
            staph        = as.numeric(has_staph),
            mrsa         = as.numeric(has_mrsa),
            coinfection  = as.numeric(has_coinfection),
            n_coinf      = as.numeric(n_coinfections)) %>%
  as.data.frame()

resist <- c(ag_R = "Aminoglycoside", carb_R = "Carbapenem",
            ceph_R = "Cephalosporin", fq_R = "Fluoroquinolone")
eco    <- c(staph = "S. aureus present", mrsa = "MRSA present",
            coinfection = "Any co-infection", n_coinf = "N co-organisms")

ms_vec <- function(df, col, id = "individual_id") {
  s <- split(df[[col]], df[[id]])
  as.numeric(unlist(lapply(s, function(x) c(x, NA))))
}
ms <- lapply(c(names(resist), names(eco)), function(v) ms_vec(eco_src, v))
names(ms) <- c(names(resist), names(eco))

cat("\n=== 1  series ===\n")
cat("  patients:", n_distinct(eco_src$individual_id),
    "| vector length:", length(ms[[1]]), "\n")
for (v in names(ms))
  cat(sprintf("  %-12s real points %4d | variance %.3f\n",
              v, sum(!is.na(ms[[v]])), var(ms[[v]], na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 2  Embedding dimension per variable
# -----------------------------------------------------------------------------

ms_E <- map_dfr(names(ms), function(v) {
  best <- map_dfr(2:6, function(E) {
    r <- tryCatch(SSR_pred_boot(A = ms[[v]], E = E, predstep = 1, tau = 1),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    tibble(variable = v, E = E, rho = mean(r$rho, na.rm = TRUE))
  })
  if (nrow(best) == 0) return(tibble(variable = v, E = 3L, rho = NA_real_))
  best %>% slice_max(rho, n = 1, with_ties = FALSE)
})
cat("\n=== 2  optimal E ===\n")
print(ms_E %>% mutate(rho = round(rho, 3)) %>% as.data.frame(), row.names = FALSE)
write_csv(ms_E, "results/msccm/eco_embedding.csv")
E_of <- setNames(ms_E$E, ms_E$variable)

# -----------------------------------------------------------------------------
# 3  CCM both directions for every ecology x resistance pair
#
# CCM_boot(A, B): cross-maps A from B's manifold. Convergence of that means
# A DRIVES B. So "ecology -> resistance" is CCM_boot(A = ecology, B = resistance).
# -----------------------------------------------------------------------------

eco_cache <- "results/msccm/eco_ccm_full.csv"
if (file.exists(eco_cache)) {
  cat("\n=== 3  reusing", eco_cache, "(delete it to recompute) ===\n")
  res <- read_csv(eco_cache, show_col_types = FALSE)
} else {
cat("\n=== 3  running CCM:", length(eco), "x", length(resist), "pairs x 2 directions ===\n")
res <- map_dfr(names(eco), function(e) {
  map_dfr(names(resist), function(r) {
    ab <- tryCatch(CCM_boot(A = ms[[e]], B = ms[[r]], E = E_of[[e]],
                            tau = 1, iterations = MSCCM_ITER), error = function(x) NULL)
    ba <- tryCatch(CCM_boot(A = ms[[r]], B = ms[[e]], E = E_of[[r]],
                            tau = 1, iterations = MSCCM_ITER), error = function(x) NULL)
    if (is.null(ab) || is.null(ba)) {
      cat("  skipped", e, "<->", r, "\n"); return(NULL) }
    tt <- tryCatch(ccmtest(ab, ba), error = function(x) c(NA_real_, NA_real_))
    cat(sprintf("  %-12s <-> %-16s rho_end %.3f / %.3f\n", e, r,
                tail(ab$rho, 1), tail(ba$rho, 1)))
    bind_rows(
      tibble(driver = unname(eco[e]), affected = unname(resist[r]),
             direction = paste(unname(eco[e]), "->", unname(resist[r])),
             E = E_of[[e]], L = ab$Lobs, rho = ab$rho, sdev = ab$sdevrho,
             p_converge = unname(tt[1])),
      tibble(driver = unname(resist[r]), affected = unname(eco[e]),
             direction = paste(unname(resist[r]), "->", unname(eco[e])),
             E = E_of[[r]], L = ba$Lobs, rho = ba$rho, sdev = ba$sdevrho,
             p_converge = unname(tt[2])))
  })
})
write_csv(res, eco_cache)
}

# -----------------------------------------------------------------------------
# 4  Summary + the diagnostic that matters: does skill beat plain correlation?
#
# Section 13 found cross-map skill flat from the smallest library and equal to
# the contemporaneous correlation - i.e. CCM recovering association, not
# dynamics. The same test has to be applied here before any edge is believed.
# -----------------------------------------------------------------------------

summ <- res %>% group_by(driver, affected, direction, E) %>%
  summarise(L_min = min(L), L_max = max(L),
            rho_start = rho[which.min(L)], rho_end = rho[which.max(L)],
            delta = rho_end - rho_start,
            p_converge = first(p_converge), .groups = "drop")

static_r <- map_dfr(names(eco), function(e) map_dfr(names(resist), function(r)
  tibble(a = unname(eco[e]), b = unname(resist[r]),
         static = suppressWarnings(cor(eco_src[[e]], eco_src[[r]],
                                       use = "pairwise.complete.obs")))))

out <- summ %>%
  mutate(key = paste(pmin(driver, affected), pmax(driver, affected))) %>%
  left_join(static_r %>% mutate(key = paste(pmin(a, b), pmax(a, b))) %>%
              select(key, static), by = "key") %>%
  mutate(beats_static = rho_end > abs(static),
         converges = !is.na(p_converge) & p_converge < 0.05 & delta > 0)

cat("\n=== 4  SUMMARY ===\n")
print(out %>% select(direction, E, rho_start, rho_end, delta, p_converge,
                     static, converges, beats_static) %>%
        mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(), row.names = FALSE)
write_csv(out, "results/msccm/eco_ccm_summary.csv")

cat("\n  converging edges (p < 0.05 and rho rising):", sum(out$converges, na.rm = TRUE),
    "of", nrow(out), "\n")
cat("  edges whose terminal skill EXCEEDS the plain correlation:",
    sum(out$beats_static, na.rm = TRUE), "\n")
cat("  An edge is only interesting if BOTH are true. Convergence alone can be\n")
cat("  produced by autocorrelation; beating the contemporaneous correlation is\n")
cat("  what distinguishes dynamical coupling from co-occurrence.\n")

# -----------------------------------------------------------------------------
# 5  Figures
# -----------------------------------------------------------------------------

fig_conv <- res %>%
  filter(driver %in% unname(eco)) %>%
  ggplot(aes(L, rho, colour = driver)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ affected, nrow = 1) +
  scale_colour_manual(values = unname(pal[c("red", "blue", "gold", "plum")]), name = NULL) +
  labs(x = "Library size (L)", y = "Cross-map skill (rho)",
       title = "Does the co-colonising community predict resistance?") +
  theme_talk(base_size = 11) +
  theme(legend.position = "top", plot.subtitle = element_blank(),
        plot.caption = element_blank(), strip.text = element_text(size = 10, face = "bold"))
ggsave("figures/msccm/fig_eco_ccm_convergence.pdf", fig_conv,
       width = 11, height = 4.2, dpi = 300)

# A scatter against the identity line was the wrong design here. In the
# resistance-class version rho ran 0.18-0.62 and tracked r, so the diagonal
# earned its place. Here rho is pinned at zero (-0.024 to 0.07, several
# negative) while r reaches only 0.16, so "below the diagonal" is trivially
# true and the reference line does no work. A forest at zero says the same
# thing without inviting anyone to read positions against a line.
fc <- res %>% filter(driver %in% unname(eco)) %>%
  group_by(driver, affected) %>%
  filter(L == max(L)) %>%
  summarise(rho = first(rho), sdev = first(sdev), .groups = "drop") %>%
  mutate(lo = rho - 1.96 * sdev, hi = rho + 1.96 * sdev,
         row = paste(driver, "->", affected),
         crosses0 = lo <= 0 & hi >= 0)

cat("\n=== 5  cross-map skill at largest library, with 95% bootstrap CI ===\n")
print(fc %>% select(row, rho, lo, hi, crosses0) %>%
        mutate(across(where(is.numeric), ~round(., 3))) %>% as.data.frame(), row.names = FALSE)
cat(sprintf("  %d of %d intervals cross zero\n", sum(fc$crosses0), nrow(fc)))
write_csv(fc, "results/msccm/eco_ccm_forest.csv")

fig_vs <- fc %>%
  mutate(row = fct_reorder(row, rho)) %>%
  ggplot(aes(rho, row, colour = driver)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 1.1) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = unname(pal[c("red", "blue", "gold", "plum")]), name = NULL) +
  labs(x = "Cross-map skill at largest library (rho, 95% CI)", y = NULL,
       title = "No ecological driver shows coupling detectable at this design") +
  theme_talk(base_size = 11) +
  theme(legend.position = "top", plot.subtitle = element_blank(),
        plot.caption = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
ggsave("figures/msccm/fig_eco_ccm_forest.pdf", fig_vs, width = 8.5, height = 6.5, dpi = 300)

cat("\n  CAVEATS to carry with this figure:\n")
cat("  1. Every static correlation here is |r| < 0.16, so this cannot separate\n")
cat("     'there is no coupling' from 'there is coupling too small to see'.\n")
cat("     Say 'no coupling detectable at this sampling design', not 'no coupling'.\n")
cat("  2. Three of the four drivers are BINARY. Delay-embedding a binary series\n")
cat("     is near-degenerate - the attractor is a lattice of 2^E points and\n")
cat("     neighbours are exact ties - which is an alternative explanation for\n")
cat("     rho ~ 0 that has nothing to do with biology.\n")

cat("\nWrote results/msccm/eco_ccm_summary.csv and figures/msccm/*.pdf\n")
