# =============================================================================
# analyses/export_figures.R  -  manuscript figures under their manuscript numbers
#
# The drawing scripts name figures by content (AMR_figures_manuscript.R ->
# figures/manuscript/, analyses/figures_new.R -> figures/<analysis>/). Figure
# numbers change between drafts, so the numbering lives in ONE table here.
# Edit DRAFT and the table when the draft is renumbered, then rerun.
#
# Output: figures/submission/<Fig#>_<name>.pdf and .png
#
# Run from the repo root:  Rscript analyses/export_figures.R
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")

DRAFT <- "20261019_AST_measurement_manuscript_tracked_v6"

FIGS <- read.csv(text = "
number,name,source
Fig1,sampling_depth,figures/manuscript/Fig1_sampling_depth
Fig2,morphotype_picking,figures/manuscript/supplement/FigS9_morphotype_picking
Fig3,colony_disagreement,figures/sir_matrix/sir_disagreement_matrix
Fig4,variance_components,figures/manuscript/Fig2_variance_components
Fig5,reliability,figures/reliability_split/reliability_split
Fig6,pseudo_trajectories,figures/pseudo_trajectories/pseudo_trajectories
Fig7,power_by_lever,figures/power_recalibrated/power_frontiers
FigS1,cohort_structure,figures/manuscript/supplement/FigS1_cohort_structure
FigS2,agent_coverage,figures/manuscript/supplement/FigS2_agent_coverage
FigS3,trajectories,figures/manuscript/supplement/FigS3_trajectories
FigS4,sd_visit_identifiability,figures/manuscript/supplement/FigS4_sd_visit_identifiability
FigS5,prior_sensitivity,figures/manuscript/supplement/FigS5_prior_sensitivity
FigS6,reliability_estimators,figures/manuscript/supplement/FigS6_colonies_by_estimator
FigS7,autocorrelation,figures/manuscript/supplement/FigS7_autocorrelation_exacerbation
FigS8,multispatial_ccm,figures/manuscript/supplement/FigS8_multispatial_ccm
FigS9,calendar_effects,figures/calendar/calendar_effects
FigS10,design_by_estimand,figures/design_estimand/design_by_estimand
FigS11,deep_downsampling,figures/deep_downsampling/deep_downsampling
", stringsAsFactors = FALSE)

out <- "figures/submission"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
old <- list.files(out, pattern = "\\.(pdf|png)$", full.names = TRUE)
if (length(old)) file.remove(old)           # numbering may have changed: start clean

for (i in seq_len(nrow(FIGS))) {
  for (ext in c("pdf", "png")) {
    src <- paste0(FIGS$source[i], ".", ext)
    if (!file.exists(src)) stop("missing ", src, " - run the drawing scripts first")
    file.copy(src, file.path(out, sprintf("%s_%s.%s", FIGS$number[i], FIGS$name[i], ext)), overwrite = TRUE)
  }
}
writeLines(c(sprintf("Figures numbered as in %s.", DRAFT), "",
             sprintf("%-7s %-26s <- %s", FIGS$number, FIGS$name, FIGS$source)),
           file.path(out, "MAPPING.txt"))
cat(sprintf("exported %d figures to %s (numbering: %s)\n", nrow(FIGS), out, DRAFT))
