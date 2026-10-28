# longitudinal_AST

Code and outputs for a study of what routine antimicrobial susceptibility testing (AST) can resolve in chronic *Pseudomonas aeruginosa* infection in cystic fibrosis. It asks two questions. How well do the one or two colonies tested per sputum culture represent a heterogeneous infection? And, given that sampling, what can routine records tell us about resistance over time?

The manuscript is in preparation. Figures are numbered as in draft v6 (19 October 2026); `figures/submission/MAPPING.txt` lists the source of each.

## Data availability

The clinical data are 1,793 *P. aeruginosa* isolates from 1,078 sputum cultures in 92 people with cystic fibrosis followed at one CF center, 2017–2022. **They are not in this repository.** They were collected under an IRB-approved protocol and cannot be shared publicly. They are available from the corresponding author, subject to IRB approval and a data use agreement. Fitted models and patient-level intermediate files are excluded for the same reason.

The deep-sampling data (75 isolates per sputum, four sputa) are the supplement of Vanderwoude et al. 2024, *mBio* (CC BY 4.0), [doi:10.1128/mbio.03109-23](https://doi.org/10.1128/mbio.03109-23). To run the deep-sample analyses, download `mbio.03109-23-s0008.xlsx` and place it in the repository root.

## Layout

| Path | Contents |
| --- | --- |
| `analytical_data_prep.R` | Cleans the raw export into `data/` |
| `data_inspection_summary.R` | Cohort description |
| `AMR_common.R` | Packages, palettes and shared helpers |
| `AMR_pipeline.R` | Main analysis: variance decomposition, reliability, autocorrelation, cross-mapping. Also holds exploratory sections not reported in the manuscript |
| `AMR_sir_ordinal.R` | Analyses using the full S/I/R call |
| `AMR_msccm_ecology.R` | Multispatial convergent cross-mapping |
| `AMR_figures_manuscript.R` | Figures 1, 2, 4 and S1–S8 (reads `data/` and `results/`) |
| `analyses/` | Supporting fits and the analyses added in revision (below), the figure scripts, and `export_figures.R` |
| `figures/submission/` | Every manuscript figure (PDF and PNG) under its manuscript number |
| `results/` | The numbers behind the figures and tables (aggregated; no patient-level rows) |
| `run_pipeline.sh` | Runs prep, summary and the main pipeline in order |

## Figures and tables

| Item | Analysis script | Results |
| --- | --- | --- |
| Fig 1, Fig 4, Table 1, Figs S1–S8 | `AMR_pipeline.R`, `AMR_sir_ordinal.R`, `AMR_msccm_ecology.R`, `analyses/fit_identifiability_prior_sensitivity.R` | `results/manuscript/` |
| Fig 2 | `AMR_figures_manuscript.R` (morphotype picking) | `results/manuscript/figS9_*` |
| Fig 3 | `analyses/fit_sir_disagreement_matrix.R` | `results/sir_matrix/` |
| Fig 5 | `analyses/fit_reliability_within_patient.R` | `results/reliability_split/` |
| Fig 6 | `analyses/fit_pseudo_trajectories.R` | `results/pseudo_trajectories/` |
| Fig 7, Table 3 | `analyses/fit_power_calibration.R`, `analyses/fit_power_frontiers.R` | `results/power_recalibrated/` |
| Table 2, Fig S10 | `analyses/fit_design_by_estimand.R` | `results/design_estimand/` |
| Fig S9, Table S3 | `analyses/fit_calendar_effects.R`, `analyses/fit_calendar_ou_null.R` | `results/calendar/` |
| Fig S11 | `analyses/fit_deep_downsampling.R` | `results/deep_downsampling/` |
| Table S1 | `analyses/convergence_table.R` | `results/manuscript/tableS1_convergence.csv` |

Result files are named by content, not by figure number: `results/manuscript/figS9_*` is the morphotype analysis now shown as Fig 2.

## Running

The scripts need the data in `data/` (see `data/README.md`). Run each from the repository root, one at a time: two brms jobs at once exhaust 24 GB of memory.

```bash
./run_pipeline.sh all                                   # prep, summary, main pipeline (hours)
Rscript AMR_sir_ordinal.R
Rscript AMR_msccm_ecology.R
Rscript analyses/fit_identifiability_prior_sensitivity.R
Rscript analyses/convergence_table.R
Rscript analyses/fit_reliability_within_patient.R       # ~20 min (patient bootstrap)
Rscript analyses/fit_sir_disagreement_matrix.R
Rscript analyses/fit_pseudo_trajectories.R
Rscript analyses/fit_deep_downsampling.R                # needs the mBio supplement
Rscript analyses/fit_design_by_estimand.R
Rscript analyses/fit_power_calibration.R
Rscript analyses/fit_power_frontiers.R                  # ~4 min
Rscript analyses/fit_calendar_effects.R                 # ~6 min
Rscript analyses/fit_calendar_ou_null.R                 # ~6 min
Rscript AMR_figures_manuscript.R
Rscript analyses/figures_new.R
Rscript analyses/export_figures.R                       # figures/submission/
```

When a draft is renumbered, edit the table in `analyses/export_figures.R` and rerun it.

## Software

R 4.6.0 on macOS. Key packages: brms 2.23.0, lme4 2.0.1, rEDM 1.15.4, tidyverse 2.0.0 (dplyr 1.2.1, ggplot2 4.0.3), patchwork 1.3.2, readxl 1.4.5, ragg 1.5.2. PDF figures use the macOS Quartz device for font embedding.
