# Convergence diagnostics for every brms fit the manuscript reports (Table S1).
# Reads cached fits only. Run from the repo root:
#   Rscript analyses/convergence_table.R
if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
suppressPackageStartupMessages({library(brms); library(posterior); library(dplyr); library(readr)})

fits <- tribble(
  ~model,                                        ~file,
  "Table 1: aminoglycoside",                     "models/noise/vardecomp_ag.rds",
  "Table 1: carbapenem",                         "models/noise/vardecomp_carb.rds",
  "Table 1: fluoroquinolone",                    "models/noise/vardecomp_fq.rds",
  "Table 1: cephalosporin",                      "models/noise/vardecomp_ceph.rds",
  "Drug x isolate, aminoglycoside",              "models/sensitivity/m_var_isolate.rds",
  "Drug x isolate, visits with >= 2 colonies",   "models/sensitivity/m_var_isolate_multi.rds",
  "Ordinal robustness: cumulative logit",        "models/sir/m_ordinal_v2.rds",
  "Ordinal robustness: R + I vs S",              "models/sir/m_binary_conservative_v2.rds",
  "Ordinal robustness: R vs S + I",              "models/sir/m_binary_strict_v2.rds")

out <- purrr::pmap_dfr(fits, function(model, file) {
  if (!file.exists(file)) return(tibble(model = model, note = "fit not found"))
  m  <- readRDS(file)
  dr <- as_draws_df(m)
  keep <- grep("^(b_|bsp_|sd_|Intercept)", variables(dr), value = TRUE)
  s  <- summarise_draws(subset_draws(dr, variable = keep), "rhat", "ess_bulk", "ess_tail")
  vis <- grepl("visit", s$variable)
  np <- nuts_params(m)
  res <- tibble(model = model, observations = nobs(m),
                chains = m$fit@sim$chains, iter = m$fit@sim$iter, warmup = m$fit@sim$warmup,
                rhat_max = max(s$rhat), rhat_max_param = s$variable[which.max(s$rhat)],
                rhat_max_excl_visit_sd = max(s$rhat[!vis]),
                ess_bulk_min = min(s$ess_bulk), ess_tail_min = min(s$ess_tail),
                divergent = sum(np$Value[np$Parameter == "divergent__"]),
                draws = nrow(dr), note = "")
  rm(m, dr); gc(verbose = FALSE)
  res
})
print(as.data.frame(out), digits = 4)
dir.create("results/manuscript", showWarnings = FALSE)
write_csv(out, "results/manuscript/tableS1_convergence.csv")
