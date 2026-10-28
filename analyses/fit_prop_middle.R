# The missing middle model: prior resistance as a PROPORTION.
# naive = any-R (0/1) | prop = k/n | ME = latent theta with its own SD.
suppressMessages({library(tidyverse); library(brms)})
if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
source("AMR_common.R")

m_me <- readRDS("models/m2_me_comp.rds")
d <- m_me$data
v <- readRDS("data/pa_visit_full.rds") %>%
  transmute(individual_id, date_collected_from_subject,
            prev_prop = prev_prop_resistant_conservative)
# me_comp_df rows carry prev_n_resistant/prev_n_tested; rebuild the proportion
d$prev_prop <- d$prev_theta_pm    # placeholder, replaced below if available
if (all(c("prev_n_resistant","prev_n_tested") %in% names(d)))
  d$prev_prop <- d$prev_n_resistant / d$prev_n_tested

cat("prev_prop: median", round(median(d$prev_prop),3),
    "| range", paste(round(range(d$prev_prop),3), collapse="-"), "\n")

m_prop <- brm(
  n_resistant | trials(n_tested) ~ prev_prop +
    prev_has_mrsa + prev_has_mssa + prev_has_any_mucoid +
    prev_aminoglycoside_prescribed + prev_carbapenem_prescribed +
    prev_cephalosporin_prescribed + prev_fluoroquinolones_prescribed +
    prev_colistin_prescribed + (1 | individual_id),
  data = d, family = beta_binomial(), prior = priors_beta_random,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = "models/m2_prop_comp", file_refit = "on_change")
print(summary(m_prop))

get1 <- function(m, term, lab) {
  s <- posterior_summary(m, variable = term)
  tibble(model = lab, est = s[1,1], lo = s[1,3], hi = s[1,4])
}
cmp <- bind_rows(
  get1(readRDS("models/m2_naive_comp.rds"), "b_prev_any_resistance", "1. Naive: any-R (0/1)"),
  get1(m_prop, "b_prev_prop", "2. Proportion: k/n"),
  get1(m_me, "bsp_meprev_theta_pmprev_theta_psd", "3. Noise-aware: latent theta"))
cat("\n=== PRIOR-RESISTANCE EFFECT UNDER THREE CODINGS ===\n")
print(cmp %>% mutate(across(where(is.numeric), ~round(.,2)),
                     width = round(hi-lo,2)) %>% as.data.frame(), row.names=FALSE)
write_csv(cmp, "results/heterogeneity/prior_R_three_codings.csv")
