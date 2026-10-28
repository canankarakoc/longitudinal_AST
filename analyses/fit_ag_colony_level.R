if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
# AG at the COLONY denominator: prior resistance coded three ways.
#   1 any  = was any colony AG-resistant (0/1)
#   2 frac = fraction of colonies AG-resistant
#   3 ME   = that fraction with its Beta(1+k, 1+n-k) posterior SD
# Attenuation predicts the ME step should scale the coefficient by ~1/reliability.
suppressMessages({library(tidyverse); library(brms)}); source("AMR_common.R")

iso <- readRDS("data/pa_isolates.rds"); vis <- readRDS("data/pa_visit_full.rds")
ag <- c("amikacin","tobramycin","gentamicin")
sir <- c("Susceptible","Intermediate","Resistant")

# one row per visit: colonies tested for AG, colonies non-susceptible to any AG
cv <- iso %>%
  rowwise() %>%
  mutate(tested = any(c_across(all_of(ag)) %in% sir),
         resR   = any(c_across(all_of(ag)) %in% c("Resistant","Intermediate"))) %>%
  ungroup() %>% filter(tested) %>%
  group_by(individual_id, visit_number, date_collected_from_subject) %>%
  summarise(n_col = n(), k_col = sum(resR), .groups = "drop")

d <- cv %>% arrange(individual_id, date_collected_from_subject) %>%
  group_by(individual_id) %>%
  mutate(pk = lag(k_col), pn = lag(n_col)) %>% ungroup() %>%
  filter(!is.na(pk)) %>%
  mutate(prev_any  = as.integer(pk > 0),
         prev_frac = pk / pn,
         prev_pm   = (1 + pk) / (2 + pn),
         prev_psd  = sqrt((1+pk)*(1+pn-pk) / ((2+pn)^2 * (3+pn)))) %>%
  left_join(vis %>% select(individual_id, visit_number, prev_has_mrsa, prev_has_mssa,
                           prev_aminoglycoside_prescribed),
            by = c("individual_id","visit_number")) %>%
  filter(if_all(c(prev_has_mrsa, prev_has_mssa, prev_aminoglycoside_prescribed), ~!is.na(.)))

cat(sprintf("N = %d visits, %d patients | mean colonies %.2f | median prev_psd %.3f\n",
            nrow(d), n_distinct(d$individual_id), mean(d$n_col), median(d$prev_psd)))
cat(sprintf("prev_frac: %.0f%% at 0, %.0f%% at 1, %.0f%% in between\n",
    100*mean(d$prev_frac==0), 100*mean(d$prev_frac==1),
    100*mean(d$prev_frac>0 & d$prev_frac<1)))

cov <- "+ prev_has_mrsa + prev_has_mssa + prev_aminoglycoside_prescribed + (1 | individual_id)"
fitit <- function(pred, tag, pri) brm(
  as.formula(paste("k_col | trials(n_col) ~", pred, cov)),
  data = d, family = beta_binomial(), prior = pri,
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95), seed = 42,
  file = paste0("models/ag_col_", tag), file_refit = "on_change")

m1 <- fitit("prev_any",  "any",  priors_beta_random)
m2 <- fitit("prev_frac", "frac", priors_beta_random)
m3 <- fitit("me(prev_pm, prev_psd)", "me", c(
  prior(normal(0,1.5), class="b"), prior(normal(0,2), class="Intercept"),
  prior(exponential(1), class="sd"), prior(gamma(2,0.1), class="phi"),
  prior(normal(0.5,0.3), class="meanme"), prior(exponential(2), class="sdme")))

g <- function(m, v, lab) { s <- posterior_summary(m, variable=v)
  tibble(coding=lab, est=s[1,1], lo=s[1,3], hi=s[1,4], width=s[1,4]-s[1,3]) }
res <- bind_rows(g(m1,"b_prev_any","1 any-R (0/1)"),
                 g(m2,"b_prev_frac","2 fraction k/n"),
                 g(m3,"bsp_meprev_pmprev_psd","3 fraction + its SD"))
cat("\n=== AG, COLONY denominator: prior resistance three ways ===\n")
print(res %>% mutate(across(where(is.numeric), ~round(.,2))) %>% as.data.frame(), row.names=FALSE)
cat(sprintf("\n  ME / fraction ratio = %.2fx  (attenuation predicts ~1/reliability = 1.45x)\n",
            res$est[3]/res$est[2]))
write_csv(res, "results/heterogeneity/ag_colony_three_codings.csv")
