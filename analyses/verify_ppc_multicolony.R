# Decisive test of the PPC misfit: is sd_isolate inflated by single-colony
# visits (where isolate and visit effects are perfectly confounded), or is the
# over-prediction structural?
# Refit the ordinal model on MULTI-COLONY visits only, where both terms are
# identified from data with actual within-visit replication, then rerun the PPC.
suppressMessages({library(tidyverse); library(brms)})
if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
set.seed(42)
sir_levels <- c("Susceptible","Intermediate","Resistant")
sir_long <- readRDS("data/sir_long.rds")

multi <- sir_long %>% group_by(visit_key) %>%
  filter(n_distinct(isolate_id) >= 2) %>% ungroup()
cat("multi-colony frame:", nrow(multi), "rows |",
    n_distinct(multi$isolate_id), "colonies |",
    n_distinct(multi$visit_key), "visits |",
    n_distinct(multi$individual_id), "patients\n")

m_mc <- brm(sir ~ drug + (1|individual_id) + (1|visit_key) + (1|isolate_id),
  data = multi, family = cumulative("logit"),
  prior = c(prior(normal(0,4), class="Intercept"), prior(normal(0,2), class="b"),
            prior(exponential(1), class="sd")),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 12), seed = 42,
  file = "models/sir/m_ordinal_multicolony", file_refit = "on_change")
print(summary(m_mc))

mk <- function(d) d %>% inner_join(d, by=c("visit_key","drug"),
  relationship="many-to-many", suffix=c(".i",".j")) %>% filter(isolate_id.i != isolate_id.j)
obs <- mk(multi %>% select(visit_key,drug,isolate_id,sir)) %>%
  group_by(index=sir.i) %>% summarise(obs=mean(sir.i!=sir.j), .groups="drop")

yrep <- posterior_predict(m_mc, ndraws = 200)
key <- multi %>% select(visit_key, drug, isolate_id)
ppc <- map_dfr(seq_len(nrow(yrep)), function(i)
  mk(key %>% mutate(sir = factor(sir_levels[yrep[i,]], levels=sir_levels, ordered=TRUE))) %>%
    group_by(index=sir.i) %>% summarise(d=mean(sir.i!=sir.j), .groups="drop")) %>%
  group_by(index) %>% summarise(sim=median(d), lo=unname(quantile(d,.025)),
                                hi=unname(quantile(d,.975)), .groups="drop")

out <- obs %>% mutate(index=as.character(index)) %>%
  left_join(ppc %>% mutate(index=as.character(index)), by="index") %>%
  mutate(ratio = obs/sim)
cat("\n=== PPC ON MULTI-COLONY-ONLY REFIT ===\n")
print(out %>% mutate(across(where(is.numeric), ~round(.,3))) %>% as.data.frame(), row.names=FALSE)
cat("\n  Full-data model gave ratios 0.79 / 0.85 / 0.84.\n")
cat("  Ratios moving toward 1.0 here => the misfit was single-colony visits\n")
cat("  inflating sd_isolate (the 10.3d identifiability problem again).\n")
cat("  Ratios still ~0.8 => structural, and nothing to do with the isolate term.\n")

sd_cmp <- bind_rows(
  tibble(model="full data", term=c("patient","visit","colony"),
         sd=c(median(as_draws_df(readRDS("models/sir/m_ordinal_v2.rds"))[["sd_individual_id__Intercept"]]),
              median(as_draws_df(readRDS("models/sir/m_ordinal_v2.rds"))[["sd_visit_key__Intercept"]]),
              median(as_draws_df(readRDS("models/sir/m_ordinal_v2.rds"))[["sd_isolate_id__Intercept"]]))),
  tibble(model="multi-colony only", term=c("patient","visit","colony"),
         sd=c(median(as_draws_df(m_mc)[["sd_individual_id__Intercept"]]),
              median(as_draws_df(m_mc)[["sd_visit_key__Intercept"]]),
              median(as_draws_df(m_mc)[["sd_isolate_id__Intercept"]]))))
cat("\n=== RANDOM-EFFECT SDs ===\n")
print(sd_cmp %>% pivot_wider(names_from=model, values_from=sd) %>%
        mutate(across(where(is.numeric), ~round(.,3))) %>% as.data.frame(), row.names=FALSE)
write_csv(out, "results/sir/ppc_multicolony_refit.csv")
write_csv(sd_cmp, "results/sir/sd_full_vs_multicolony.csv")
