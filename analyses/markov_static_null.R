if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")
# Does the Markov transition matrix survive a null with NO dynamics?
# Null: each patient has a FIXED resistance frequency theta_i. Every visit draws
# colonies independently at that visit's observed depth. State = "any resistant".
# Any apparent transition is then pure sampling. If the null reproduces the
# observed matrix, the transition rates are measuring the lab, not the patient.
suppressMessages(library(tidyverse)); set.seed(42)
NSIM <- 500
v <- readRDS("data/pa_visit_full.rds")
cls <- c(ag="Aminoglycoside", carb="Carbapenem", fq="Fluoroquinolone", ceph="Cephalosporin")

tmat <- function(state, id) {
  tibble(id, state) %>% group_by(id) %>% mutate(prev = lag(state)) %>% ungroup() %>%
    filter(!is.na(prev)) %>%
    summarise(RR = sum(prev==1 & state==1)/max(sum(prev==1),1),
              SS = sum(prev==0 & state==0)/max(sum(prev==0),1),
              SR = sum(prev==0 & state==1)/max(sum(prev==0),1),
              RS = sum(prev==1 & state==0)/max(sum(prev==1),1))
}

res <- map_dfr(names(cls), function(c) {
  k <- v[[paste0(c,"_n_resistant")]]; n <- v[[paste0(c,"_n_tested")]]
  ok <- !is.na(k) & !is.na(n) & n > 0
  d <- tibble(id = v$individual_id[ok], k = k[ok], n = n[ok],
              vd = v$date_collected_from_subject[ok]) %>% arrange(id, vd)
  obs <- tmat(as.integer(d$k > 0), d$id)
  # fixed per-patient theta, shrunk slightly so patients at 0 or 1 are not absorbing
  th <- d %>% group_by(id) %>% summarise(theta = (sum(k)+0.5)/(sum(n)+1), .groups="drop")
  d2 <- d %>% left_join(th, by = "id")
  sim <- map_dfr(seq_len(NSIM), function(i)
    tmat(as.integer(rbinom(nrow(d2), d2$n, d2$theta) > 0), d2$id))
  tibble(class = unname(cls[c]),
         obs_RR = obs$RR, null_RR = mean(sim$RR),
         null_RR_lo = quantile(sim$RR,.025), null_RR_hi = quantile(sim$RR,.975),
         obs_SS = obs$SS, null_SS = mean(sim$SS),
         obs_SR = obs$SR, null_SR = mean(sim$SR),
         obs_RS = obs$RS, null_RS = mean(sim$RS))
})
cat("=== STATIC-COHORT NULL: no dynamics, only sampling ===\n")
cat("    (theta fixed per patient; colonies redrawn at the observed depths)\n\n")
print(res %>% transmute(class,
        `stay R obs` = round(obs_RR,2), `stay R null` = round(null_RR,2),
        `null 95%` = sprintf("[%.2f, %.2f]", null_RR_lo, null_RR_hi),
        `stay S obs` = round(obs_SS,2), `stay S null` = round(null_SS,2),
        `acquire obs` = round(obs_SR,2), `acquire null` = round(null_SR,2),
        `lose obs` = round(obs_RS,2), `lose null` = round(null_RS,2)) %>%
      as.data.frame(), row.names = FALSE)
write_csv(res, "results/noise/markov_static_null.csv")
