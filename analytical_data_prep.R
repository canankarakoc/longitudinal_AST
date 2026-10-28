# =============================================================================
# CF Emory Cohort Antibiotic Resistance Analysis 
# Author: Canan Karakoç
# Created: December 11, 2025
# Last Modified: June 9, 2026
# 
# Data prep for analysis and summary statistics
# =============================================================================

# =============================================================================
# PART 1: SETUP AND DATA LOADING
# =============================================================================

if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")

# Load packages
library(tidyverse)
library(readxl)
library(stringr)
library(magrittr)
library(janitor)
library(brms)
library(bayesplot)
library(patchwork)
library(lubridate)
library(markovchain)
library(rEDM)
library(lme4)

# Conflict resolution
#
# The conflicted package is deliberately NOT used here. It breaks brms mid-fit:
# rstan calls utils::apropos() while compiling a Stan model, apropos() calls
# exists() across every attached package, and conflicted aborts on the first
# ambiguous name it meets (ddirichlet, extraDistr vs brms). The failure shows up
# as "[conflicted] ddirichlet found in 2 packages" from inside brm(), nowhere
# near the library() calls.
#
# It matters here even though this script fits no models: if prep and
# AMR_pipeline.R are sourced in the SAME R session, conflicted stays attached
# and takes the pipeline's models down with it.
#
# Of the 14 preferences this block used to declare, only filter/select/lag are
# ever called in this file, and load order already resolves those to dplyr.
# Asserted rather than assumed:
.verbs <- c("filter", "select", "lag")
.owner <- vapply(.verbs, function(f)
  environmentName(environment(get(f))), character(1))
if (!all(.owner == "dplyr")) {
  stop("Namespace clash: ", paste(sprintf("%s -> %s", .verbs, .owner), collapse = ", "),
       "\n  Expected all three to resolve to dplyr. Qualify with dplyr:: or fix load order.")
}
rm(.verbs, .owner)

# Set theme for all plots
theme_set(theme_minimal(base_size = 14))

# Create output directories
for (dir in c("figures", "data", "models", "results")) {
  if (!dir.exists(dir)) dir.create(dir)
}

# Data paths
# The raw export belongs in data/ (never committed); the repository root is
# still checked for older checkouts.
clincal_microbe_data_csv <- "CDC_114_subjects_bacterial_isolates_04-28-2023.csv"
if (file.exists(file.path("data", clincal_microbe_data_csv)))
  clincal_microbe_data_csv <- file.path("data", clincal_microbe_data_csv)

# =============================================================================
# PART 2: DATA LOADING AND PREPROCESSING
# =============================================================================
# =============================================================================
# 2.1 CLEAN DATA LOADING - EXPLICIT HANDLING OF MISSING DATA
# 
# - 5 rows: Missing identification AND date → completely empty (DELETE)
# - 15 rows: Has identification, missing date → NO susceptibility data (DELETE)
# =============================================================================

clinical_microbe_data_raw <- read_csv(clincal_microbe_data_csv, col_types = cols(.default = "c")) %>%
  clean_names()

cat("Raw data loaded:", nrow(clinical_microbe_data_raw), "rows\n\n")

# -----------------------------------------------------------------------------
# Remove rows with missing identification (5 empty rows)
# -----------------------------------------------------------------------------
clinical_microbe_data <- clinical_microbe_data_raw %>%
  filter(
    !is.na(identification_of_isolate) & 
      identification_of_isolate != "" &
      identification_of_isolate != "NA"
  )

n_removed_no_id <- nrow(clinical_microbe_data_raw) - nrow(clinical_microbe_data)
cat("Removed", n_removed_no_id, "rows with missing identification (empty placeholder rows)\n")

# -----------------------------------------------------------------------------
# Remove rows with missing dates (15 rows - all have NO susceptibility data)
# -----------------------------------------------------------------------------

# First fix date typos. The corrections are record-level dates, so they live
# with the data in data/date_corrections.csv (columns raw, corrected), not here.
date_fix <- read_csv("data/date_corrections.csv", col_types = cols(.default = "c"))
clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    date_collected_from_subject = coalesce(
      date_fix$corrected[match(date_collected_from_subject, date_fix$raw)],
      date_collected_from_subject),
    date_collected_from_subject = as.Date(date_collected_from_subject)
  )

# Log what we're removing (for documentation)
missing_date_isolates <- clinical_microbe_data %>%
  filter(is.na(date_collected_from_subject)) %>%
  count(identification_of_isolate, name = "n_removed")

cat("\nRemoving isolates with missing dates (none have susceptibility data):\n")
print(missing_date_isolates)

# Remove missing dates
clinical_microbe_data <- clinical_microbe_data %>%
  filter(!is.na(date_collected_from_subject))

n_removed_no_date <- sum(missing_date_isolates$n_removed)
cat("\nRemoved", n_removed_no_date, "rows with missing dates\n")

# -----------------------------------------------------------------------------
# PRE-2017 FILTER  (flag - set to FALSE to keep the early visits)
#
# The cohort is 2017-2022 apart from a handful of 2012 and 2015 records whose
# dates are not yet confirmed. They are dropped by default until they are.
#
# This has to happen HERE, before section 2.5 builds the time variables:
# days_since_first_visit, visit_number and every prev_* lag are derived from
# each patient's first surviving visit, so filtering later in the analysis
# layer would leave those patients' clocks still anchored in 2012.
#
# Cost of dropping (measured on the 2026-08 data): 13 PA visits, 16 PA isolates
# (0.9%), 0 patients - all 10 affected patients continue after 2017, so only
# their earliest observation is lost. Benefit: removes 8 of the 9 longest
# inter-visit gaps (1609-1848 days, all 2012/2015 -> 2017 jumps), which
# otherwise dominate the Markov gap distribution and stretch the follow-up
# axis by 5-10 years for those patients. It does NOT solve the long-gap
# problem outright: 12 gaps over two years remain inside 2017-2022.
# -----------------------------------------------------------------------------

DROP_PRE_2017 <- TRUE   # <- flag

if (DROP_PRE_2017) {
  n_before  <- nrow(clinical_microbe_data)
  pat_before <- n_distinct(clinical_microbe_data$record_id)
  dropped <- clinical_microbe_data %>%
    filter(date_collected_from_subject < as.Date("2017-01-01"))
  clinical_microbe_data <- clinical_microbe_data %>%
    filter(date_collected_from_subject >= as.Date("2017-01-01"))
  cat("\nDROP_PRE_2017 = TRUE: removed", n_before - nrow(clinical_microbe_data),
      "isolate records dated before 2017-01-01\n")
  if (nrow(dropped) > 0) {
    cat("  years dropped:\n")
    print(dropped %>% count(year = lubridate::year(date_collected_from_subject)))
    cat("  patients affected:", n_distinct(dropped$record_id),
        "| patients lost entirely:",
        pat_before - n_distinct(clinical_microbe_data$record_id), "\n")
  }
} else {
  cat("\nDROP_PRE_2017 = FALSE: early (2012/2015) visits retained.\n")
}

cat("Total isolates:", nrow(clinical_microbe_data), "\n")
cat("Unique patients:", n_distinct(clinical_microbe_data$record_id), "\n\n")

# -----------------------------------------------------------------------------
# 2.2 Consolidate duplicate susceptibility columns
# -----------------------------------------------------------------------------

antibiotics_to_coalesce <- c("amikacin", "aztreonam", "cefepime", "ceftazidime",
                             "colistin", "gentamicin", "imipenem", "levofloxacin",
                             "meropenem", "piperacillin_tazobactam", "tobramycin",
                             "linezolid", "trimethoprim_sulfa", "minocycline",
                             "tigecycline", "ciprofloxacin", "clarithromycin",
                             "clofazamine", "moxifloxacin", "ceftriaxone")

clinical_microbe_data <- clinical_microbe_data %>%
  rename_with(~ "clinical_status", 
              .cols = starts_with("clinical_status_clinically_stable"))

for (ab in antibiotics_to_coalesce) {
  pattern <- paste0("^", ab, "(_\\d+)?$")
  if (any(str_detect(names(clinical_microbe_data), pattern))) {
    clinical_microbe_data <- clinical_microbe_data %>%
      mutate(!!ab := coalesce(!!! select(., matches(pattern))))
  }
}

# Remove duplicate columns
clinical_microbe_data <- clinical_microbe_data %>%
  select(-matches("^(amikacin|aztreonam|cefepime|ceftazidime|colistin|gentamicin|imipenem|levofloxacin|meropenem|piperacillin_tazobactam|tobramycin|linezolid|trimethoprim_sulfa|minocycline|tigecycline|ciprofloxacin|clarithromycin|clofazamine|moxifloxacin|ceftriaxone)_\\d+$")) %>%
  rename(linezolid_start_date = l_inezolid_start_date,
         linezolid_end_date = l_inezolid_end_date)

# -----------------------------------------------------------------------------
# 2.3 Collapse administration route columns
# -----------------------------------------------------------------------------

collapse_admin_routes_detailed <- function(df, drug_prefix) {
  iv_col <- paste0(drug_prefix, "_administration_route_choice_intravenous")
  oral_col <- paste0(drug_prefix, "_administration_route_choice_oral")
  inhaled_col <- paste0(drug_prefix, "_administration_route_choice_inhaled")
  
  cols_exist <- c(iv_col, oral_col, inhaled_col) %in% names(df)
  if (!any(cols_exist)) return(df)
  
  df %>%
    mutate(
      # Combined route string (existing)
      !!paste0(drug_prefix, "_route") := case_when(
        .data[[iv_col]] == "Checked" & .data[[oral_col]] == "Checked" & .data[[inhaled_col]] == "Checked" ~ "IV+Oral+Inhaled",
        .data[[iv_col]] == "Checked" & .data[[oral_col]] == "Checked" ~ "IV+Oral",
        .data[[iv_col]] == "Checked" & .data[[inhaled_col]] == "Checked" ~ "IV+Inhaled",
        .data[[oral_col]] == "Checked" & .data[[inhaled_col]] == "Checked" ~ "Oral+Inhaled",
        .data[[iv_col]] == "Checked" ~ "IV",
        .data[[oral_col]] == "Checked" ~ "Oral",
        .data[[inhaled_col]] == "Checked" ~ "Inhaled",
        TRUE ~ NA_character_
      ),
      
      # NEW: Binary indicators for each route (for modeling)
      !!paste0(drug_prefix, "_iv") := as.integer(.data[[iv_col]] == "Checked"),
      !!paste0(drug_prefix, "_oral") := as.integer(.data[[oral_col]] == "Checked"),
      !!paste0(drug_prefix, "_inhaled") := as.integer(.data[[inhaled_col]] == "Checked"),
      
      # NEW: Route intensity score (IV > Oral > Inhaled for systemic exposure)
      # This is a simplification - adjust weights based on PK knowledge
      !!paste0(drug_prefix, "_route_intensity") := case_when(
        .data[[iv_col]] == "Checked" ~ 3L,  # Highest systemic exposure
        .data[[oral_col]] == "Checked" ~ 2L,  # Moderate systemic exposure
        .data[[inhaled_col]] == "Checked" ~ 1L,  # Primarily local (lung) exposure
        TRUE ~ 0L
      )
    ) %>%
    select(-any_of(c(iv_col, oral_col, inhaled_col)))
}

drug_prefixes <- c("aminoglycoside", "antifungal", "carbapenem", "cephalosporin",
                   "fluoroquinolones", "lincosamide", "macrolide", "penicillin",
                   "tetracyclin", "aztreonam", "ceftazidime_avibactam", 
                   "chloramphenicol", "colistin", "ethambutol", "isoniazid",
                   "linezolid", "polymyxin_b", "rifabutin", "rifampin",
                   "trimethoprim_sulfamethoxazole_bactrim", "vancomycin")

for (prefix in drug_prefixes) {
  clinical_microbe_data <- collapse_admin_routes_detailed(clinical_microbe_data, prefix)
}

# -----------------------------------------------------------------------------
# 2.4 Drop unused/rare drug columns
# -----------------------------------------------------------------------------

drugs_to_drop <- c("antifungal", "chloramphenicol", "isoniazid", "polymyxin_b")
rare_drugs <- c("ethambutol", "macrolide", "rifabutin", "rifampin", 
                "lincosamide", "ceftazidime_avibactam")

clinical_microbe_data <- clinical_microbe_data %>%
  select(-matches(paste0("^(", paste(c(drugs_to_drop, rare_drugs), collapse = "|"), ")_")))

# -----------------------------------------------------------------------------
# 2.5 Create time variables
# -----------------------------------------------------------------------------

clinical_microbe_data <- clinical_microbe_data %>%
  group_by(record_id) %>%
  arrange(date_collected_from_subject) %>%
  mutate(
    days_since_first_visit = as.numeric(
      date_collected_from_subject - min(date_collected_from_subject, na.rm = TRUE)
    ),
    visit_number = dense_rank(date_collected_from_subject)
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# 2.6 Process co-infection data
# -----------------------------------------------------------------------------

clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    co_organisms = paste(
      if_else(microorganisms_choice_mssa == "Checked", "MSSA", ""),
      if_else(microorganisms_choice_mrsa == "Checked", "MRSA", ""),
      if_else(microorganisms_choice_pseudomonas_aeruginosa == "Checked", "Pseudomonas aeruginosa", ""),
      if_else(microorganisms_choice_stenotrophomonas_maltophilia == "Checked", "Stenotrophomonas maltophilia", ""),
      if_else(microorganisms_choice_aspergillus_fumigatus == "Checked", "Aspergillus fumigatus", ""),
      if_else(microorganisms_choice_mycobacterium_abscessus_mab == "Checked", "mycobacterium abscessus", ""),
      if_else(microorganisms_choice_mycobacterium_avium_complex_mac == "Checked", "Mycobacterium avium complex", ""),
      if_else(microorganisms_choice_achromobacter_xylosoxidans == "Checked", "Achromobacter xylosoxidans", ""),
      if_else(microorganisms_choice_burkholderia_cepacia_complex == "Checked", "Burkholderia cepacia complex", ""),
      if_else(microorganisms_choice_haemophilus_influenzae == "Checked", "Haemophilus influenzae", ""),
      if_else(microorganisms_choice_other == "Checked", other_organism, ""),
      sep = ", "),
    # -------------------------------------------------------------------------
    # FIX (5): clean co_organisms by SPLITTING rather than by regex substitution.
    # The old chain used str_replace_all(",\\s*,", ",") in a single pass, which
    # cannot collapse runs of 3+ commas ("A, , , B" -> "A, , B"), so unchecked
    # boxes left empty slots behind and n_coinfections counted them as organisms.
    # Splitting on comma is correct here: `other_organism` is free text that
    # itself holds comma-separated lists of multiple organisms.
    # The isolate's own species is removed by exact element match (the old
    # str_replace treated identification_of_isolate as a regex and matched
    # substrings).
    # -------------------------------------------------------------------------
    co_organisms = map2_chr(co_organisms, identification_of_isolate, function(s, self) {
      parts <- str_trim(str_split(s, ",")[[1]])
      parts <- parts[nzchar(parts) & parts != "NA"]
      parts <- unique(parts[parts != self])
      if (length(parts) == 0) NA_character_ else paste(parts, collapse = ", ")
    }),
    n_coinfections = if_else(is.na(co_organisms), 0L, str_count(co_organisms, ", ") + 1L)
  ) %>%
  select(-starts_with("microorganisms_choice_"), -other_organism)

# -----------------------------------------------------------------------------
# 2.7 Clean up isolate identification and phenotype
# -----------------------------------------------------------------------------

clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    identification_of_isolate = if_else(
      identification_of_isolate == "Other" & !is.na(if_other_please_type_the_name_of_the_isolate),
      if_other_please_type_the_name_of_the_isolate,
      identification_of_isolate
    )
  ) %>%
  mutate(
    # Original phenotype categories
    pa_phenotype = case_when(
      pseudomonas_aeruginosa_co_culture_phenotype_choice_mucoid == "Checked" & 
        pseudomonas_aeruginosa_co_culture_phenotype_choice_non_mucoid == "Checked" ~ "Mixed",
      pseudomonas_aeruginosa_co_culture_phenotype_choice_mucoid == "Checked" ~ "Mucoid",
      pseudomonas_aeruginosa_co_culture_phenotype_choice_non_mucoid == "Checked" ~ "Non-mucoid",
      TRUE ~ NA_character_
    ),
    
    # Binary variables for modeling
    has_any_mucoid = as.integer(pa_phenotype %in% c("Mucoid", "Mixed")),
    has_pure_mucoid = as.integer(pa_phenotype == "Mucoid"),
    has_mixed_phenotype = as.integer(pa_phenotype == "Mixed"),
    has_nonmucoid_only = as.integer(pa_phenotype == "Non-mucoid")
  ) %>%
  select(
    -if_other_please_type_the_name_of_the_isolate,  # ADD: remove after using
    -starts_with("pseudomonas_aeruginosa_co_culture_phenotype_choice_")
  )


# -----------------------------------------------------------------------------
# 2.8 Create anonymous patient ID and drop irrelevant columns
# -----------------------------------------------------------------------------

clinical_microbe_data <- clinical_microbe_data %>%
  mutate(individual_id = as.numeric(factor(record_id)))

# Persist the CFBR number <-> individual_id mapping BEFORE record_id is dropped.
# record_id holds the CFBR number (format "CFBR-EM-nnnn"), which is how other labs
# refer to the same source subject. Once the original csv is deleted this mapping
# is unrecoverable, so write it out every run. See make_cfbr_crosswalk.R for the
# standalone regenerator.
cfbr_id_crosswalk <- clinical_microbe_data %>%
  group_by(individual_id, cfbr_number = record_id) %>%
  summarise(n_isolates  = n(),
            n_visits    = n_distinct(date_collected_from_subject),
            first_visit = min(date_collected_from_subject),
            last_visit  = max(date_collected_from_subject),
            .groups = "drop") %>%
  arrange(individual_id)

stopifnot(!any(duplicated(cfbr_id_crosswalk$individual_id)),
          !any(duplicated(cfbr_id_crosswalk$cfbr_number)))

write_csv(cfbr_id_crosswalk, "data/cfbr_id_crosswalk.csv")
saveRDS(cfbr_id_crosswalk, "data/cfbr_id_crosswalk.rds")
cat("Saved CFBR crosswalk for", nrow(cfbr_id_crosswalk), "patients\n")

clinical_microbe_data <- clinical_microbe_data %>%
  select(-record_id) %>%
  select(individual_id, everything()) %>%
  select(-contains(c("lab_distributing_isolate", "investigator_lab", "variant_id", 
                     "freezer_location", "aliquot_id", "box", "well_location", "complete")))


# Drop columns with >90% NA
# never drop PA-relevant susceptibility columns even if sparsely
# tested. Colistin / aztreonam are often >90% untested but must stay in the
# overall resistance denominator. Protect them (and clinical_status) explicitly.
protected_from_na_drop <- c(
  "amikacin", "aztreonam", "cefepime", "ceftazidime", "colistin",
  "gentamicin", "imipenem", "levofloxacin", "meropenem",
  "piperacillin_tazobactam", "tobramycin", "ciprofloxacin",
  "clinical_status"
)
na_counts <- sapply(clinical_microbe_data, function(x) sum(is.na(x)))
cols_to_drop <- setdiff(
  names(na_counts)[na_counts / nrow(clinical_microbe_data) > 0.90],
  protected_from_na_drop
)
clinical_microbe_data <- clinical_microbe_data %>% select(-all_of(cols_to_drop))

# -----------------------------------------------------------------------------
# 2.9 Create resistance outcome variables
# -----------------------------------------------------------------------------

# All susceptibility columns in the dataset (for all organisms)
susceptibility_cols <- c(
  "amikacin", "aztreonam", "cefepime", "ceftazidime", "colistin",
  "gentamicin", "imipenem", "levofloxacin", "meropenem",
  "piperacillin_tazobactam", "tobramycin", "ciprofloxacin",
  "trimethoprim_sulfa", "linezolid", "tetracycline", "vancomycin",
  "ticarcillin_clavulanate", "minocycline", "tigecycline",
  "amoxicillin_clavulanate", "azithromycin", "cefoxitin",
  "clarithromycin", "clindamycin", "erythromycin", "doxycycline", 
  "moxifloxacin", "oxacillin", "penicillin", "daptomycin"
)

# PA-specific antibiotics (clinically relevant for Pseudomonas)
susceptibility_cols_pa <- c(
  "amikacin", "aztreonam", "cefepime", "ceftazidime", "colistin",
  "gentamicin", "imipenem", "levofloxacin", "meropenem",
  "piperacillin_tazobactam", "tobramycin", "ciprofloxacin"
)

# Class-specific columns
carb_cols <- c("imipenem", "meropenem")
ag_cols <- c("amikacin", "tobramycin", "gentamicin")
fq_cols <- c("ciprofloxacin", "levofloxacin")
ceph_cols <- c("cefepime", "ceftazidime")

# Check which columns exist
existing_pa_cols <- susceptibility_cols_pa[susceptibility_cols_pa %in% names(clinical_microbe_data)]

# every susceptibility column that survived the NA drop, not just the
# 12 PA ones. Previously the non-PA columns kept their raw values, so
# "Not Tested" survived as a literal level and the "Susceptibility" typo
# (present only in the numbered duplicate columns, hence still reachable after
# coalescing) was never corrected in them.
existing_susc_cols <- susceptibility_cols[susceptibility_cols %in% names(clinical_microbe_data)]


# correct clinical status 

# If "Undetermined" and NA are the same, consolidate them
clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    clinical_status = case_when(
      clinical_status == "Indetermined" ~ NA_character_,
      clinical_status == "" ~ NA_character_,  # Also catch empty strings
      TRUE ~ clinical_status
    ),
    # Make it a proper factor with meaningful levels
    clinical_status = factor(clinical_status, levels = c(
      "Clinically Stable",
      "APE Outpatient",
      "APE IV Antibiotics"
    ))
  )

# =============================================================================
# FIX: CLEAN SUSCEPTIBILITY VALUES
# =============================================================================
# Apply fix: Standardize susceptibility values
clinical_microbe_data <- clinical_microbe_data %>%
  mutate(across(
    all_of(existing_susc_cols),
    ~case_when(
      . == "Susceptibility" ~ "Susceptible",  # FIX TYPO
      . == "Susceptible" ~ "Susceptible",
      . == "Resistant" ~ "Resistant",
      . == "Intermediate" ~ "Intermediate",
      . %in% c("Not Tested", "Unchecked", "NA", "") ~ NA_character_,
      TRUE ~ NA_character_  # Anything else becomes NA
    )
  ))

# the previous rowwise pass that computed n_resistant / n_tested
# here was redundant - the RECALCULATE block below supersedes it. Its only
# unique outputs (prop_resistance, n_resistant_pa, n_tested_pa) are now created
# inside RECALCULATE so nothing downstream loses a column.

# =============================================================================
# FIX: RECALCULATE RESISTANCE COUNTS
# =============================================================================

# Helper function to count values across columns
count_values <- function(df, cols, values) {
  existing <- cols[cols %in% names(df)]
  if (length(existing) == 0) return(rep(0L, nrow(df)))
  
  df %>%
    select(all_of(existing)) %>%
    mutate(across(everything(), ~. %in% values)) %>%
    rowSums(na.rm = TRUE)
}

# Calculate resistance metrics
clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    # Overall PA resistance
    n_resistant = count_values(., existing_pa_cols, c("Resistant", "Intermediate")),
    n_susceptible = count_values(., existing_pa_cols, "Susceptible"),
    n_tested = n_resistant + n_susceptible,
    
    # Three coding approaches
    n_resistant_strict = count_values(., existing_pa_cols, "Resistant"),  # I = not resistant
    n_resistant_conservative = count_values(., existing_pa_cols, c("Resistant", "Intermediate")),  # I = resistant
    n_intermediate = count_values(., existing_pa_cols, "Intermediate"),
    
    # Proportions
    prop_resistant_strict = if_else(n_tested > 0, n_resistant_strict / n_tested, NA_real_),
    prop_resistant_conservative = if_else(n_tested > 0, n_resistant_conservative / n_tested, NA_real_),
    prop_resistant_weighted = if_else(n_tested > 0, 
                                      (n_resistant_strict + 0.5 * n_intermediate) / n_tested, 
                                      NA_real_),
    
    # Binary indicators
    # NA, not 0, when nothing was tested - "not tested" is not
    # "not resistant".: state the coding convention in the expression
    # so it cannot drift when n_resistant is redefined further down.
    any_resistance = if_else(n_tested > 0,
                             as.integer(n_resistant_conservative > 0),
                             NA_integer_)

    # prop_resistance / n_resistant_pa / n_tested_pa removed. They were
    # duplicates of prop_resistant_conservative / n_resistant / n_tested kept
    # "for backward compatibility", but pa_isolates redefines n_resistant to the
    # STRICT count below without refreshing them, so they went stale (958 rows
    # disagreed with prop_resistant_strict). Nothing downstream reads them.

    # multi_drug_resistant intentionally NOT computed here: the old
    # `n_resistant >= 3` counted individual DRUGS, not antibiotic CLASSES.
    # Class-based MDR is computed after the class-specific block (point 2).
  )

# =============================================================================
# CLASS-SPECIFIC RESISTANCE
# =============================================================================

# Carbapenem
carb_existing <- carb_cols[carb_cols %in% names(clinical_microbe_data)]
clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    carb_n_resistant = count_values(., carb_existing, c("Resistant", "Intermediate")),
    carb_n_susceptible = count_values(., carb_existing, "Susceptible"),
    carb_n_tested = carb_n_resistant + carb_n_susceptible,
    carb_prop_resistant = if_else(carb_n_tested > 0, carb_n_resistant / carb_n_tested, NA_real_),
    carb_any_resistant_conservative = if_else(carb_n_tested > 0, as.integer(carb_n_resistant > 0), NA_integer_)  # FIX (4)
  )

# Aminoglycoside
ag_existing <- ag_cols[ag_cols %in% names(clinical_microbe_data)]
clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    ag_n_resistant = count_values(., ag_existing, c("Resistant", "Intermediate")),
    ag_n_susceptible = count_values(., ag_existing, "Susceptible"),
    ag_n_tested = ag_n_resistant + ag_n_susceptible,
    ag_prop_resistant = if_else(ag_n_tested > 0, ag_n_resistant / ag_n_tested, NA_real_),
    ag_any_resistant_conservative = if_else(ag_n_tested > 0, as.integer(ag_n_resistant > 0), NA_integer_)  # FIX (4)
  )

# Fluoroquinolone
fq_existing <- fq_cols[fq_cols %in% names(clinical_microbe_data)]
clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    fq_n_resistant = count_values(., fq_existing, c("Resistant", "Intermediate")),
    fq_n_susceptible = count_values(., fq_existing, "Susceptible"),
    fq_n_tested = fq_n_resistant + fq_n_susceptible,
    fq_prop_resistant = if_else(fq_n_tested > 0, fq_n_resistant / fq_n_tested, NA_real_),
    fq_any_resistant_conservative = if_else(fq_n_tested > 0, as.integer(fq_n_resistant > 0), NA_integer_)  # FIX (4)
  )

# Cephalosporin
ceph_existing <- ceph_cols[ceph_cols %in% names(clinical_microbe_data)]
clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    ceph_n_resistant = count_values(., ceph_existing, c("Resistant", "Intermediate")),
    ceph_n_susceptible = count_values(., ceph_existing, "Susceptible"),
    ceph_n_tested = ceph_n_resistant + ceph_n_susceptible,
    ceph_prop_resistant = if_else(ceph_n_tested > 0, ceph_n_resistant / ceph_n_tested, NA_real_),
    ceph_any_resistant_conservative = if_else(ceph_n_tested > 0, as.integer(ceph_n_resistant > 0), NA_integer_)  # FIX (4)
  )

# =============================================================================
# CLASS-BASED MULTI-DRUG RESISTANCE (isolate level)
# MDR-PA = non-susceptible (R or I) to >=1 agent in >=3 antipseudomonal CLASSES
# (Magiorakos et al. 2012), NOT >=3 individual drugs. Classes considered:
# aminoglycoside, carbapenem, cephalosporin, fluoroquinolone, penicillin
# (pip/tazo), monobactam (aztreonam), polymyxin (colistin).
# MDR is left NA when fewer than 3 classes were tested (cannot be determined).
# =============================================================================
nonsusc_iso    <- function(x) as.integer(x %in% c("Resistant", "Intermediate"))
tested_any_iso <- function(x) as.integer(x %in% c("Resistant", "Intermediate", "Susceptible"))

clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    pen_nonsusc_iso  = nonsusc_iso(piperacillin_tazobactam),
    pen_tested_iso   = tested_any_iso(piperacillin_tazobactam),
    mono_nonsusc_iso = nonsusc_iso(aztreonam),
    mono_tested_iso  = tested_any_iso(aztreonam),
    poly_nonsusc_iso = nonsusc_iso(colistin),
    poly_tested_iso  = tested_any_iso(colistin),
    
    n_classes_tested = (ag_n_tested > 0) + (carb_n_tested > 0) +
      (ceph_n_tested > 0) + (fq_n_tested > 0) +
      pen_tested_iso + mono_tested_iso + poly_tested_iso,
    # FIX (4): the *_any_resistant_conservative flags are now NA for untested
    # classes, so coalesce to 0 here - an untested class contributes no
    # non-susceptible class to the MDR count. Whether MDR is determinable at
    # all is handled by n_classes_tested below.
    n_classes_resistant = coalesce(ag_any_resistant_conservative, 0L) +
      coalesce(carb_any_resistant_conservative, 0L) +
      coalesce(ceph_any_resistant_conservative, 0L) +
      coalesce(fq_any_resistant_conservative, 0L) +
      pen_nonsusc_iso + mono_nonsusc_iso + poly_nonsusc_iso,
    
    multi_drug_resistant = case_when(
      n_classes_tested < 3     ~ NA_integer_,
      n_classes_resistant >= 3 ~ 1L,
      TRUE                     ~ 0L
    )
  ) %>%
  select(-pen_nonsusc_iso, -pen_tested_iso, -mono_nonsusc_iso,
         -mono_tested_iso, -poly_nonsusc_iso, -poly_tested_iso)

# -----------------------------------------------------------------------------
# 2.10 Create binary prescription variables
# -----------------------------------------------------------------------------

prescription_cols <- names(clinical_microbe_data)[str_detect(names(clinical_microbe_data), "_prescribed$")]

clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    across(all_of(prescription_cols), ~if_else(. == "Yes", 1L, 0L)),
    clinical_status = factor(clinical_status)
  ) %>%
  rowwise() %>%
  mutate(n_drug_classes = sum(c_across(all_of(prescription_cols)), na.rm = TRUE)) %>%
  ungroup()

cat("Isolate-level data:", nrow(clinical_microbe_data), "rows from", 
    n_distinct(clinical_microbe_data$individual_id), "patients\n")

# =============================================================================
# PART 3: CREATE FINAL DATASETS
# =============================================================================

# Function to calculate duration from start/end dates
# For pilot dataset it is not useful but for a large dataset it is a useful information
calculate_treatment_duration <- function(df, drug_prefix) {
  start_col <- paste0(drug_prefix, "_start_date")
  end_col <- paste0(drug_prefix, "_end_date")
  
  if (!start_col %in% names(df)) {
    message(paste0("Note: ", start_col, " not found, skipping duration for ", drug_prefix))
    return(df)
  }
  
  if (!end_col %in% names(df)) {
    message(paste0("Note: ", end_col, " not found, skipping duration for ", drug_prefix))
    return(df)
  }
  
  # Handle different naming conventions
  start_col_alt <- str_replace(start_col, "_", " ") %>% str_to_title() %>% str_replace_all(" ", "_")
  end_col_alt <- str_replace(end_col, "_", " ") %>% str_to_title() %>% str_replace_all(" ", "_")
  
  if (!any(c(start_col, start_col_alt) %in% names(df))) return(df)
  
  # Find the actual column names
  start_actual <- intersect(c(start_col, start_col_alt), names(df))[1]
  end_actual <- intersect(c(end_col, end_col_alt), names(df))[1]
  
  if (is.na(start_actual) || is.na(end_actual)) return(df)
  
  df %>%
    mutate(
      # Parse dates
      .start_date = as.Date(.data[[start_actual]]),
      .end_date = as.Date(.data[[end_actual]]),
      
      # Calculate duration in days
      !!paste0(drug_prefix, "_duration_days") := as.numeric(
        difftime(.end_date, .start_date, units = "days")
      ),
      
      # Categorize duration
      !!paste0(drug_prefix, "_duration_cat") := case_when(
        is.na(.data[[paste0(drug_prefix, "_duration_days")]]) ~ NA_character_,
        .data[[paste0(drug_prefix, "_duration_days")]] <= 7 ~ "short",
        .data[[paste0(drug_prefix, "_duration_days")]] <= 14 ~ "medium",
        .data[[paste0(drug_prefix, "_duration_days")]] <= 28 ~ "long",
        TRUE ~ "very_long"
      )
    ) %>%
    select(-.start_date, -.end_date)
}

# Apply to all drug classes
drug_prefixes_for_duration <- c(
  "aminoglycoside", "carbapenem", "cephalosporin", "fluoroquinolone",
  "penicillin", "aztreonam", "colistin", "tetracyclin", "vancomycin",
  "trimethoprim_sulfamethoxazole_bactrim"
)

for (prefix in drug_prefixes_for_duration) {
  clinical_microbe_data <- calculate_treatment_duration(clinical_microbe_data, prefix)
}

#Currently not enough exposure data -REVISIT durations!

clinical_microbe_data <- clinical_microbe_data %>%
  mutate(
    # -------------------------------------------------------------------------
    # AMINOGLYCOSIDE EXPOSURE (route only - duration data too sparse)
    # -------------------------------------------------------------------------
    aminoglycoside_exposure_score = case_when(
      is.na(aminoglycoside_prescribed) ~ NA_integer_,  # FIX (3): no rx data != no exposure
      aminoglycoside_prescribed == 0 ~ 0L,
      aminoglycoside_iv == 1 ~ 3L,      # IV: highest systemic
      aminoglycoside_inhaled == 1 ~ 2L, # Inhaled: high local, low systemic
      aminoglycoside_oral == 1 ~ 1L,    # Oral: rare for aminoglycosides
      TRUE ~ 1L                          # Prescribed, route unknown
    ),
    
    # -------------------------------------------------------------------------
    # CARBAPENEM EXPOSURE
    # -------------------------------------------------------------------------
    carbapenem_exposure_score = case_when(
      is.na(carbapenem_prescribed) ~ NA_integer_,  # FIX (3): no rx data != no exposure
      carbapenem_prescribed == 0 ~ 0L,
      carbapenem_iv == 1 ~ 3L,
      carbapenem_oral == 1 ~ 2L,  # Only ertapenem is oral
      TRUE ~ 1L
    ),
    
    # -------------------------------------------------------------------------
    # CEPHALOSPORIN EXPOSURE
    # -------------------------------------------------------------------------
    cephalosporin_exposure_score = case_when(
      is.na(cephalosporin_prescribed) ~ NA_integer_,  # FIX (3): no rx data != no exposure
      cephalosporin_prescribed == 0 ~ 0L,
      cephalosporin_iv == 1 ~ 3L,
      cephalosporin_oral == 1 ~ 2L,
      TRUE ~ 1L
    ),
    
    # -------------------------------------------------------------------------
    # FLUOROQUINOLONE EXPOSURE
    # Note: FQs have excellent oral bioavailability, so oral ≈ IV
    # -------------------------------------------------------------------------
    fluoroquinolone_exposure_score = case_when(
      is.na(fluoroquinolones_prescribed) ~ NA_integer_,  # FIX (3): no rx data != no exposure
      fluoroquinolones_prescribed == 0 ~ 0L,
      fluoroquinolones_iv == 1 ~ 3L,
      fluoroquinolones_oral == 1 ~ 3L,    # Oral FQ has ~90% bioavailability
      fluoroquinolones_inhaled == 1 ~ 2L,
      TRUE ~ 1L
    ),
    
    # -------------------------------------------------------------------------
    # PENICILLIN COMBO EXPOSURE (pip/tazo)
    # -------------------------------------------------------------------------
    penicillin_exposure_score = case_when(
      is.na(penicillin_prescribed) ~ NA_integer_,  # FIX (3): no rx data != no exposure
      penicillin_prescribed == 0 ~ 0L,
      penicillin_iv == 1 ~ 3L,
      penicillin_oral == 1 ~ 2L,
      TRUE ~ 1L
    ),
    
    # -------------------------------------------------------------------------
    # AZTREONAM EXPOSURE (important in CF - Cayston is inhaled)
    # -------------------------------------------------------------------------
    aztreonam_exposure_score = case_when(
      is.na(aztreonam_prescribed) ~ NA_integer_,  # FIX (3): no rx data != no exposure
      aztreonam_prescribed == 0 ~ 0L,
      aztreonam_iv == 1 ~ 3L,
      aztreonam_inhaled == 1 ~ 2L,  # Cayston
      TRUE ~ 1L
    ),
    
    # -------------------------------------------------------------------------
    # COLISTIN EXPOSURE (last-resort)
    # -------------------------------------------------------------------------
    colistin_exposure_score = case_when(
      is.na(colistin_prescribed) ~ NA_integer_,  # FIX (3): no rx data != no exposure
      colistin_prescribed == 0 ~ 0L,
      colistin_iv == 1 ~ 3L,       # IV colistin is nephrotoxic
      colistin_inhaled == 1 ~ 2L,  # Inhaled more common for maintenance
      TRUE ~ 1L
    ),
    
    # -------------------------------------------------------------------------
    # TOTAL ANTIBIOTIC PRESSURE
    # FIX (3): plain `+`, so the total is NA whenever any component class is
    # unknown - an unknown class means the total really is unknown. Isolates
    # with no prescription data at all (the majority before 2017) now get NA
    # here instead of 7, which is what the old `TRUE ~ 1L` fall-through gave
    # them: one point per class for data that was never collected.
    # -------------------------------------------------------------------------
    total_antibiotic_pressure =
      aminoglycoside_exposure_score + 
      carbapenem_exposure_score + 
      cephalosporin_exposure_score + 
      fluoroquinolone_exposure_score + 
      penicillin_exposure_score +
      aztreonam_exposure_score +
      colistin_exposure_score
  )

# -----------------------------------------------------------------------------
# ISOLATE-LEVEL (only isolate-level calculations)
# -----------------------------------------------------------------------------

pa_isolates <- clinical_microbe_data %>%
  filter(identification_of_isolate == "Pseudomonas aeruginosa") %>%
  mutate(
    # =========================================================================
    # CLASS-SPECIFIC RESISTANCE AT ISOLATE LEVEL
    #
    # NAMING CONVENTION, READ BEFORE EDITING:
    # The *_n_resistant columns below are STRICT (Resistant only). They shadow
    # the CONSERVATIVE (R + I) columns of the same name built in the
    # clinical_microbe_data block above. That is deliberate - the visit-level
    # summarise sums these and then adds *_n_intermediate to reach the
    # conservative count - but it means "n_resistant" means different things in
    # clinical_microbe_data_cleaned.* and in pa_isolates.*.
    # Explicit *_n_resistant_strict aliases are written at the end of this
    # mutate so downstream code never has to infer the convention.
    # =========================================================================

    # Aminoglycosides
    ag_n_susceptible = rowSums(across(c(amikacin, tobramycin, gentamicin), 
                                      ~. == "Susceptible"), na.rm = TRUE),
    ag_n_intermediate = rowSums(across(c(amikacin, tobramycin, gentamicin), 
                                       ~. == "Intermediate"), na.rm = TRUE),
    ag_n_resistant = rowSums(across(c(amikacin, tobramycin, gentamicin), 
                                    ~. == "Resistant"), na.rm = TRUE),
    ag_n_tested = ag_n_susceptible + ag_n_intermediate + ag_n_resistant,
    
    # Carbapenems
    carb_n_susceptible = rowSums(across(c(meropenem, imipenem), 
                                        ~. == "Susceptible"), na.rm = TRUE),
    carb_n_intermediate = rowSums(across(c(meropenem, imipenem), 
                                         ~. == "Intermediate"), na.rm = TRUE),
    carb_n_resistant = rowSums(across(c(meropenem, imipenem), 
                                      ~. == "Resistant"), na.rm = TRUE),
    carb_n_tested = carb_n_susceptible + carb_n_intermediate + carb_n_resistant,
    
    # Cephalosporins
    ceph_n_susceptible = rowSums(across(c(cefepime, ceftazidime), 
                                        ~. == "Susceptible"), na.rm = TRUE),
    ceph_n_intermediate = rowSums(across(c(cefepime, ceftazidime), 
                                         ~. == "Intermediate"), na.rm = TRUE),
    ceph_n_resistant = rowSums(across(c(cefepime, ceftazidime), 
                                      ~. == "Resistant"), na.rm = TRUE),
    ceph_n_tested = ceph_n_susceptible + ceph_n_intermediate + ceph_n_resistant,
    
    # Fluoroquinolones
    fq_n_susceptible = rowSums(across(c(ciprofloxacin, levofloxacin), 
                                      ~. == "Susceptible"), na.rm = TRUE),
    fq_n_intermediate = rowSums(across(c(ciprofloxacin, levofloxacin), 
                                       ~. == "Intermediate"), na.rm = TRUE),
    fq_n_resistant = rowSums(across(c(ciprofloxacin, levofloxacin), 
                                    ~. == "Resistant"), na.rm = TRUE),
    fq_n_tested = fq_n_susceptible + fq_n_intermediate + fq_n_resistant,
    
    # Single-drug classes (binary at isolate)
    pen_susceptible = as.integer(piperacillin_tazobactam == "Susceptible"),
    pen_intermediate = as.integer(piperacillin_tazobactam == "Intermediate"),
    pen_resistant = as.integer(piperacillin_tazobactam == "Resistant"),
    pen_tested = as.integer(piperacillin_tazobactam %in% c("Susceptible", "Intermediate", "Resistant")),
    
    mono_susceptible = as.integer(aztreonam == "Susceptible"),
    mono_intermediate = as.integer(aztreonam == "Intermediate"),
    mono_resistant = as.integer(aztreonam == "Resistant"),
    mono_tested = as.integer(aztreonam %in% c("Susceptible", "Intermediate", "Resistant")),
    
    poly_susceptible = as.integer(colistin == "Susceptible"),
    poly_intermediate = as.integer(colistin == "Intermediate"),
    poly_resistant = as.integer(colistin == "Resistant"),
    poly_tested = as.integer(colistin %in% c("Susceptible", "Intermediate", "Resistant")),
    
    # Overall resistance counts (isolate level)
    n_susceptible = rowSums(across(all_of(susceptibility_cols_pa), 
                                   ~. == "Susceptible"), na.rm = TRUE),
    n_intermediate = rowSums(across(all_of(susceptibility_cols_pa), 
                                    ~. == "Intermediate"), na.rm = TRUE),
    n_resistant = rowSums(across(all_of(susceptibility_cols_pa), 
                                 ~. == "Resistant"), na.rm = TRUE),
    n_tested = rowSums(across(all_of(susceptibility_cols_pa), 
                              ~. %in% c("Susceptible", "Intermediate", "Resistant")), 
                       na.rm = TRUE),
    
    # Proportions at isolate level
    prop_resistant_conservative = if_else(n_tested > 0, 
                                          (n_resistant + n_intermediate) / n_tested, 
                                          NA_real_),
    prop_resistant_strict = if_else(n_tested > 0, 
                                    n_resistant / n_tested, 
                                    NA_real_),
    prop_resistant_weighted = if_else(n_tested > 0,
                                      (n_resistant + 0.5 * n_intermediate) / n_tested,
                                      NA_real_),

    # -------------------------------------------------------------------------
    # FIX (7): explicit, self-describing aliases for the counts recomputed above,
    # plus refreshed conservative counts so nothing carried over from
    # clinical_microbe_data can go stale against the redefined n_resistant.
    # -------------------------------------------------------------------------
    n_resistant_strict       = n_resistant,
    ag_n_resistant_strict    = ag_n_resistant,
    carb_n_resistant_strict  = carb_n_resistant,
    ceph_n_resistant_strict  = ceph_n_resistant,
    fq_n_resistant_strict    = fq_n_resistant,

    n_resistant_conservative = n_resistant + n_intermediate,
    ag_n_resistant_cons      = ag_n_resistant   + ag_n_intermediate,
    carb_n_resistant_cons    = carb_n_resistant + carb_n_intermediate,
    ceph_n_resistant_cons    = ceph_n_resistant + ceph_n_intermediate,
    fq_n_resistant_cons      = fq_n_resistant   + fq_n_intermediate,

    # These four were computed from the CONSERVATIVE counts up in
    # clinical_microbe_data and were left untouched here; recompute them from
    # the isolate-level numbers so they agree with everything around them.
    ag_prop_resistant   = if_else(ag_n_tested   > 0, ag_n_resistant_cons   / ag_n_tested,   NA_real_),
    carb_prop_resistant = if_else(carb_n_tested > 0, carb_n_resistant_cons / carb_n_tested, NA_real_),
    ceph_prop_resistant = if_else(ceph_n_tested > 0, ceph_n_resistant_cons / ceph_n_tested, NA_real_),
    fq_prop_resistant   = if_else(fq_n_tested   > 0, fq_n_resistant_cons   / fq_n_tested,   NA_real_),

    ag_any_resistant_conservative   = if_else(ag_n_tested   > 0, as.integer(ag_n_resistant_cons   > 0), NA_integer_),
    carb_any_resistant_conservative = if_else(carb_n_tested > 0, as.integer(carb_n_resistant_cons > 0), NA_integer_),
    ceph_any_resistant_conservative = if_else(ceph_n_tested > 0, as.integer(ceph_n_resistant_cons > 0), NA_integer_),
    fq_any_resistant_conservative   = if_else(fq_n_tested   > 0, as.integer(fq_n_resistant_cons   > 0), NA_integer_),

    any_resistance = if_else(n_tested > 0, as.integer(n_resistant_conservative > 0), NA_integer_)
  )

# -----------------------------------------------------------------------------
# VISIT-LEVEL AGGREGATION (all aggregation in one place)
# -----------------------------------------------------------------------------

safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

pa_visit_full <- pa_isolates %>%
  group_by(individual_id, date_collected_from_subject) %>%
  summarise(
    # =========================================================================
    # CLASS-SPECIFIC RESISTANCE (aggregated)
    # =========================================================================
    
    # Aminoglycosides
    ag_n_susceptible = sum(ag_n_susceptible, na.rm = TRUE),
    ag_n_intermediate = sum(ag_n_intermediate, na.rm = TRUE),
    ag_n_resistant = sum(ag_n_resistant, na.rm = TRUE),
    ag_n_tested = sum(ag_n_tested, na.rm = TRUE),
    ag_prop_resistant_conservative = if_else(ag_n_tested > 0, 
                                             (ag_n_resistant + ag_n_intermediate) / ag_n_tested, 
                                             NA_real_),
    ag_any_resistant_conservative = if_else(ag_n_tested > 0, as.integer((ag_n_resistant + ag_n_intermediate) > 0), NA_integer_),  # FIX (4)
    
    # Carbapenems
    carb_n_susceptible = sum(carb_n_susceptible, na.rm = TRUE),
    carb_n_intermediate = sum(carb_n_intermediate, na.rm = TRUE),
    carb_n_resistant = sum(carb_n_resistant, na.rm = TRUE),
    carb_n_tested = sum(carb_n_tested, na.rm = TRUE),
    carb_prop_resistant_conservative = if_else(carb_n_tested > 0, 
                                               (carb_n_resistant + carb_n_intermediate) / carb_n_tested, 
                                               NA_real_),
    carb_any_resistant_conservative = if_else(carb_n_tested > 0, as.integer((carb_n_resistant + carb_n_intermediate) > 0), NA_integer_),  # FIX (4)
    
    # Cephalosporins
    ceph_n_susceptible = sum(ceph_n_susceptible, na.rm = TRUE),
    ceph_n_intermediate = sum(ceph_n_intermediate, na.rm = TRUE),
    ceph_n_resistant = sum(ceph_n_resistant, na.rm = TRUE),
    ceph_n_tested = sum(ceph_n_tested, na.rm = TRUE),
    ceph_prop_resistant_conservative = if_else(ceph_n_tested > 0, 
                                               (ceph_n_resistant + ceph_n_intermediate) / ceph_n_tested, 
                                               NA_real_),
    ceph_any_resistant_conservative = if_else(ceph_n_tested > 0, as.integer((ceph_n_resistant + ceph_n_intermediate) > 0), NA_integer_),  # FIX (4)
    
    # Fluoroquinolones
    fq_n_susceptible = sum(fq_n_susceptible, na.rm = TRUE),
    fq_n_intermediate = sum(fq_n_intermediate, na.rm = TRUE),
    fq_n_resistant = sum(fq_n_resistant, na.rm = TRUE),
    fq_n_tested = sum(fq_n_tested, na.rm = TRUE),
    fq_prop_resistant_conservative = if_else(fq_n_tested > 0, 
                                             (fq_n_resistant + fq_n_intermediate) / fq_n_tested, 
                                             NA_real_),
    fq_any_resistant_conservative = if_else(fq_n_tested > 0, as.integer((fq_n_resistant + fq_n_intermediate) > 0), NA_integer_),  # FIX (4)
    
    # Penicillin combo
    pen_n_susceptible = sum(pen_susceptible, na.rm = TRUE),
    pen_n_intermediate = sum(pen_intermediate, na.rm = TRUE),
    pen_n_resistant = sum(pen_resistant, na.rm = TRUE),
    pen_n_tested = sum(pen_tested, na.rm = TRUE),
    pen_prop_resistant_conservative = if_else(pen_n_tested > 0, 
                                              (pen_n_resistant + pen_n_intermediate) / pen_n_tested, 
                                              NA_real_),
    pen_any_resistant_conservative = if_else(pen_n_tested > 0, as.integer((pen_n_resistant + pen_n_intermediate) > 0), NA_integer_),  # FIX (4)
    
    # Monobactam
    mono_n_susceptible = sum(mono_susceptible, na.rm = TRUE),
    mono_n_intermediate = sum(mono_intermediate, na.rm = TRUE),
    mono_n_resistant = sum(mono_resistant, na.rm = TRUE),
    mono_n_tested = sum(mono_tested, na.rm = TRUE),
    mono_prop_resistant_conservative = if_else(mono_n_tested > 0, 
                                               (mono_n_resistant + mono_n_intermediate) / mono_n_tested, 
                                               NA_real_),
    mono_any_resistant_conservative = if_else(mono_n_tested > 0, as.integer((mono_n_resistant + mono_n_intermediate) > 0), NA_integer_),  # FIX (4)
    
    # Polymyxin
    poly_n_susceptible = sum(poly_susceptible, na.rm = TRUE),
    poly_n_intermediate = sum(poly_intermediate, na.rm = TRUE),
    poly_n_resistant = sum(poly_resistant, na.rm = TRUE),
    poly_n_tested = sum(poly_tested, na.rm = TRUE),
    poly_prop_resistant_conservative = if_else(poly_n_tested > 0, 
                                               (poly_n_resistant + poly_n_intermediate) / poly_n_tested, 
                                               NA_real_),
    poly_any_resistant_conservative = if_else(poly_n_tested > 0, as.integer((poly_n_resistant + poly_n_intermediate) > 0), NA_integer_),  # FIX (4)
    
    # =========================================================================
    # OVERALL RESISTANCE (aggregated)
    # =========================================================================
    n_resistant = sum(n_resistant, na.rm = TRUE),
    n_intermediate = sum(n_intermediate, na.rm = TRUE),
    n_susceptible = sum(n_susceptible, na.rm = TRUE),
    n_tested = sum(n_tested, na.rm = TRUE),
    
    prop_resistant_strict = if_else(n_tested > 0, n_resistant / n_tested, NA_real_),
    prop_resistant_conservative = if_else(n_tested > 0, 
                                          (n_resistant + n_intermediate) / n_tested, 
                                          NA_real_),
    prop_resistant_weighted = if_else(n_tested > 0,
                                      (n_resistant + 0.5 * n_intermediate) / n_tested,
                                      NA_real_),
    
    # FIX (4): NA when no PA drug was tested at this visit
    any_resistance = if_else(n_tested > 0, as.integer(n_resistant > 0), NA_integer_),
    any_resistance_conservative = if_else(n_tested > 0,
                                          as.integer((n_resistant + n_intermediate) > 0),
                                          NA_integer_),
    
    # =========================================================================
    # ROUTE INFORMATION (aggregated - use any() because multiple isolates)
    # =========================================================================
    aminoglycoside_iv = as.integer(any(aminoglycoside_iv == 1, na.rm = TRUE)),
    aminoglycoside_oral = as.integer(any(aminoglycoside_oral == 1, na.rm = TRUE)),
    aminoglycoside_inhaled = as.integer(any(aminoglycoside_inhaled == 1, na.rm = TRUE)),
    
    carbapenem_iv = as.integer(any(carbapenem_iv == 1, na.rm = TRUE)),
    carbapenem_oral = as.integer(any(carbapenem_oral == 1, na.rm = TRUE)),
    carbapenem_inhaled = as.integer(any(carbapenem_inhaled == 1, na.rm = TRUE)),
    
    cephalosporin_iv = as.integer(any(cephalosporin_iv == 1, na.rm = TRUE)),
    cephalosporin_oral = as.integer(any(cephalosporin_oral == 1, na.rm = TRUE)),
    cephalosporin_inhaled = as.integer(any(cephalosporin_inhaled == 1, na.rm = TRUE)),
    
    fluoroquinolones_iv = as.integer(any(fluoroquinolones_iv == 1, na.rm = TRUE)),
    fluoroquinolones_oral = as.integer(any(fluoroquinolones_oral == 1, na.rm = TRUE)),
    fluoroquinolones_inhaled = as.integer(any(fluoroquinolones_inhaled == 1, na.rm = TRUE)),
    
    penicillin_iv = as.integer(any(penicillin_iv == 1, na.rm = TRUE)),
    penicillin_oral = as.integer(any(penicillin_oral == 1, na.rm = TRUE)),
    penicillin_inhaled = as.integer(any(penicillin_inhaled == 1, na.rm = TRUE)),
    
    aztreonam_iv = as.integer(any(aztreonam_iv == 1, na.rm = TRUE)),
    aztreonam_inhaled = as.integer(any(aztreonam_inhaled == 1, na.rm = TRUE)),
    
    colistin_iv = as.integer(any(colistin_iv == 1, na.rm = TRUE)),
    colistin_inhaled = as.integer(any(colistin_inhaled == 1, na.rm = TRUE)),
    
    # =========================================================================
    # EXPOSURE SCORES (use max - highest exposure level seen)
    # =========================================================================
    aminoglycoside_exposure_score = safe_max(aminoglycoside_exposure_score),
    carbapenem_exposure_score = safe_max(carbapenem_exposure_score),
    cephalosporin_exposure_score = safe_max(cephalosporin_exposure_score),
    fluoroquinolone_exposure_score = safe_max(fluoroquinolone_exposure_score),
    penicillin_exposure_score = safe_max(penicillin_exposure_score),
    aztreonam_exposure_score = safe_max(aztreonam_exposure_score),
    colistin_exposure_score = safe_max(colistin_exposure_score),
    total_antibiotic_pressure = safe_max(total_antibiotic_pressure),
    
    # =========================================================================
    # PHENOTYPE
    # =========================================================================
    # FIX (1): a visit where NO isolate has a recorded phenotype must stay
    # unknown. all(x == 1, na.rm = TRUE) on an all-NA group is all(logical(0)),
    # which is TRUE - so every phenotype-less visit was being labelled
    # "Non-mucoid" (642 of 1091 visits; ~86% of the Non-mucoid category).
    # phenotype_recorded gates all four indicators.
    phenotype_recorded = any(!is.na(pa_phenotype)),
    has_any_mucoid = if_else(phenotype_recorded,
                             as.integer(any(has_any_mucoid == 1, na.rm = TRUE)), NA_integer_),
    has_pure_mucoid = if_else(phenotype_recorded,
                              as.integer(any(has_pure_mucoid == 1, na.rm = TRUE)), NA_integer_),
    has_mixed_phenotype = if_else(phenotype_recorded,
                                  as.integer(any(has_mixed_phenotype == 1, na.rm = TRUE)), NA_integer_),
    # Only if ALL isolates with a recorded phenotype are non-mucoid
    has_nonmucoid_only = if_else(phenotype_recorded,
                                 as.integer(all(has_nonmucoid_only == 1, na.rm = TRUE)), NA_integer_),

    # Categorical version for display
    pa_phenotype = case_when(
      has_mixed_phenotype == 1 ~ "Mixed",
      has_pure_mucoid == 1 ~ "Mucoid", 
      has_nonmucoid_only == 1 ~ "Non-mucoid",
      TRUE ~ NA_character_  # Changed from "Unknown"
    ),
    
    # =========================================================================
    # ECOLOGY
    # =========================================================================
    # NB: has_staph/mrsa/mssa read the per-isolate co_organisms vector, before
    # co_organisms is reassigned below (summarise evaluates in order).
    has_staph = as.integer(any(str_detect(co_organisms, "MRSA|MSSA"), na.rm = TRUE)),
    has_mrsa = as.integer(any(str_detect(co_organisms, "MRSA"), na.rm = TRUE)),
    has_mssa = as.integer(any(str_detect(co_organisms, "MSSA"), na.rm = TRUE)),
    # FIX (5): take the union of ORGANISMS across the visit's isolates, not the
    # union of isolate-level strings. The old paste(unique(...), collapse = "; ")
    # glued whole comma-separated lists together, so a visit string like
    # "MRSA, Stenotrophomonas maltophilia" held two organisms but no ";" - and
    # n_coinfections, which counted ";", scored it as one.
    co_organisms = {
      parts <- unlist(str_split(na.omit(co_organisms), ",\\s*"))
      parts <- unique(str_trim(parts[nzchar(parts)]))
      if (length(parts) == 0) NA_character_ else paste(parts, collapse = "; ")
    },
    
    # =========================================================================
    # PRESCRIPTIONS
    # =========================================================================
    # A prescription recorded on ANY isolate at the visit counts for the visit -
    # matches the any()-based route aggregation above. n_drug_classes is
    # recomputed from these aggregated indicators after the summarise.
    #
    # but keep NA when NO isolate at the visit carries prescription
    # data. any(. == 1, na.rm = TRUE) on an all-NA group returns FALSE, which
    # silently recoded 521 of 1091 visits (all the pre-2017 ones, where the form
    # was not collected) from "unknown" to "not prescribed" and roughly halved
    # every reported prescription rate.
    across(ends_with("_prescribed"),
           ~if_else(all(is.na(.)), NA_integer_, as.integer(any(. == 1, na.rm = TRUE)))),
    
    # =========================================================================
    # CLINICAL
    # =========================================================================
    clinical_status = first(clinical_status),
    
    # =========================================================================
    # TIME
    # =========================================================================
    visit_number = first(visit_number),
    days_since_first_visit = first(days_since_first_visit),
    
    # Meta
    n_isolates = n(),
    
    .groups = "drop"
  ) %>%
  mutate(
    pa_phenotype = factor(pa_phenotype, levels = c("Mucoid", "Mixed", "Non-mucoid")),

    # -----------------------------------------------------------------------
    # Co-infection burden recomputed from the aggregated organism list so it
    # agrees with has_staph/has_mrsa/has_mssa (which use any() across isolates).
    # the separator is "; " now that co_organisms holds one organism
    # per element rather than one isolate-string per element.
    # -----------------------------------------------------------------------
    n_coinfections  = if_else(is.na(co_organisms), 0L, str_count(co_organisms, "; ") + 1L),
    has_coinfection = as.integer(n_coinfections > 0),

    # Drug-class count from the any()-aggregated prescriptions.
    # FIX (2): NA when the visit has no prescription data at all, rather than 0.
    n_drug_classes  = if_else(
      rowSums(!is.na(across(ends_with("_prescribed")))) == 0,
      NA_integer_,
      as.integer(rowSums(across(ends_with("_prescribed")), na.rm = TRUE))
    ),
    
    # -----------------------------------------------------------------------
    # RESISTANCE CONVENTION
    # Conservative (R + I) is now the primary modeled quantity, matching the
    # *_any_resistant_conservative binaries and the fig1 variance decomposition.
    # Strict (R only) counts are preserved as *_n_resistant_strict for
    # sensitivity. The *_prop_* columns were already computed correctly in the
    # summarise above and are unchanged.
    # >>> This changes the numerators used by M1, M2 and the heterogeneity
    # >>> beta-binomial from strict to conservative. Re-fit and re-verify.
    # -----------------------------------------------------------------------
    ag_n_resistant_strict   = ag_n_resistant,
    carb_n_resistant_strict = carb_n_resistant,
    ceph_n_resistant_strict = ceph_n_resistant,
    fq_n_resistant_strict   = fq_n_resistant,
    pen_n_resistant_strict  = pen_n_resistant,
    mono_n_resistant_strict = mono_n_resistant,
    poly_n_resistant_strict = poly_n_resistant,
    n_resistant_strict      = n_resistant,
    
    ag_n_resistant   = ag_n_resistant   + ag_n_intermediate,
    carb_n_resistant = carb_n_resistant + carb_n_intermediate,
    ceph_n_resistant = ceph_n_resistant + ceph_n_intermediate,
    fq_n_resistant   = fq_n_resistant   + fq_n_intermediate,
    pen_n_resistant  = pen_n_resistant  + pen_n_intermediate,
    mono_n_resistant = mono_n_resistant + mono_n_intermediate,
    poly_n_resistant = poly_n_resistant + poly_n_intermediate,
    n_resistant      = n_resistant      + n_intermediate,
    
    # now conservative; FIX (4): still NA when nothing was tested
    any_resistance = if_else(n_tested > 0, as.integer(n_resistant > 0), NA_integer_),
    
    # -----------------------------------------------------------------------
    # class-based MDR-PA at the visit level (was absent here).
    # >=3 antipseudomonal classes non-susceptible; NA if <3 classes tested.
    # -----------------------------------------------------------------------
    n_classes_tested = (ag_n_tested > 0) + (carb_n_tested > 0) +
      (ceph_n_tested > 0) + (fq_n_tested > 0) + (pen_n_tested > 0) +
      (mono_n_tested > 0) + (poly_n_tested > 0),
    # FIX (4): coalesce untested classes to 0 (see the isolate-level MDR block).
    n_classes_resistant = coalesce(ag_any_resistant_conservative, 0L) +
      coalesce(carb_any_resistant_conservative, 0L) +
      coalesce(ceph_any_resistant_conservative, 0L) +
      coalesce(fq_any_resistant_conservative, 0L) +
      coalesce(pen_any_resistant_conservative, 0L) +
      coalesce(mono_any_resistant_conservative, 0L) +
      coalesce(poly_any_resistant_conservative, 0L),
    multi_drug_resistant = case_when(
      n_classes_tested < 3     ~ NA_integer_,
      n_classes_resistant >= 3 ~ 1L,
      TRUE                     ~ 0L
    )
  )


# -----------------------------------------------------------------------------
# ALL LAGGED VARIABLES (in one place, after proper arrangement)
# -----------------------------------------------------------------------------

pa_visit_full <- pa_visit_full %>%
  group_by(individual_id) %>%
  arrange(date_collected_from_subject) %>%
  mutate(
    # =========================================================================
    # TIME
    # =========================================================================
    days_since_prev_visit = as.numeric(
      date_collected_from_subject - lag(date_collected_from_subject)
    ),
    
    # =========================================================================
    # LAGGED CLASS-SPECIFIC RESISTANCE
    # =========================================================================
    prev_ag_prop_resistant = lag(ag_prop_resistant_conservative),
    prev_ag_any_resistant = lag(ag_any_resistant_conservative),
    
    prev_carb_prop_resistant = lag(carb_prop_resistant_conservative),
    prev_carb_any_resistant = lag(carb_any_resistant_conservative),
    
    prev_ceph_prop_resistant = lag(ceph_prop_resistant_conservative),
    prev_ceph_any_resistant = lag(ceph_any_resistant_conservative),
    
    prev_fq_prop_resistant = lag(fq_prop_resistant_conservative),
    prev_fq_any_resistant = lag(fq_any_resistant_conservative),
    
    prev_pen_prop_resistant = lag(pen_prop_resistant_conservative),
    prev_pen_any_resistant = lag(pen_any_resistant_conservative),
    
    prev_mono_prop_resistant = lag(mono_prop_resistant_conservative),
    prev_mono_any_resistant = lag(mono_any_resistant_conservative),
    
    prev_poly_prop_resistant = lag(poly_prop_resistant_conservative),
    prev_poly_any_resistant = lag(poly_any_resistant_conservative),
    
    # =========================================================================
    # LAGGED OVERALL RESISTANCE
    # =========================================================================
    prev_prop_resistant_strict = lag(prop_resistant_strict),
    prev_prop_resistant_conservative = lag(prop_resistant_conservative),
    prev_prop_resistant_weighted = lag(prop_resistant_weighted),
    prev_any_resistance = lag(any_resistance),
    prev_any_resistance_conservative = lag(any_resistance_conservative),
    prev_n_resistant = lag(n_resistant),
    prev_n_intermediate = lag(n_intermediate),
    prev_n_tested = lag(n_tested),
    
    # =========================================================================
    # LAGGED ROUTE
    # =========================================================================
    prev_aminoglycoside_iv = lag(aminoglycoside_iv),
    prev_aminoglycoside_inhaled = lag(aminoglycoside_inhaled),
    
    prev_carbapenem_iv = lag(carbapenem_iv),
    
    prev_cephalosporin_iv = lag(cephalosporin_iv),
    
    prev_fluoroquinolones_iv = lag(fluoroquinolones_iv),
    prev_fluoroquinolones_oral = lag(fluoroquinolones_oral),
    
    prev_penicillin_iv = lag(penicillin_iv),
    
    prev_aztreonam_iv = lag(aztreonam_iv),
    prev_aztreonam_inhaled = lag(aztreonam_inhaled),
    
    prev_colistin_iv = lag(colistin_iv),
    prev_colistin_inhaled = lag(colistin_inhaled),
    
    # =========================================================================
    # LAGGED EXPOSURE SCORES
    # =========================================================================
    prev_aminoglycoside_exposure = lag(aminoglycoside_exposure_score),
    prev_carbapenem_exposure = lag(carbapenem_exposure_score),
    prev_cephalosporin_exposure = lag(cephalosporin_exposure_score),
    prev_fluoroquinolone_exposure = lag(fluoroquinolone_exposure_score),
    prev_penicillin_exposure = lag(penicillin_exposure_score),
    prev_total_antibiotic_pressure = lag(total_antibiotic_pressure),
    
    # =========================================================================
    # LAGGED PHENOTYPE
    # =========================================================================
    prev_has_any_mucoid = lag(has_any_mucoid),
    prev_has_pure_mucoid = lag(has_pure_mucoid),
    prev_has_mixed_phenotype = lag(has_mixed_phenotype),
    prev_pa_phenotype = lag(pa_phenotype),
    
    # =========================================================================
    # LAGGED ECOLOGY
    # =========================================================================
    prev_has_staph = lag(has_staph),
    prev_has_mrsa = lag(has_mrsa),
    prev_has_mssa = lag(has_mssa),
    prev_n_coinfections = lag(n_coinfections),
    prev_has_coinfection = lag(has_coinfection),
    
    # Staph transitions (for competition analysis)
    staph_gained = if_else(prev_has_staph == 0 & has_staph == 1, 1L, 0L),
    staph_lost = if_else(prev_has_staph == 1 & has_staph == 0, 1L, 0L),
    
    mrsa_gained = if_else(prev_has_mrsa == 0 & has_mrsa == 1, 1L, 0L),
    mrsa_lost = if_else(prev_has_mrsa == 1 & has_mrsa == 0, 1L, 0L),
    
    mssa_gained = if_else(prev_has_mssa == 0 & has_mssa == 1, 1L, 0L),
    mssa_lost = if_else(prev_has_mssa == 1 & has_mssa == 0, 1L, 0L),
    
    # MRSA ↔ MSSA switching (MGE dynamics?)
    mssa_to_mrsa = if_else(prev_has_mssa == 1 & prev_has_mrsa == 0 & 
                             has_mrsa == 1, 1L, 0L),
    mrsa_to_mssa = if_else(prev_has_mrsa == 1 & prev_has_mssa == 0 & 
                             has_mssa == 1, 1L, 0L),
    
    # Staph type at previous visit
    # FIX (6): the first visit of a patient has no previous visit, so its
    # prior-Staph state is unknown, not "None". The old `TRUE ~ "None"` caught
    # the NA lags and labelled all 92 patients' first visits "None", and
    # filter(!is.na(prev_staph_type)) downstream did not drop them - it padded
    # the "None" reference group with rows that have no prior visit at all.
    prev_staph_type = case_when(
      is.na(prev_has_mrsa) & is.na(prev_has_mssa) ~ NA_character_,
      prev_has_mrsa == 1 & prev_has_mssa == 1 ~ "Both",
      prev_has_mrsa == 1 ~ "MRSA",
      prev_has_mssa == 1 ~ "MSSA",
      TRUE ~ "None"
    ),
    
    # =========================================================================
    # LAGGED PRESCRIPTIONS
    # =========================================================================
    across(ends_with("_prescribed"), ~lag(.), .names = "prev_{.col}"),
    prev_n_drug_classes = lag(n_drug_classes),
    
    # =========================================================================
    # LAGGED CLINICAL
    # =========================================================================
    prev_clinical_status = lag(clinical_status),
    
    # =========================================================================
    # BASELINE & CHANGE
    # =========================================================================
    delta_resistance_weighted = prop_resistant_weighted - lag(prop_resistant_weighted),
    delta_resistance_conservative = prop_resistant_conservative - lag(prop_resistant_conservative)
    
  ) %>%
  ungroup()

# SAVE 

# After all cleaning
write_csv(clinical_microbe_data, "data/clinical_microbe_data_cleaned.csv")
# After PA filtering + isolate-level calculations
write_csv(pa_isolates, "data/pa_isolates.csv")
# After visit aggregation + lags
write_csv(pa_visit_full, "data/pa_visit_full.csv")

saveRDS(clinical_microbe_data, "data/clinical_microbe_data_cleaned.rds")
saveRDS(pa_isolates, "data/pa_isolates.rds")
saveRDS(pa_visit_full, "data/pa_visit_full.rds")