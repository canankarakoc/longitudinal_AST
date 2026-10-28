# =============================================================================
# PART 5: INSPECT DATA
# =============================================================================
# =============================================================================
# DESCRIPTIVE STATISTICS AND EXPLORATORY FIGURES
# =============================================================================
# =============================================================================
# SETUP
# =============================================================================
if (!file.exists("AMR_common.R")) stop("Run from the repository root (the folder that holds AMR_common.R).")

# Packages. This script was split out of the original prep file and lost its
# library() calls with the split - it referenced tibble(), ggplot(), year() and
# so on with nothing loaded, so it failed at the first tidyverse call.
library(tidyverse)
library(lubridate)
library(patchwork)
library(lme4)

# Load validated data
# read the .rds, not the .csv. The csv round-trip drops every factor -
# clinical_status came back as character with the level order set above
# ("Clinically Stable", "APE Outpatient", "APE IV Antibiotics") lost, and
# pa_phenotype likewise - so plots and tables fell back to alphabetical order.
# The csv files are still written above, for sharing outside R.
clinical_microbe_data <- readRDS("data/clinical_microbe_data_cleaned.rds")
pa_isolates           <- readRDS("data/pa_isolates.rds")
pa_visit_full         <- readRDS("data/pa_visit_full.rds")

# Create output directories
dir.create("figures", showWarnings = FALSE)
dir.create("figures/summary", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# Constants and color schemes
# -----------------------------------------------------------------------------

MAX_DAYS <- 1825  # 5 years for visualizations
REG_DAYS <- max(clinical_microbe_data$days_since_first_visit)

# Antibiotic class definitions (for pa_isolates level analysis)
ab_classes <- tibble(
  antibiotic = c("amikacin", "tobramycin", "gentamicin",
                 "meropenem", "imipenem",
                 "cefepime", "ceftazidime",
                 "ciprofloxacin", "levofloxacin",
                 "piperacillin_tazobactam",
                 "aztreonam",
                 "colistin"),
  class = c(rep("Aminoglycoside", 3),
            rep("Carbapenem", 2),
            rep("Cephalosporin", 2),
            rep("Fluoroquinolone", 2),
            "Penicillin Combo",
            "Monobactam",
            "Polymyxin")
)

# Color palettes
class_colors <- c(
  "Aminoglycoside" = "#E41A1C",
  "Carbapenem" = "#377EB8",
  "Cephalosporin" = "#4DAF4A",
  "Fluoroquinolone" = "#984EA3",
  "Penicillin Combo" = "#FF7F00",
  "Monobactam" = "#A65628",
  "Colistin" = "#F781BF"
)

clinical_colors <- c(
  "Clinically Stable" = "#4DAF4A",
  "APE Outpatient" = "#FF7F00",
  "APE IV Antibiotics" = "#E41A1C",
  "Indeterminate" = "#999999"
)

trajectory_colors <- c(
  "Persistently High" = "#E41A1C",
  "Increasing" = "#FF7F00",
  "Fluctuating" = "#984EA3",
  "Stable/Moderate" = "#4DAF4A",
  "Decreasing" = "#377EB8",
  "Persistently Low" = "#A6CEE3"
)

# =============================================================================
# PART 5: OVERALL DATA STRUCTURE
# =============================================================================
cat("Total isolate records:", nrow(clinical_microbe_data), "\n")
cat("Unique patients:", n_distinct(clinical_microbe_data$individual_id), "\n")
cat("Unique species identified:", n_distinct(clinical_microbe_data$identification_of_isolate), "\n")
cat("Date range:", 
    format(min(as.Date(clinical_microbe_data$date_collected_from_subject), na.rm = TRUE), "%Y-%m-%d"), 
    "to",
    format(max(as.Date(clinical_microbe_data$date_collected_from_subject), na.rm = TRUE), "%Y-%m-%d"), "\n")
cat("Follow-up duration: ~", 
    round(as.numeric(difftime(max(as.Date(clinical_microbe_data$date_collected_from_subject), na.rm = TRUE),
                              min(as.Date(clinical_microbe_data$date_collected_from_subject), na.rm = TRUE),
                              units = "days")) / 365, 1), "years\n\n")

# -----------------------------------------------------------------------------
# 5.1: Patient-level summary
# -----------------------------------------------------------------------------
patient_summary <- clinical_microbe_data %>%
  group_by(individual_id) %>%
  summarise(
    n_samples = n_distinct(date_collected_from_subject),
    n_isolates_total = n(),
    n_species = n_distinct(identification_of_isolate),
    follow_up_days = max(days_since_first_visit, na.rm = TRUE),
    has_pa = any(identification_of_isolate == "Pseudomonas aeruginosa"),
    .groups = "drop"
  )

cat("Samples per patient:\n")
print(summary(patient_summary$n_samples))
cat("\nIsolates per patient:\n")
print(summary(patient_summary$n_isolates_total))
cat("\nSpecies per patient:\n")
print(summary(patient_summary$n_species))
cat("\nPatients with P. aeruginosa:", sum(patient_summary$has_pa), 
    "(", round(mean(patient_summary$has_pa)*100, 1), "%)\n\n")

# -----------------------------------------------------------------------------
# 5.2: Sample-level summary
# -----------------------------------------------------------------------------
sample_summary <- clinical_microbe_data %>%
  group_by(individual_id, date_collected_from_subject, sample_type) %>%
  summarise(
    n_isolates = n(),
    n_species = n_distinct(identification_of_isolate),
    species_list = paste(unique(identification_of_isolate), collapse = "; "),
    has_multiple_species = n_species > 1,
    .groups = "drop"
  )

cat("Isolates per sample:\n")
print(table(sample_summary$n_isolates))
cat("\nSamples with multiple species:", sum(sample_summary$has_multiple_species),
    "(", round(mean(sample_summary$has_multiple_species)*100, 1), "%)\n\n")

# -----------------------------------------------------------------------------
# 5.3: Species distribution
# -----------------------------------------------------------------------------
species_counts <- clinical_microbe_data %>%
  count(identification_of_isolate, sort = TRUE)

species_prevalence <- clinical_microbe_data %>%
  group_by(identification_of_isolate) %>%
  summarise(
    n_isolates = n(),
    n_patients = n_distinct(individual_id),
    pct_patients = round(n_patients / n_distinct(clinical_microbe_data$individual_id) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_patients))

cat("Top 10 species by patient prevalence:\n")
print(species_prevalence %>% head(10))

# =============================================================================
# PART 5.4: P. AERUGINOSA PATTERNS
# =============================================================================

cat("Total PA isolates:", nrow(pa_isolates), "\n")
cat("Patients with PA:", n_distinct(pa_isolates$individual_id), "\n")
cat("PA visits:", nrow(pa_visit_full), "\n")
cat("Mean isolates per visit:", round(nrow(pa_isolates) / nrow(pa_visit_full), 2), "\n\n")

# PA isolates per patient
pa_per_patient <- pa_isolates %>%
  group_by(individual_id) %>%
  summarise(
    n_pa_isolates = n(),
    n_pa_samples = n_distinct(date_collected_from_subject),
    .groups = "drop"
  )

cat("PA samples per patient:\n")
print(summary(pa_per_patient$n_pa_samples))

# PA phenotype distribution (isolate level)
cat("\nPA Phenotype distribution (isolate level):\n")
print(pa_isolates %>% 
        filter(!is.na(pa_phenotype)) %>%
        count(pa_phenotype, sort = TRUE))

# PA phenotype distribution (visit level)
cat("\nPA Phenotype distribution (visit level):\n")
print(pa_visit_full %>% count(pa_phenotype, sort = TRUE))

# Multiple PA isolates per visit
pa_multi_isolate <- pa_visit_full %>%
  filter(n_isolates > 1) %>%
  select(individual_id, date_collected_from_subject, n_isolates, pa_phenotype)

cat("\nVisits with multiple PA isolates:", nrow(pa_multi_isolate), "\n")
if (nrow(pa_multi_isolate) > 0) {
  cat("Distribution:\n")
  print(table(pa_multi_isolate$n_isolates))
}

# =============================================================================
# PART 5.5: RESISTANCE PATTERNS
# =============================================================================

# Compare three coding approaches
resistance_comparison <- pa_visit_full %>%
  filter(!is.na(prop_resistant_strict)) %>%
  summarise(
    n_visits = n(),
    `Strict (I≠R)` = round(mean(prop_resistant_strict, na.rm = TRUE) * 100, 1),
    `Conservative (I=R)` = round(mean(prop_resistant_conservative, na.rm = TRUE) * 100, 1),
    `Weighted (I=0.5)` = round(mean(prop_resistant_weighted, na.rm = TRUE) * 100, 1)
  )

cat("Mean PA resistance by coding approach:\n")
print(resistance_comparison)
cat("\nDifference (Conservative - Strict):", 
    resistance_comparison$`Conservative (I=R)` - resistance_comparison$`Strict (I≠R)`, 
    "percentage points\n\n")

# Class-specific resistance (visit level, conservative coding)
class_resistance <- pa_visit_full %>%
  summarise(
    Aminoglycoside = round(mean(ag_any_resistant_conservative, na.rm = TRUE) * 100, 1),
    Carbapenem = round(mean(carb_any_resistant_conservative, na.rm = TRUE) * 100, 1),
    Cephalosporin = round(mean(ceph_any_resistant_conservative, na.rm = TRUE) * 100, 1),
    Fluoroquinolone = round(mean(fq_any_resistant_conservative, na.rm = TRUE) * 100, 1),
    `Penicillin Combo` = round(mean(pen_any_resistant_conservative, na.rm = TRUE) * 100, 1),
    Monobactam = round(mean(mono_any_resistant_conservative, na.rm = TRUE) * 100, 1),
    Polymyxin = round(mean(poly_any_resistant_conservative, na.rm = TRUE) * 100, 1)
  ) %>%
  pivot_longer(everything(), names_to = "Class", values_to = "% Resistant")

cat("Resistance by antibiotic class (visit level, conservative coding):\n")
print(class_resistance)

# Multi-drug resistance

mdr_summary <- pa_isolates %>%
  filter(!is.na(multi_drug_resistant)) %>%
  count(multi_drug_resistant) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat("\nMulti-drug resistance (≥3 drugs) in PA isolates:\n")
print(mdr_summary)
#61%!
# =============================================================================
# PART 5.6: CO-COLONIZATION PATTERNS
# =============================================================================

# Overall co-infection prevalence
coinfection_summary <- clinical_microbe_data %>%
  group_by(individual_id, date_collected_from_subject) %>%
  slice(1) %>%
  ungroup() %>%
  count(n_coinfections) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat("Co-infection counts (all samples):\n")
print(coinfection_summary)

# Staph co-colonization in PA visits
cat("\nStaph co-colonization in PA visits:\n")
cat("  MRSA:", sum(pa_visit_full$has_mrsa, na.rm = TRUE), 
    "(", round(mean(pa_visit_full$has_mrsa, na.rm = TRUE)*100, 1), "%)\n")
cat("  MSSA:", sum(pa_visit_full$has_mssa, na.rm = TRUE), 
    "(", round(mean(pa_visit_full$has_mssa, na.rm = TRUE)*100, 1), "%)\n")
cat("  Any Staph:", sum(pa_visit_full$has_staph, na.rm = TRUE), 
    "(", round(mean(pa_visit_full$has_staph, na.rm = TRUE)*100, 1), "%)\n\n")

# Most common co-colonizers with PA (separator is "; ")
pa_coinfections <- pa_visit_full %>%
  filter(!is.na(co_organisms), co_organisms != "") %>%
  separate_rows(co_organisms, sep = "; ") %>%
  mutate(co_organisms = str_trim(co_organisms)) %>%
  filter(co_organisms != "", co_organisms != "NA") %>%
  count(co_organisms, sort = TRUE)

cat("Most common co-colonizing organisms with PA:\n")
print(pa_coinfections %>% head(10))

# =============================================================================
# PART 5.7: PRESCRIPTION AND ROUTE DATA
# =============================================================================

# Note: ~52% missing due to data collection starting 2017
presc_availability <- pa_visit_full %>%
  summarise(
    total_visits = n(),
    has_prescription_data = sum(!is.na(aminoglycoside_prescribed)),
    pct_with_data = round(has_prescription_data / total_visits * 100, 1)
  )

cat("Prescription data availability:\n")
cat("  Total visits:", presc_availability$total_visits, "\n")
cat("  With prescription data:", presc_availability$has_prescription_data, 
    "(", presc_availability$pct_with_data, "%)\n")
cat("  Note: Data collection started 2017\n\n")

# Prescription rates (when data available)
presc_rates <- pa_visit_full %>%
  filter(!is.na(aminoglycoside_prescribed)) %>%
  summarise(
    n = n(),
    Aminoglycoside = round(mean(aminoglycoside_prescribed, na.rm = TRUE) * 100, 1),
    Carbapenem = round(mean(carbapenem_prescribed, na.rm = TRUE) * 100, 1),
    Cephalosporin = round(mean(cephalosporin_prescribed, na.rm = TRUE) * 100, 1),
    Fluoroquinolone = round(mean(fluoroquinolones_prescribed, na.rm = TRUE) * 100, 1),
    Colistin = round(mean(colistin_prescribed, na.rm = TRUE) * 100, 1)
  )

print(presc_rates)

# Route summary
route_summary <- pa_visit_full %>%
  summarise(
    `AG IV` = sum(aminoglycoside_iv == 1, na.rm = TRUE),
    `AG Inhaled` = sum(aminoglycoside_inhaled == 1, na.rm = TRUE),
    `Carb IV` = sum(carbapenem_iv == 1, na.rm = TRUE),
    `FQ IV` = sum(fluoroquinolones_iv == 1, na.rm = TRUE),
    `FQ Oral` = sum(fluoroquinolones_oral == 1, na.rm = TRUE),
    `Colistin IV` = sum(colistin_iv == 1, na.rm = TRUE),
    `Colistin Inhaled` = sum(colistin_inhaled == 1, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "Route", values_to = "N Visits")
print(route_summary)

# =============================================================================
# PART 5.8: TEMPORAL PATTERNS
# =============================================================================

# Samples per year
samples_per_year <- clinical_microbe_data %>%
  mutate(year = year(as.Date(date_collected_from_subject))) %>%
  group_by(year) %>%
  summarise(
    n_samples = n_distinct(individual_id, date_collected_from_subject),
    n_patients = n_distinct(individual_id),
    n_isolates = n(),
    .groups = "drop"
  )

print(samples_per_year)

# Visit frequency
visits_per_patient <- pa_visit_full %>%
  group_by(individual_id) %>%
  summarise(n_visits = n(), .groups = "drop")
cat("\nVisits per patient in PA data:\n")
print(summary(visits_per_patient$n_visits))

cat("\nDays between visits:\n")
print(summary(pa_visit_full$days_since_prev_visit[!is.na(pa_visit_full$days_since_prev_visit)]))

# =============================================================================
# PART 5.9: CLINICAL CONTEXT
# =============================================================================

# Clinical status distribution
cat("Clinical status distribution (PA visits):\n")
print(pa_visit_full %>% count(clinical_status, sort = TRUE))

# =============================================================================
# ECOLOGICAL INTERACTIONS
# =============================================================================

# Summarize Staph transitions
pa_visit_full %>%
  filter(!is.na(prev_has_staph)) %>%
  summarise(
    N = n(),
    `Staph Gained` = sum(staph_gained, na.rm = TRUE),
    `Staph Lost` = sum(staph_lost, na.rm = TRUE),
    `MRSA Gained` = sum(mrsa_gained, na.rm = TRUE),
    `MRSA Lost` = sum(mrsa_lost, na.rm = TRUE),
    `MSSA Gained` = sum(mssa_gained, na.rm = TRUE),
    `MSSA Lost` = sum(mssa_lost, na.rm = TRUE),
    `MSSA→MRSA` = sum(mssa_to_mrsa, na.rm = TRUE),
    `MRSA→MSSA` = sum(mrsa_to_mssa, na.rm = TRUE)
  ) %>%
  print()

# Cross-tab: Previous Staph type vs PA resistance
cat("\n\nPA resistance by previous Staph type:\n")
pa_visit_full %>%
  filter(!is.na(prev_staph_type), !is.na(any_resistance)) %>%
  group_by(prev_staph_type) %>%
  summarise(
    N = n(),
    `Mean Resistance` = mean(prop_resistant_conservative, na.rm = TRUE),
    `Any Resistance (%)` = mean(any_resistance, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(`Mean Resistance`)) %>%
  print()


# =============================================================================
# PRESCRIPTION SUMMARY - Why we chose these drug classes
# =============================================================================

# =============================================================================
# DRUG PRESCRIPTION SUMMARY
# =============================================================================
# -----------------------------------------------------------------------------
# Filter to current visit columns only (exclude prev_)
# -----------------------------------------------------------------------------

all_rx_cols <- names(pa_visit_full)[str_detect(names(pa_visit_full), "_prescribed$")]
current_rx_cols <- all_rx_cols[!str_detect(all_rx_cols, "^prev_")]

print(current_rx_cols)

# -----------------------------------------------------------------------------
# Summary by drug class (current visit)
# -----------------------------------------------------------------------------

rx_summary <- pa_visit_full %>%
  summarise(
    across(all_of(current_rx_cols), 
           list(
             n_prescribed = ~sum(. == 1, na.rm = TRUE),
             n_not_prescribed = ~sum(. == 0, na.rm = TRUE),
             pct_prescribed = ~mean(. == 1, na.rm = TRUE) * 100
           ),
           .names = "{.col}_{.fn}")
  ) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  separate(metric, into = c("drug", "stat"), sep = "_prescribed_") %>%
  pivot_wider(names_from = stat, values_from = value) %>%
  mutate(
    drug = str_replace_all(drug, "_", " ") %>% str_to_title(),
    total_with_data = n_prescribed + n_not_prescribed
  ) %>%
  arrange(desc(pct_prescribed))

# -----------------------------------------------------------------------------
# Visualization: Clean figure for data description
# -----------------------------------------------------------------------------

fig_rx_frequency <- rx_summary %>%
  mutate(Drug = fct_reorder(drug, pct_prescribed)) %>%
  ggplot(aes(x = Drug, y = pct_prescribed, fill = pct_prescribed)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(pct_prescribed, 1), "%\n(n=", n_prescribed, ")")), 
            hjust = -0.1, size = 4) +
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C", guide = "none") +
  scale_y_continuous(limits = c(0, max(rx_summary$pct_prescribed) * 1.4),
                     labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(
    title = "Antibiotic Prescription Frequency",
    subtitle = paste0("N = ", rx_summary$total_with_data[1], " visits"),
    x = NULL,
    y = "% of Visits with Prescription"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 12)
  )

ggsave("figures/summary/summaryprescription_frequency_5classes.pdf", fig_rx_frequency,
       width = 9, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# Summary table for manuscript
# -----------------------------------------------------------------------------

write_csv(rx_summary %>%
            select(Drug = drug, 
                   `N Prescribed` = n_prescribed, 
                   `N Total` = total_with_data, 
                   `% Prescribed` = pct_prescribed) %>%
            mutate(`% Prescribed` = round(`% Prescribed`, 1)),
          "results/prescription_summary_5classes.csv")

# -----------------------------------------------------------------------------
# Co-prescription patterns (which drugs are prescribed together?)
# -----------------------------------------------------------------------------

# Focus on the 5 main classes
main_rx <- c("carbapenem_prescribed", "aminoglycoside_prescribed", 
             "fluoroquinolones_prescribed", "cephalosporin_prescribed", 
             "colistin_prescribed")

# Filter to visits with any prescription data
rx_data <- pa_visit_full %>%
  select(all_of(main_rx)) %>%
  filter(if_any(everything(), ~!is.na(.)))

# Count number of drug classes prescribed per visit
rx_data <- rx_data %>%
  mutate(
    n_drugs = rowSums(across(everything(), ~replace_na(., 0))),
    any_prescribed = as.integer(n_drugs > 0)
  )

cat("Number of drug classes prescribed per visit:\n")
print(table(rx_data$n_drugs))

cat("\n% of visits with any prescription:", 
    round(mean(rx_data$any_prescribed) * 100, 1), "%\n")

# Co-prescription matrix
co_rx_matrix <- rx_data %>%
  select(all_of(main_rx)) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  {crossprod(as.matrix(.))} %>%
  as.data.frame()

rownames(co_rx_matrix) <- c("Carbapenem", "Aminoglycoside", "Fluoroquinolone", 
                            "Cephalosporin", "Colistin")
colnames(co_rx_matrix) <- rownames(co_rx_matrix)

cat("\nCo-prescription counts (diagonal = single drug):\n")
print(co_rx_matrix)

# Convert to percentages (% of visits where both drugs prescribed)
n_visits <- nrow(rx_data)
co_rx_pct <- round(co_rx_matrix / n_visits * 100, 1)

cat("\nCo-prescription percentages:\n")
print(co_rx_pct)

# -----------------------------------------------------------------------------
# Visualization: Co-prescription heatmap
# -----------------------------------------------------------------------------

co_rx_long <- co_rx_matrix %>%
  rownames_to_column("Drug1") %>%
  pivot_longer(-Drug1, names_to = "Drug2", values_to = "count") %>%
  mutate(
    pct = count / n_visits * 100,
    label = paste0(count, "\n(", round(pct, 1), "%)")
  )

fig_co_rx <- co_rx_long %>%
  ggplot(aes(x = Drug1, y = Drug2, fill = pct)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = "white", high = "#E41A1C", name = "% of visits") +
  labs(
    title = "Antibiotic Co-Prescription Patterns",
    subtitle = "Diagonal = single drug frequency | Off-diagonal = co-prescription",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave("figures/summary/summarycoprescription_heatmap.pdf", fig_co_rx,
       width = 8, height = 7, dpi = 300)

# -----------------------------------------------------------------------------
# Prescription by clinical status
# -----------------------------------------------------------------------------

rx_by_status <- pa_visit_full %>%
  filter(!is.na(clinical_status)) %>%
  group_by(clinical_status) %>%
  summarise(
    n_visits = n(),
    across(all_of(main_rx), ~mean(. == 1, na.rm = TRUE) * 100, .names = "pct_{.col}")
  ) %>%
  arrange(desc(n_visits))

print(rx_by_status %>%
        mutate(across(starts_with("pct_"), ~round(., 1))))

# Long format for plotting
rx_by_status_long <- rx_by_status %>%
  pivot_longer(starts_with("pct_"), names_to = "drug", values_to = "pct") %>%
  mutate(
    drug = str_remove(drug, "pct_") %>% str_remove("_prescribed") %>%
      str_replace_all("_", " ") %>% str_to_title()
  )

fig_rx_by_status <- rx_by_status_long %>%
  ggplot(aes(x = clinical_status, y = pct, fill = drug)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Antibiotic Prescription by Clinical Status",
    x = NULL, y = "% Prescribed", fill = "Drug Class"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave("figures/summary/summaryprescription_by_clinical_status.pdf", fig_rx_by_status,
       width = 12, height = 7, dpi = 300)


# Get all prescription columns
rx_cols <- names(clinical_microbe_data)[grep("prescribed", names(clinical_microbe_data))]

# Count prescriptions for each antibiotic
rx_summary <- clinical_microbe_data %>%
  summarise(across(all_of(rx_cols),
                   list(n  = ~sum(. == "Yes" | . == 1, na.rm = TRUE),
                        # FIX (2): isolates where the drug's prescription field
                        # was actually collected, for use as the denominator
                        nd = ~sum(!is.na(.))),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(), names_to = c("Antibiotic", ".value"), names_sep = "__") %>%
  rename(N_Prescribed = n, N_Rx_Data = nd) %>%
  mutate(
    Antibiotic = str_remove(Antibiotic, "_prescribed"),
    Antibiotic = str_replace_all(Antibiotic, "_", " "),
    Antibiotic = str_to_title(Antibiotic)
  ) %>%
  arrange(desc(N_Prescribed)) %>%
  mutate(
    Pct = round(N_Prescribed / N_Rx_Data * 100, 1),
    Drug_Class = case_when(
      str_detect(Antibiotic, "Aminoglycoside") ~ "Aminoglycoside",
      str_detect(Antibiotic, "Carbapenem") ~ "Carbapenem",
      str_detect(Antibiotic, "Cephalosporin") ~ "Cephalosporin",
      str_detect(Antibiotic, "Fluoroquinolone") ~ "Fluoroquinolone",
      str_detect(Antibiotic, "Colistin") ~ "Polymyxin",
      str_detect(Antibiotic, "Aztreonam") ~ "Monobactam",
      str_detect(Antibiotic, "Penicillin") ~ "Penicillin",
      str_detect(Antibiotic, "Tetracyclin") ~ "Tetracycline",
      str_detect(Antibiotic, "Linezolid") ~ "Oxazolidinone",
      str_detect(Antibiotic, "Vancomycin") ~ "Glycopeptide",
      str_detect(Antibiotic, "Trimethoprim|Bactrim") ~ "Sulfonamide",
      TRUE ~ "Other"
    ),
    PA_Relevant = if_else(Drug_Class %in% c("Aminoglycoside", "Carbapenem", 
                                            "Cephalosporin", "Fluoroquinolone", 
                                            "Polymyxin", "Monobactam"), 
                          "Yes", "No")
  )

cat("All antibiotic prescriptions:\n")
print(rx_summary, n = 20)

# Summary by drug class
rx_by_class <- rx_summary %>%
  group_by(Drug_Class, PA_Relevant) %>%
  summarise(
    Total_Prescribed = sum(N_Prescribed),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Prescribed))

cat("\n\nPrescriptions by drug class:\n")
print(rx_by_class)

# Justification table for manuscript
pa_drugs <- rx_summary %>%
  filter(PA_Relevant == "Yes") %>%
  select(Drug_Class, Antibiotic, N_Prescribed, Pct)

cat("Drug classes selected for analysis (PA-relevant with sufficient prescriptions):\n\n")
print(pa_drugs)

# Save for manuscript
write_csv(rx_summary, "results/prescription_summary_all.csv")
write_csv(pa_drugs, "results/prescription_summary_pa_relevant.csv")

# =============================================================================
# PART 6: FIGURES - DATA STRUCTURE
# =============================================================================
# Figure: Samples per patient
fig_samples_per_patient <- ggplot(patient_summary, aes(x = n_samples)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(patient_summary$n_samples), 
             color = "coral", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = median(patient_summary$n_samples) + 0.5, y = Inf, 
           vjust = 2, label = paste0("Median = ", median(patient_summary$n_samples)),
           color = "coral", size = 4) +
  labs(title = "Number of Sample Collection Timepoints per Patient",
       x = "Number of Samples", y = "Number of Patients") +
  theme_minimal(base_size = 12)

ggsave("figures/summary/desc_samples_per_patient.pdf", fig_samples_per_patient, 
       width = 8, height = 5, dpi = 300)

# Figure: Species distribution
fig_species <- clinical_microbe_data %>%
  count(identification_of_isolate, sort = TRUE) %>%
  head(12) %>%
  mutate(identification_of_isolate = fct_reorder(identification_of_isolate, n)) %>%
  ggplot(aes(x = n, y = identification_of_isolate)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = scales::comma(n)), hjust = -0.1, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Distribution of Bacterial Isolates",
       subtitle = paste0("Total: ", nrow(clinical_microbe_data), " isolates from ", 
                         n_distinct(clinical_microbe_data$individual_id), " patients"),
       x = "Number of Isolates", y = NULL) +
  theme_minimal(base_size = 12)

ggsave("figures/summary/fig1_species_distribution.pdf", fig_species, width = 9, height = 6, dpi = 300)

# Figure: Isolates per sample
fig_isolates_per_sample <- sample_summary %>%
  count(n_isolates) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = factor(n_isolates), y = n)) +
  geom_col(fill = "darkgreen", alpha = 0.8) +
  geom_text(aes(label = paste0(n, "\n(", round(pct, 1), "%)")), vjust = -0.5, size = 3.5) +
  labs(title = "Distribution of Isolates per Sample",
       x = "Number of Isolates", y = "Number of Samples") +
  theme_minimal(base_size = 12)

ggsave("figures/summary/desc_isolates_per_sample.pdf", fig_isolates_per_sample, 
       width = 8, height = 5, dpi = 300)

# Figure: Temporal coverage
fig_temporal_coverage <- clinical_microbe_data %>%
  mutate(year_month = floor_date(as.Date(date_collected_from_subject), "month")) %>%
  group_by(year_month) %>%
  summarise(
    n_samples = n_distinct(individual_id, date_collected_from_subject),
    n_patients = n_distinct(individual_id),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = year_month)) +
  geom_col(aes(y = n_samples), fill = "steelblue", alpha = 0.7) +
  geom_line(aes(y = n_patients * 3), color = "coral", linewidth = 1.2) +
  scale_y_continuous(
    name = "Number of Samples",
    sec.axis = sec_axis(~./3, name = "Number of Patients")
  ) +
  labs(title = "Temporal Coverage of Sample Collection",
       subtitle = "Bars = samples; Line = unique patients",
       x = "Date") +
  theme_minimal(base_size = 12) +
  theme(axis.title.y.right = element_text(color = "coral"),
        axis.text.y.right = element_text(color = "coral"))

ggsave("figures/summary/desc_temporal_coverage.pdf", fig_temporal_coverage, 
       width = 10, height = 5, dpi = 300)

# =============================================================================
# PART 6.1: FIGURES - RESISTANCE PATTERNS
# =============================================================================

# Figure: Resistance distribution - comparing coding approaches
fig2a <- pa_visit_full %>%
  filter(!is.na(prop_resistant_strict)) %>%
  ggplot(aes(x = prop_resistant_strict * 100)) +
  geom_histogram(bins = 20, fill = "#377EB8", color = "white", alpha = 0.8) +
  labs(title = "A. Strict (I ≠ R)", x = "% Resistant", y = "Count") +
  theme_minimal(base_size = 11)

fig2b <- pa_visit_full %>%
  filter(!is.na(prop_resistant_conservative)) %>%
  ggplot(aes(x = prop_resistant_conservative * 100)) +
  geom_histogram(bins = 20, fill = "#E41A1C", color = "white", alpha = 0.8) +
  labs(title = "B. Conservative (I = R)", x = "% Resistant", y = "Count") +
  theme_minimal(base_size = 11)

fig2c <- pa_visit_full %>%
  filter(!is.na(prop_resistant_weighted)) %>%
  ggplot(aes(x = prop_resistant_weighted * 100)) +
  geom_histogram(bins = 20, fill = "#4DAF4A", color = "white", alpha = 0.8) +
  labs(title = "C. Weighted (I = 0.5)", x = "% Resistant", y = "Count") +
  theme_minimal(base_size = 11)

fig2 <- (fig2a | fig2b | fig2c) +
  plot_annotation(
    title = "P. aeruginosa Resistance Distribution",
  )

ggsave("figures/summary/fig2_resistance_distribution.pdf", fig2, width = 12, height = 4, dpi = 300)


# Figure: Class-specific resistance bar chart
fig_class_bar <- pa_visit_full %>%
  summarise(
    Aminoglycoside = mean(ag_any_resistant_conservative, na.rm = TRUE),
    Carbapenem = mean(carb_any_resistant_conservative, na.rm = TRUE),
    Cephalosporin = mean(ceph_any_resistant_conservative, na.rm = TRUE),
    Fluoroquinolone = mean(fq_any_resistant_conservative, na.rm = TRUE),
    Colistin = mean(poly_any_resistant_conservative, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "class", values_to = "pct_resistant") %>%
  mutate(class = fct_reorder(class, pct_resistant)) %>%
  ggplot(aes(x = class, y = pct_resistant, fill = class)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(round(pct_resistant * 100, 1), "%")), vjust = -0.5, size = 4) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_fill_manual(values = class_colors, guide = "none") +
  labs(title = "P. aeruginosa Resistance by Antibiotic Class",
       subtitle = "Conservative coding (I = R)",
       x = NULL, y = "% Visits with Resistance, Prescribed at least 10% of the visits") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/summary/fig3_class_resistance.pdf", fig_class_bar, width = 9, height = 6, dpi = 300)

# Figure: Co-colonization
fig_coinfection <- pa_coinfections %>%
  head(10) %>%
  mutate(co_organisms = fct_reorder(co_organisms, n)) %>%
  ggplot(aes(x = n, y = co_organisms)) +
  geom_col(fill = "darkorange", alpha = 0.8) +
  geom_text(aes(label = n), hjust = -0.1, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Most Common Co-colonizing Organisms with P. aeruginosa",
       x = "Number of Visits", y = NULL) +
  theme_minimal(base_size = 12)

ggsave("figures/summary/fig4_coinfection.pdf", fig_coinfection, width = 9, height = 5, dpi = 300)

# Figure: Confounding by indication
fig_confounding <- pa_visit_full %>%
  mutate(prior_resistant = if_else(prev_any_resistance == 1, 
                                   "Prior Resistant", "Prior Susceptible")) %>%
  filter(!is.na(prior_resistant), !is.na(prev_carbapenem_prescribed)) %>%
  group_by(prior_resistant) %>%
  summarise(
    Carbapenem = mean(prev_carbapenem_prescribed, na.rm = TRUE),
    Cephalosporin = mean(prev_cephalosporin_prescribed, na.rm = TRUE),
    Aminoglycoside = mean(prev_aminoglycoside_prescribed, na.rm = TRUE),
    Fluoroquinolone = mean(prev_fluoroquinolones_prescribed, na.rm = TRUE),
    Colistin = mean(prev_colistin_prescribed, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(-prior_resistant, names_to = "Drug Class", values_to = "Prescription Rate") %>%
  ggplot(aes(x = `Drug Class`, y = `Prescription Rate`, fill = prior_resistant)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("Prior Resistant" = "coral", "Prior Susceptible" = "steelblue")) +
  labs(title = "Prescription Patterns by Prior Resistance Status",
       subtitle = "Evidence of confounding by indication",
       y = "Prescription Rate", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 15, hjust = 1))

ggsave("figures/summary/fig5_confounding.pdf", fig_confounding, width = 9, height = 5, dpi = 300)

# =============================================================================
# PART 6.2: FIGURES - TEMPORAL TRAJECTORIES
# =============================================================================

# Population means over time (6-month bins) - comparing coding approaches
temporal_resistance <- pa_visit_full %>%
  #filter(days_since_first_visit <= MAX_DAYS) %>%
  mutate(days_bin = floor(days_since_first_visit / 180) * 180) %>%
  group_by(days_bin) %>%
  summarise(
    mean_strict = mean(prop_resistant_strict, na.rm = TRUE) * 100,
    mean_conservative = mean(prop_resistant_conservative, na.rm = TRUE) * 100,
    mean_weighted = mean(prop_resistant_weighted, na.rm = TRUE) * 100,
    se_conservative = sd(prop_resistant_conservative, na.rm = TRUE) / sqrt(n()) * 100,
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

max(pa_visit_full$days_since_first_visit)

fig_temporal <- ggplot(temporal_resistance, aes(x = days_bin)) +
  geom_ribbon(aes(ymin = mean_conservative - 1.96*se_conservative,
                  ymax = mean_conservative + 1.96*se_conservative),
              fill = "#E41A1C", alpha = 0.2) +
  geom_line(aes(y = mean_strict, color = "Strict (I≠R)"), linewidth = 1.2) +
  geom_line(aes(y = mean_conservative, color = "Conservative (I=R)"), linewidth = 1.2) +
  geom_line(aes(y = mean_weighted, color = "Weighted (I=0.5)"), linewidth = 1.2) +
  geom_point(aes(y = mean_conservative, size = n), color = "#E41A1C", alpha = 0.6) +
  scale_color_manual(values = c("Strict (I≠R)" = "#377EB8", 
                                "Conservative (I=R)" = "#E41A1C",
                                "Weighted (I=0.5)" = "#4DAF4A")) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  scale_x_continuous(breaks = seq(0, REG_DAYS, by = 365),
                     labels = function(x) paste0(round(x/365), "y")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "P. aeruginosa Resistance Over Time",
       subtitle = paste0("Population means (6-month bins, n≥10, max ", round(REG_DAYS/365), " years)"),
       x = "Time Since First Visit", y = "% Resistant", color = "Coding") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figures/summary/fig6_temporal_resistance.pdf", fig_temporal, width = 10, height = 6, dpi = 300)

# Class-specific trajectories over time (using pre-calculated class columns)
class_temporal <- pa_visit_full %>%
  #filter(days_since_first_visit <= MAX_DAYS) %>%
  mutate(days_bin = floor(days_since_first_visit / 180) * 180) %>%
  group_by(days_bin) %>%
  summarise(
    Aminoglycoside = mean(ag_any_resistant_conservative, na.rm = TRUE) * 100,
    Carbapenem = mean(carb_any_resistant_conservative, na.rm = TRUE) * 100,
    Cephalosporin = mean(ceph_any_resistant_conservative, na.rm = TRUE) * 100,
    Fluoroquinolone = mean(fq_any_resistant_conservative, na.rm = TRUE) * 100,
    Colistin = mean(poly_any_resistant_conservative, na.rm = TRUE) * 100,
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 10) %>%
  pivot_longer(cols = c(Aminoglycoside, Carbapenem, Cephalosporin, Fluoroquinolone,
                        Colistin),
               names_to = "class", values_to = "pct_resistant")

fig_class_temporal <- ggplot(class_temporal, aes(x = days_bin, y = pct_resistant, color = class)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_color_manual(values = class_colors) +
  scale_x_continuous(breaks = seq(0, REG_DAYS, by = 365),
                     labels = function(x) paste0(round(x/365), "y")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "P. aeruginosa Resistance by Antibiotic Class Over Time",
       subtitle = "Population means (6-month bins, n≥10)",
       x = "Time Since First Visit", y = "% Visits with Resistance", color = "Class") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave("figures/summary/fig7_class_temporal.pdf", fig_class_temporal, width = 10, height = 6, dpi = 300)
cat("Saved: fig7_class_temporal.pdf\n")

# Faceted class trajectories with individual patient lines
patients_3plus <- pa_visit_full %>%
  #filter(days_since_first_visit <= MAX_DAYS) %>%
  group_by(individual_id) %>%
  filter(n() >= 3) %>%
  ungroup()

# Create long format for class-specific data
pa_class_long <- patients_3plus %>%
  select(individual_id, date_collected_from_subject, visit_number,
         days_since_first_visit, clinical_status,
         ag_prop_resistant_conservative, carb_prop_resistant_conservative,
         ceph_prop_resistant_conservative, fq_prop_resistant_conservative,
         pen_prop_resistant_conservative, mono_prop_resistant_conservative) %>%
  pivot_longer(cols = ends_with("_prop_resistant_conservative"),
               names_to = "class", values_to = "pct_resistant") %>%
  mutate(
    pct_resistant = pct_resistant * 100,
    class = case_when(
      str_detect(class, "^ag_") ~ "Aminoglycoside",
      str_detect(class, "^carb_") ~ "Carbapenem",
      str_detect(class, "^ceph_") ~ "Cephalosporin",
      str_detect(class, "^fq_") ~ "Fluoroquinolone",
      str_detect(class, "^pen_") ~ "Penicillin Combo",
      str_detect(class, "^mono_") ~ "Monobactam"
    )
  ) %>%
  filter(!is.na(pct_resistant))

# Population means for overlay
pop_mean_class <- pa_class_long %>%
  mutate(days_bin = floor(days_since_first_visit / 180) * 180) %>%
  group_by(class, days_bin) %>%
  summarise(
    mean_resistant = mean(pct_resistant, na.rm = TRUE),
    se = sd(pct_resistant, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 5)

fig_class_facet <- ggplot() +
  geom_line(data = pa_class_long,
            aes(x = days_since_first_visit, y = pct_resistant, group = individual_id),
            alpha = 0.15, linewidth = 0.3, color = "grey50") +
  geom_ribbon(data = pop_mean_class,
              aes(x = days_bin, 
                  ymin = pmax(mean_resistant - 1.96*se, 0),
                  ymax = pmin(mean_resistant + 1.96*se, 100),
                  fill = class),
              alpha = 0.3) +
  geom_line(data = pop_mean_class,
            aes(x = days_bin, y = mean_resistant, color = class),
            linewidth = 1.5) +
  geom_point(data = pop_mean_class,
             aes(x = days_bin, y = mean_resistant, color = class, size = n),
             alpha = 0.8) +
  facet_wrap(~class, ncol = 3) +
  scale_color_manual(values = class_colors, guide = "none") +
  scale_fill_manual(values = class_colors, guide = "none") +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, REG_DAYS, by = 365),
                     labels = function(x) paste0(round(x/365), "y")) +
  labs(title = "P. aeruginosa Resistance by Antibiotic Class",
       subtitle = paste0("Grey = individual patients; Colored = population mean ± 95% CI (n≥5, max ", 
                         round(REG_DAYS/365), " years)"),
       x = "Time Since First Visit", y = "% Resistant") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold", size = 10))

ggsave("figures/summary/fig8_class_trajectories_facet.pdf", fig_class_facet, width = 12, height = 10, dpi = 300)

# =============================================================================
# PART 6.3: FIGURES - ECOLOGICAL DYNAMICS (STAPH)
# =============================================================================

staph_prevalence <- pa_visit_full %>%
  #filter(days_since_first_visit <= MAX_DAYS) %>%
  mutate(days_bin = floor(days_since_first_visit / 180) * 180) %>%
  group_by(days_bin) %>%
  summarise(
    pct_mrsa = 100 * mean(has_mrsa, na.rm = TRUE),
    pct_mssa = 100 * mean(has_mssa, na.rm = TRUE),
    pct_any_staph = 100 * mean(has_staph, na.rm = TRUE),
    mean_resistance = 100 * mean(prop_resistant_conservative, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 10)

fig_staph <- ggplot(staph_prevalence, aes(x = days_bin)) +
  geom_line(aes(y = mean_resistance, color = "PA Resistance"), linewidth = 1.2) +
  geom_line(aes(y = pct_any_staph, color = "Any Staph"), linewidth = 1.2, linetype = "dashed") +
  geom_line(aes(y = pct_mrsa, color = "MRSA"), linewidth = 1) +
  geom_line(aes(y = pct_mssa, color = "MSSA"), linewidth = 1) +
  scale_color_manual(values = c("PA Resistance" = "#E41A1C", "Any Staph" = "#4DAF4A",
                                "MRSA" = "#984EA3", "MSSA" = "#377EB8")) +
  scale_x_continuous(breaks = seq(0, REG_DAYS, by = 365),
                     labels = function(x) paste0(round(x/365), "y")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "P. aeruginosa Resistance and Staph Co-colonization Over Time",
       x = "Time Since First Visit", y = "Percentage", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figures/summary/fig9_staph_resistance.pdf", fig_staph, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 6.4: FIGURES - ROUTE OF ADMINISTRATION
# =============================================================================

route_data <- pa_visit_full %>%
  summarise(
    `AG IV` = sum(aminoglycoside_iv == 1, na.rm = TRUE),
    `AG Inhaled` = sum(aminoglycoside_inhaled == 1, na.rm = TRUE),
    `FQ IV` = sum(fluoroquinolones_iv == 1, na.rm = TRUE),
    `FQ Oral` = sum(fluoroquinolones_oral == 1, na.rm = TRUE),
    `Carb IV` = sum(carbapenem_iv == 1, na.rm = TRUE),
    `Colistin IV` = sum(colistin_iv == 1, na.rm = TRUE),
    `Colistin Inhaled` = sum(colistin_inhaled == 1, na.rm = TRUE),
    `Aztreonam IV` = sum(aztreonam_iv == 1, na.rm = TRUE),
    `Aztreonam Inhaled` = sum(aztreonam_inhaled == 1, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "drug_route", values_to = "count") %>%
  separate(drug_route, into = c("drug", "route"), sep = " ", extra = "merge") %>%
  mutate(drug = factor(drug, levels = c("AG", "Carb", "FQ", "Colistin", "Aztreonam")))

fig_routes <- ggplot(route_data, aes(x = drug, y = count, fill = route)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("IV" = "#E41A1C", "Oral" = "#377EB8", "Inhaled" = "#4DAF4A")) +
  labs(title = "Route of Administration by Drug Class",
       subtitle = "Number of PA visits with each route",
       x = NULL, y = "Number of Visits", fill = "Route") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figures/summary/fig10_routes.pdf", fig_routes, width = 10, height = 6, dpi = 300)

# =============================================================================
# PART 7: HETEROGENEITY ANALYSIS (Updated: 5 main drug classes)
# =============================================================================

# -----------------------------------------------------------------------------
# 7.1: Variance decomposition (ICC) - Overall
# -----------------------------------------------------------------------------

var_model_overall <- lmer(
  prop_resistant_conservative ~ 1 + (1 | individual_id),
  data = pa_visit_full %>% filter(!is.na(prop_resistant_conservative))
)

var_components <- as.data.frame(VarCorr(var_model_overall))
between_var <- var_components$vcov[1]
within_var <- attr(VarCorr(var_model_overall), "sc")^2
total_var <- between_var + within_var
icc_overall <- between_var / total_var

cat("Variance decomposition (overall resistance, conservative coding):\n")
cat("  Between-patient variance:", round(between_var, 4), "\n")
cat("  Within-patient variance:", round(within_var, 4), "\n")
cat("  ICC:", round(icc_overall * 100, 1), "% of variance is between-patient\n\n")

# -----------------------------------------------------------------------------
# 7.2: ICC by antibiotic class (5 main classes)
# -----------------------------------------------------------------------------

icc_by_class <- tibble(
  class = c("Carbapenem", "Aminoglycoside", "Fluoroquinolone", 
            "Cephalosporin", "Colistin"),
  prop_col = c("carb_prop_resistant_conservative", "ag_prop_resistant_conservative",
               "fq_prop_resistant_conservative", "ceph_prop_resistant_conservative",
               "poly_prop_resistant_conservative")
) %>%
  mutate(
    icc = map_dbl(prop_col, function(col) {
      # Check if column exists
      if (!col %in% names(pa_visit_full)) {
        cat("Column not found:", col, "\n")
        return(NA_real_)
      }
      
      df <- pa_visit_full[!is.na(pa_visit_full[[col]]), ]
      cat(col, ": N =", nrow(df), "\n")
      
      if (nrow(df) < 50) return(NA_real_)
      
      m <- tryCatch({
        lmer(as.formula(paste0(col, " ~ 1 + (1 | individual_id)")), data = df)
      }, error = function(e) {
        cat("  Error:", e$message, "\n")
        NULL
      })
      
      if (is.null(m)) return(NA_real_)
      
      vc <- as.data.frame(VarCorr(m))
      between <- vc$vcov[1]
      within <- attr(VarCorr(m), "sc")^2
      between / (between + within)
    })
  )

cat("\nICC by antibiotic class (5 main classes):\n")
print(icc_by_class %>% 
        select(class, icc) %>% 
        mutate(icc_pct = paste0(round(icc * 100, 1), "%")) %>%
        arrange(desc(icc)))

# -----------------------------------------------------------------------------
# 7.3: Patient trajectory classification (unchanged)
# -----------------------------------------------------------------------------

patient_trajectories <- pa_visit_full %>%
  filter(!is.na(prop_resistant_conservative)) %>%
  group_by(individual_id) %>%
  filter(n() >= 3) %>%
  summarise(
    n_visits = n(),
    follow_up_days = max(days_since_first_visit),
    
    # Level
    mean_resistance = mean(prop_resistant_conservative, na.rm = TRUE) * 100,
    
    # Variability
    sd_resistance = sd(prop_resistant_conservative, na.rm = TRUE) * 100,
    cv = if_else(mean_resistance > 5, sd_resistance / mean_resistance * 100, NA_real_),
    
    # Change
    baseline = first(prop_resistant_conservative[days_since_first_visit == min(days_since_first_visit)]) * 100,
    final = last(prop_resistant_conservative[days_since_first_visit == max(days_since_first_visit)]) * 100,
    change = final - baseline,
    
    # Trajectory type
    trajectory_type = case_when(
      change > 20 ~ "Increasing",
      change < -20 ~ "Decreasing",
      sd_resistance > 20 ~ "Fluctuating",
      mean_resistance > 60 ~ "Persistently High",
      mean_resistance < 30 ~ "Persistently Low",
      TRUE ~ "Stable/Moderate"
    ),
    .groups = "drop"
  )

cat("\nTrajectory types:\n")
print(table(patient_trajectories$trajectory_type))

# -----------------------------------------------------------------------------
# 7.4: CV over time by class (5 main classes)
# -----------------------------------------------------------------------------

# Create long format for 5 classes
pa_class_long <- pa_visit_full %>%
  select(individual_id, days_since_first_visit,
         Carbapenem = carb_prop_resistant_conservative,
         Aminoglycoside = ag_prop_resistant_conservative,
         Fluoroquinolone = fq_prop_resistant_conservative,
         Cephalosporin = ceph_prop_resistant_conservative,
         Colistin = poly_prop_resistant_conservative) %>%
  pivot_longer(
    cols = c(Carbapenem, Aminoglycoside, Fluoroquinolone, Cephalosporin, Colistin),
    names_to = "class",
    values_to = "prop_resistant"
  ) %>%
  mutate(pct_resistant = prop_resistant * 100) %>%
  filter(!is.na(pct_resistant))

cv_by_class_time <- pa_class_long %>%
  mutate(days_bin = floor(days_since_first_visit / 180) * 180) %>%
  group_by(class, days_bin) %>%
  summarise(
    mean_resistance = mean(pct_resistant, na.rm = TRUE),
    sd_resistance = sd(pct_resistant, na.rm = TRUE),
    cv = if_else(mean_resistance > 5, sd_resistance / mean_resistance * 100, NA_real_),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 5)

# =============================================================================
# PART 8: HETEROGENEITY FIGURES (Updated: 5 main drug classes)
# =============================================================================

# Color palette for 5 classes
class_colors <- c(
  "Carbapenem" = "#E41A1C",
  "Aminoglycoside" = "#377EB8",
  "Fluoroquinolone" = "#4DAF4A",
  "Cephalosporin" = "#984EA3",
  "Colistin" = "#FF7F00"
)

# Trajectory colors
trajectory_colors <- c(
  "Persistently High" = "#E41A1C",
  "Increasing" = "#FF7F00",
  "Fluctuating" = "#FFFF33",
  "Stable/Moderate" = "#4DAF4A",
  "Decreasing" = "#377EB8",
  "Persistently Low" = "#984EA3"
)

# -----------------------------------------------------------------------------
# Figure: CV over time (5 classes)
# -----------------------------------------------------------------------------

fig_cv_time <- ggplot(cv_by_class_time, aes(x = days_bin, y = cv, color = class)) +
  geom_line(linewidth = 1.2) +
  geom_point(aes(size = n), alpha = 0.8) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  annotate("text", x = max(cv_by_class_time$days_bin, na.rm = TRUE) * 0.9, y = 105, 
           label = "CV = 100% (high heterogeneity)", size = 3, color = "grey40") +
  scale_color_manual(values = class_colors) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  scale_x_continuous(breaks = seq(0, max(cv_by_class_time$days_bin, na.rm = TRUE), by = 365),
                     labels = function(x) paste0(round(x/365), "y")) +
  labs(title = "Coefficient of Variation in Resistance Over Time",
       subtitle = "Higher CV indicates greater between-patient heterogeneity | 5 main drug classes",
       x = "Time Since First Visit", y = "CV (%)", color = "Antibiotic Class") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("figures/summary/fig_cv_over_time_5classes.pdf", fig_cv_time, width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Figure: ICC by class (5 classes)
# -----------------------------------------------------------------------------

fig_icc_class <- icc_by_class %>%
  filter(!is.na(icc)) %>%
  mutate(class = fct_reorder(class, icc)) %>%
  ggplot(aes(x = class, y = icc, fill = class)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(round(icc * 100, 1), "%")), vjust = -0.5, size = 4) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.7, y = 0.52, label = "50% threshold", size = 3, color = "grey40") +
  scale_fill_manual(values = class_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(title = "Intraclass Correlation (ICC) by Antibiotic Class",
       subtitle = "Proportion of variance due to between-patient differences | 5 main classes",
       x = NULL, y = "ICC") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/summary/fig_icc_by_class_5classes.pdf", fig_icc_class, width = 9, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Figure: Patient CV distribution (unchanged)
# -----------------------------------------------------------------------------

fig_patient_cv <- ggplot(patient_trajectories, aes(x = cv)) +
  geom_histogram(bins = 25, fill = "steelblue", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(patient_trajectories$cv, na.rm = TRUE), 
             color = "coral", linewidth = 1.2, linetype = "dashed") +
  annotate("text", x = median(patient_trajectories$cv, na.rm = TRUE) + 5, 
           y = Inf, vjust = 2, hjust = 0,
           label = paste0("Median = ", round(median(patient_trajectories$cv, na.rm = TRUE), 1), "%"),
           color = "coral", size = 4) +
  labs(title = "Distribution of Within-Patient Coefficient of Variation",
       subtitle = "Each bar represents patients with similar CV values",
       x = "CV (%) - within-patient variability over time",
       y = "Number of Patients") +
  theme_minimal(base_size = 12)

ggsave("figures/summary/fig_patient_cv_distribution.pdf", fig_patient_cv, width = 8, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# Figure: Trajectory types (unchanged)
# -----------------------------------------------------------------------------

fig_trajectory_types <- patient_trajectories %>%
  count(trajectory_type) %>%
  mutate(
    pct = n / sum(n) * 100,
    trajectory_type = factor(trajectory_type, 
                             levels = c("Persistently High", "Increasing", "Fluctuating",
                                        "Stable/Moderate", "Decreasing", "Persistently Low"))
  ) %>%
  ggplot(aes(x = trajectory_type, y = n, fill = trajectory_type)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(n, " (", round(pct, 1), "%)")), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = trajectory_colors, guide = "none") +
  labs(title = "Patient Trajectory Phenotypes",
       subtitle = paste0("Based on resistance patterns over time (n=", nrow(patient_trajectories), 
                         " patients with ≥3 visits)"),
       x = NULL, y = "Number of Patients") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figures/summary/fig_trajectory_types.pdf", fig_trajectory_types, width = 9, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# Figure: Example trajectories (unchanged)
# -----------------------------------------------------------------------------

set.seed(42)
example_patients <- patient_trajectories %>%
  group_by(trajectory_type) %>%
  slice_sample(n = 3) %>%
  ungroup() %>%
  pull(individual_id)

example_data <- pa_visit_full %>%
  filter(individual_id %in% example_patients, 
         !is.na(prop_resistant_conservative)) %>%
  left_join(patient_trajectories %>% select(individual_id, trajectory_type), by = "individual_id")

fig_example_trajectories <- ggplot(example_data,
                                   aes(x = days_since_first_visit, 
                                       y = prop_resistant_conservative * 100,
                                       group = individual_id,
                                       color = trajectory_type)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~trajectory_type, ncol = 3) +
  scale_color_manual(values = trajectory_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, max(example_data$days_since_first_visit, na.rm = TRUE), by = 730),
                     labels = function(x) paste0(round(x/365), "y")) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Example Patient Trajectories by Phenotype",
       x = "Time Since First Visit", y = "% Resistance") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))

ggsave("figures/summary/fig_example_trajectories.pdf", fig_example_trajectories, 
       width = 12, height = 8, dpi = 300)

# =============================================================================
# PART 8 continues: HETEROGENEITY BY CLINICAL FACTORS (unchanged)
# =============================================================================

# By clinical status
cv_by_status <- pa_visit_full %>%
  filter(!is.na(prop_resistant_conservative), !is.na(clinical_status)) %>%
  group_by(clinical_status) %>%
  summarise(
    mean_resistance = mean(prop_resistant_conservative, na.rm = TRUE) * 100,
    sd_resistance = sd(prop_resistant_conservative, na.rm = TRUE) * 100,
    cv = sd_resistance / mean_resistance * 100,
    n = n(),
    .groups = "drop"
  )

cat("\n\nCV by clinical status:\n")
print(cv_by_status)

# By Staph co-colonization
cv_by_staph <- pa_visit_full %>%
  filter(!is.na(prop_resistant_conservative)) %>%
  mutate(staph_status = if_else(has_staph == 1, "With Staph", "No Staph")) %>%
  group_by(staph_status) %>%
  summarise(
    mean_resistance = mean(prop_resistant_conservative, na.rm = TRUE) * 100,
    sd_resistance = sd(prop_resistant_conservative, na.rm = TRUE) * 100,
    cv = sd_resistance / mean_resistance * 100,
    n = n(),
    .groups = "drop"
  )

cat("\nCV by Staph co-colonization:\n")
print(cv_by_staph)

# Figure: Clinical heterogeneity
clinical_heterogeneity <- bind_rows(
  cv_by_status %>% mutate(factor = "Clinical Status", group = as.character(clinical_status)),
  cv_by_staph %>% mutate(factor = "Staph Co-colonization", group = staph_status)
) %>%
  filter(!is.na(cv), n >= 10)

fig_clinical_heterogeneity <- ggplot(clinical_heterogeneity, 
                                     aes(x = group, y = cv, fill = factor)) +
  geom_col(alpha = 0.8, position = "dodge") +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3) +
  facet_wrap(~factor, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Heterogeneity (CV) by Clinical Factors",
       subtitle = "Comparing variability across patient subgroups",
       x = NULL, y = "Coefficient of Variation (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(face = "bold"))

ggsave("figures/summary/fig_clinical_heterogeneity.pdf", fig_clinical_heterogeneity, 
       width = 10, height = 5, dpi = 300)

# =============================================================================
# PART 9: SAVE RESULTS
# =============================================================================

# Heterogeneity summary table
heterogeneity_summary <- tibble(
  Metric = c(
    "Overall ICC (conservative)",
    "Between-patient variance",
    "Within-patient variance",
    "Median patient CV",
    "Patients with high CV (>50%)",
    "Patients with stable trajectories",
    "Patients with increasing resistance",
    "Patients with decreasing resistance"
  ),
  Value = c(
    paste0(round(icc_overall * 100, 1), "%"),
    round(between_var, 4),
    round(within_var, 4),
    paste0(round(median(patient_trajectories$cv, na.rm = TRUE), 1), "%"),
    paste0(sum(patient_trajectories$cv > 50, na.rm = TRUE), 
           " (", round(mean(patient_trajectories$cv > 50, na.rm = TRUE) * 100, 1), "%)"),
    paste0(sum(patient_trajectories$trajectory_type %in% c("Stable/Moderate", "Persistently Low", "Persistently High")),
           " (", round(mean(patient_trajectories$trajectory_type %in% c("Stable/Moderate", "Persistently Low", "Persistently High")) * 100, 1), "%)"),
    paste0(sum(patient_trajectories$trajectory_type == "Increasing"),
           " (", round(mean(patient_trajectories$trajectory_type == "Increasing") * 100, 1), "%)"),
    paste0(sum(patient_trajectories$trajectory_type == "Decreasing"),
           " (", round(mean(patient_trajectories$trajectory_type == "Decreasing") * 100, 1), "%)")
  ),
  Interpretation = c(
    "Proportion of variance due to between-patient differences",
    "Variance attributable to patient-level factors",
    "Variance attributable to within-patient changes over time",
    "Typical within-patient variability",
    "Patients with highly variable resistance over time",
    "Patients with consistent resistance patterns",
    "Patients showing resistance emergence",
    "Patients showing resistance loss"
  )
)

write_csv(heterogeneity_summary, "results/heterogeneity_summary.csv")
write_csv(patient_trajectories, "results/patient_trajectory_phenotypes.csv")
write_csv(icc_by_class, "results/icc_by_antibiotic_class_5classes.csv")

# -----------------------------------------------------------------------------
# COMBINED HETEROGENEITY FIGURE
# -----------------------------------------------------------------------------

combined_heterogeneity <- (fig_cv_time | fig_icc_class) /
  (fig_patient_cv | fig_trajectory_types) +
  plot_annotation(
    title = "Heterogeneity in P. aeruginosa Resistance Patterns",
    subtitle = "5 main drug classes: Carbapenem, Aminoglycoside, Fluoroquinolone, Cephalosporin, Colistin",
    tag_levels = "A",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave("figures/summary/fig_combined_heterogeneity_5classes.pdf", combined_heterogeneity, 
       width = 14, height = 10, dpi = 300)

# =============================================================================
# RESISTANCE TESTING DATA AVAILABILITY BY DRUG CLASS
# =============================================================================

# -----------------------------------------------------------------------------
# Summarize resistance status counts for each drug class
# -----------------------------------------------------------------------------

# Define the resistance variables for each class
resistance_vars <- tibble(
  
  drug_class = c("Carbapenem", "Aminoglycoside", "Fluoroquinolone", 
                 "Cephalosporin", "Colistin"),
  resist_var = c("carb_any_resistant_conservative", "ag_any_resistant_conservative",
                 "fq_any_resistant_conservative", "ceph_any_resistant_conservative",
                 "poly_any_resistant_conservative"),
  prev_resist_var = c("prev_carb_any_resistant", "prev_ag_any_resistant",
                      "prev_fq_any_resistant", "prev_ceph_any_resistant",
                      "prev_poly_any_resistant"),
  rx_var = c("carbapenem_prescribed", "aminoglycoside_prescribed",
             "fluoroquinolones_prescribed", "cephalosporin_prescribed",
             "colistin_prescribed")
)

# Calculate counts for current visit resistance
resistance_counts <- pmap_dfr(resistance_vars, function(drug_class, resist_var, prev_resist_var, rx_var) {
  
  pa_visit_full %>%
    summarise(
      Drug_Class = drug_class,
      
      # Current visit resistance
      N_Resistant = sum(.data[[resist_var]] == 1, na.rm = TRUE),
      N_Susceptible = sum(.data[[resist_var]] == 0, na.rm = TRUE),
      N_NA = sum(is.na(.data[[resist_var]])),
      N_Total = n(),
      
      # Percentages
      Pct_Resistant = round(N_Resistant / (N_Resistant + N_Susceptible) * 100, 1),
      Pct_Tested = round((N_Resistant + N_Susceptible) / N_Total * 100, 1),
      
      # Prescription frequency.
      # FIX (2): the denominator is visits where prescription data was actually
      # collected, not all visits. Prescriptions are NA (not 0) for the pre-2017
      # visits now, so N_Total would dilute the rate by roughly half.
      N_Prescribed = sum(.data[[rx_var]] == 1, na.rm = TRUE),
      N_Rx_Data = sum(!is.na(.data[[rx_var]])),
      Pct_Prescribed = round(N_Prescribed / N_Rx_Data * 100, 1)
    )
})

print(resistance_counts %>%
        select(Drug_Class, N_Resistant, N_Susceptible, N_NA, Pct_Tested, 
               Pct_Resistant, N_Prescribed, Pct_Prescribed) %>%
        arrange(desc(Pct_Tested)))

# -----------------------------------------------------------------------------
# Calculate transition counts (for Markov/GLMM stratified models)
# -----------------------------------------------------------------------------

transition_counts <- pmap_dfr(resistance_vars, function(drug_class, resist_var, prev_resist_var, rx_var) {
  
  # Need both previous and current resistance to calculate transitions
  if (!prev_resist_var %in% names(pa_visit_full)) {
    return(tibble(Drug_Class = drug_class, Note = "No lagged variable"))
  }
  
  pa_visit_full %>%
    filter(!is.na(.data[[prev_resist_var]]), !is.na(.data[[resist_var]])) %>%
    summarise(
      Drug_Class = drug_class,
      
      # Starting susceptible
      N_S_start = sum(.data[[prev_resist_var]] == 0),
      N_S_to_S = sum(.data[[prev_resist_var]] == 0 & .data[[resist_var]] == 0),
      N_S_to_R = sum(.data[[prev_resist_var]] == 0 & .data[[resist_var]] == 1),
      
      # Starting resistant
      N_R_start = sum(.data[[prev_resist_var]] == 1),
      N_R_to_S = sum(.data[[prev_resist_var]] == 1 & .data[[resist_var]] == 0),
      N_R_to_R = sum(.data[[prev_resist_var]] == 1 & .data[[resist_var]] == 1),
      
      # Total with transition data
      N_Transitions = N_S_start + N_R_start
    )
})

print(transition_counts %>%
        select(Drug_Class, N_Transitions, N_S_start, N_S_to_R, N_R_start, N_R_to_S) %>%
        mutate(
          Acquisition_Events = N_S_to_R,
          Reversion_Events = N_R_to_S,
          Acquisition_Rate = paste0(round(N_S_to_R / N_S_start * 100, 1), "%"),
          Reversion_Rate = paste0(round(N_R_to_S / N_R_start * 100, 1), "%")
        ) %>%
        arrange(desc(N_Transitions)))

# -----------------------------------------------------------------------------
# Assessment for model viability
# -----------------------------------------------------------------------------

model_viability <- transition_counts %>%
  mutate(
    # Minimum thresholds for reliable inference
    Acquisition_OK = N_S_to_R >= 20,
    Persistence_OK = N_R_to_R >= 20,
    Reversion_OK = N_R_to_S >= 10,
    
    # Overall assessment
    GLMM_Stratified = case_when(
      N_S_to_R >= 20 & N_R_start >= 30 ~ "✓ Viable",
      N_S_to_R >= 10 & N_R_start >= 20 ~ "⚠ Marginal",
      TRUE ~ "✗ Insufficient"
    ),
    Markov_Analysis = case_when(
      N_Transitions >= 100 & N_R_start >= 20 ~ "✓ Viable",
      N_Transitions >= 50 & N_R_start >= 10 ~ "⚠ Marginal",
      TRUE ~ "✗ Insufficient"
    ),
    Cross_Resistance = case_when(
      N_Transitions >= 100 ~ "✓ Viable",
      N_Transitions >= 50 ~ "⚠ Marginal",
      TRUE ~ "✗ Insufficient"
    )
  ) %>%
  select(Drug_Class, N_Transitions, N_S_to_R, N_R_start, 
         GLMM_Stratified, Markov_Analysis, Cross_Resistance)

print(model_viability)

# -----------------------------------------------------------------------------
# Visualization: Data availability
# -----------------------------------------------------------------------------

fig_data_availability <- resistance_counts %>%
  pivot_longer(cols = c(N_Resistant, N_Susceptible, N_NA),
               names_to = "Status", values_to = "Count") %>%
  mutate(
    Status = recode(Status,
                    "N_Resistant" = "Resistant",
                    "N_Susceptible" = "Susceptible",
                    "N_NA" = "Not Tested"),
    Status = factor(Status, levels = c("Not Tested", "Susceptible", "Resistant")),
    Drug_Class = fct_reorder(Drug_Class, -Pct_Tested)
  ) %>%
  ggplot(aes(x = Drug_Class, y = Count, fill = Status)) +
  geom_col(width = 0.7) +
  geom_text(data = resistance_counts,
            aes(x = Drug_Class, y = N_Total + 20, 
                label = paste0(Pct_Tested, "% tested")),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("Not Tested" = "grey70", 
                               "Susceptible" = "#377EB8", 
                               "Resistant" = "#E41A1C")) +
  labs(
    title = "Antibiotic Susceptibility Testing Data Availability",
    subtitle = "Colistin: heavily prescribed but rarely tested for resistance",
    x = NULL, y = "Number of Visits", fill = "Resistance Status"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave("figures/summary/resistance_data_availability.pdf", fig_data_availability,
       width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# Save summary
# -----------------------------------------------------------------------------

write_csv(resistance_counts, "results/resistance_data_summary.csv")
write_csv(model_viability, "results/model_viability_assessment.csv")