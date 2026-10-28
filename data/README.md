# data/

Patient data are not shared (IRB). Everything in this folder except this README is excluded from version control by `.gitignore`.

To run the analyses, place these files here.

| File | What it is |
| --- | --- |
| `CDC_114_subjects_bacterial_isolates_04-28-2023.csv` | Raw isolate export (input to `analytical_data_prep.R`) |
| `date_corrections.csv` | Typo fixes for collection dates, with columns `raw` and `corrected` |

`analytical_data_prep.R` then writes:

| File | What it is |
| --- | --- |
| `clinical_microbe_data_cleaned.{csv,rds}` | All cleaned isolate records |
| `pa_isolates.{csv,rds}` | *P. aeruginosa* isolates, one row per colony |
| `pa_visit_full.{csv,rds}` | One row per culture (patient × collection date) |
| `cfbr_id_crosswalk.{csv,rds}` | Links study IDs to registry IDs. **Linkage file: never copy it out of the secure location.** |

`sir_long.rds` is written by `AMR_sir_ordinal.R`.
