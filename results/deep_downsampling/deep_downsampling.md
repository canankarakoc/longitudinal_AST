# Downsampling the Vanderwoude et al. 2024 deep samples

75 isolates per sputum, triplicate discs; D = 5000 draws per n. CLSI calls on mean zone.

## Reference populations

| sample | drug | S | I | R | NS | plurality | heterogeneous |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P1 (cohort id 7) | amikacin | 43% | 47% | 11% | 57% | I | yes |
| P1 (cohort id 7) | ceftazidime | 100% | 0% | 0% | 0% | S |  |
| P1 (cohort id 7) | ciprofloxacin | 0% | 5% | 95% | 100% | R | yes |
| P1 (cohort id 7) | meropenem | 100% | 0% | 0% | 0% | S |  |
| P1 (cohort id 7) | piperacillin-tazobactam | 100% | 0% | 0% | 0% | S |  |
| P1 (cohort id 7) | tobramycin | 0% | 0% | 100% | 100% | R |  |
| P2 (not in cohort) | amikacin | 100% | 0% | 0% | 0% | S |  |
| P2 (not in cohort) | ceftazidime | 100% | 0% | 0% | 0% | S |  |
| P2 (not in cohort) | ciprofloxacin | 100% | 0% | 0% | 0% | S |  |
| P2 (not in cohort) | meropenem | 100% | 0% | 0% | 0% | S |  |
| P2 (not in cohort) | piperacillin-tazobactam | 100% | 0% | 0% | 0% | S |  |
| P2 (not in cohort) | tobramycin | 100% | 0% | 0% | 0% | S |  |
| P3 (cohort id 38) | amikacin | 100% | 0% | 0% | 0% | S |  |
| P3 (cohort id 38) | ceftazidime | 93% | 7% | 0% | 7% | S | yes |
| P3 (cohort id 38) | ciprofloxacin | 0% | 24% | 76% | 100% | R | yes |
| P3 (cohort id 38) | meropenem | 100% | 0% | 0% | 0% | S |  |
| P3 (cohort id 38) | piperacillin-tazobactam | 96% | 4% | 0% | 4% | S | yes |
| P3 (cohort id 38) | tobramycin | 100% | 0% | 0% | 0% | S |  |
| P4 (not in cohort) | amikacin | 100% | 0% | 0% | 0% | S |  |
| P4 (not in cohort) | ceftazidime | 100% | 0% | 0% | 0% | S |  |
| P4 (not in cohort) | ciprofloxacin | 99% | 1% | 0% | 1% | S |  |
| P4 (not in cohort) | meropenem | 100% | 0% | 0% | 0% | S |  |
| P4 (not in cohort) | piperacillin-tazobactam | 100% | 0% | 0% | 0% | S |  |
| P4 (not in cohort) | tobramycin | 100% | 0% | 0% | 0% | S |  |

## Assay vs colony discordance (every cell)

| sample | drug | same isolate, two discs | two isolates, mean zone | two isolates, one disc each |
| --- | --- | --- | --- | --- |
| P1 (cohort id 7) | amikacin | 27% | 60% | 63% |
| P1 (cohort id 7) | ceftazidime | 0% | 0% | 0% |
| P1 (cohort id 7) | ciprofloxacin | 6% | 10% | 8% |
| P1 (cohort id 7) | meropenem | 0% | 0% | 0% |
| P1 (cohort id 7) | piperacillin-tazobactam | 0% | 0% | 0% |
| P1 (cohort id 7) | tobramycin | 0% | 0% | 0% |
| P2 (not in cohort) | amikacin | 0% | 0% | 0% |
| P2 (not in cohort) | ceftazidime | 0% | 0% | 0% |
| P2 (not in cohort) | ciprofloxacin | 0% | 0% | 0% |
| P2 (not in cohort) | meropenem | 0% | 0% | 0% |
| P2 (not in cohort) | piperacillin-tazobactam | 0% | 0% | 0% |
| P2 (not in cohort) | tobramycin | 0% | 0% | 0% |
| P3 (cohort id 38) | amikacin | 1% | 0% | 1% |
| P3 (cohort id 38) | ceftazidime | 4% | 13% | 11% |
| P3 (cohort id 38) | ciprofloxacin | 9% | 37% | 32% |
| P3 (cohort id 38) | meropenem | 0% | 0% | 0% |
| P3 (cohort id 38) | piperacillin-tazobactam | 4% | 8% | 9% |
| P3 (cohort id 38) | tobramycin | 0% | 0% | 0% |
| P4 (not in cohort) | amikacin | 3% | 0% | 3% |
| P4 (not in cohort) | ceftazidime | 0% | 0% | 0% |
| P4 (not in cohort) | ciprofloxacin | 3% | 3% | 6% |
| P4 (not in cohort) | meropenem | 0% | 0% | 0% |
| P4 (not in cohort) | piperacillin-tazobactam | 0% | 0% | 0% |
| P4 (not in cohort) | tobramycin | 0% | 0% | 0% |

## Heterogeneous cells: estimands by colonies tested (colony sampling only / + single-disc noise)

| sample | drug | n | detect | |freq err| | dominant correct | switch, plurality | switch, most-resistant |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P1 (cohort id 7) | amikacin |  1 | 57% | 0.49 | 47% / 39% | 58% / 63% | 58% / 63% |
| P1 (cohort id 7) | amikacin |  2 | 82% | 0.28 | 47% / 39% | 59% / 64% | 55% / 63% |
| P1 (cohort id 7) | amikacin |  4 | 97% | 0.20 | 50% / 40% | 53% / 58% | 52% / 54% |
| P1 (cohort id 7) | amikacin |  8 | 100% | 0.14 | 55% / 41% | 51% / 54% | 48% / 34% |
| P1 (cohort id 7) | amikacin | 15 | 100% | 0.09 | 56% / 39% | 50% / 50% | 25% / 9% |
| P1 (cohort id 7) | ciprofloxacin |  1 | 100% | 0.00 | 94% / 96% | 11% / 7% | 11% / 7% |
| P1 (cohort id 7) | ciprofloxacin |  2 | 100% | 0.00 | 94% / 96% | 11% / 8% | 1% / 0% |
| P1 (cohort id 7) | ciprofloxacin |  4 | 100% | 0.00 | 99% / 100% | 1% / 1% | 0% / 0% |
| P1 (cohort id 7) | ciprofloxacin |  8 | 100% | 0.00 | 100% / 100% | 0% / 0% | 0% / 0% |
| P1 (cohort id 7) | ciprofloxacin | 15 | 100% | 0.00 | 100% / 100% | 0% / 0% | 0% / 0% |
| P3 (cohort id 38) | ceftazidime |  1 | 7% | 0.13 | 93% / 94% | 13% / 11% | 13% / 11% |
| P3 (cohort id 38) | ceftazidime |  2 | 13% | 0.12 | 93% / 95% | 12% / 10% | 22% / 20% |
| P3 (cohort id 38) | ceftazidime |  4 | 25% | 0.10 | 99% / 99% | 2% / 2% | 37% / 35% |
| P3 (cohort id 38) | ceftazidime |  8 | 44% | 0.08 | 100% / 100% | 0% / 0% | 49% / 52% |
| P3 (cohort id 38) | ceftazidime | 15 | 70% | 0.04 | 100% / 100% | 0% / 0% | 43% / 60% |
| P3 (cohort id 38) | ciprofloxacin |  1 | 100% | 0.00 | 76% / 80% | 37% / 32% | 37% / 32% |
| P3 (cohort id 38) | ciprofloxacin |  2 | 100% | 0.00 | 75% / 81% | 37% / 31% | 11% / 7% |
| P3 (cohort id 38) | ciprofloxacin |  4 | 100% | 0.00 | 85% / 90% | 24% / 18% | 1% / 0% |
| P3 (cohort id 38) | ciprofloxacin |  8 | 100% | 0.00 | 95% / 97% | 10% / 5% | 0% / 0% |
| P3 (cohort id 38) | ciprofloxacin | 15 | 100% | 0.00 | 99% / 100% | 1% / 0% | 0% / 0% |
| P3 (cohort id 38) | piperacillin-tazobactam |  1 | 4% | 0.07 | 96% / 96% | 8% / 9% | 8% / 9% |
| P3 (cohort id 38) | piperacillin-tazobactam |  2 | 8% | 0.08 | 96% / 95% | 8% / 9% | 15% / 16% |
| P3 (cohort id 38) | piperacillin-tazobactam |  4 | 15% | 0.07 | 100% / 100% | 1% / 1% | 27% / 28% |
| P3 (cohort id 38) | piperacillin-tazobactam |  8 | 30% | 0.06 | 100% / 100% | 0% / 0% | 41% / 43% |
| P3 (cohort id 38) | piperacillin-tazobactam | 15 | 50% | 0.04 | 100% / 100% | 0% / 0% | 50% / 50% |

## Colonies needed (heterogeneous cells)

detect >= 95%; dominant call right >= 90%; apparent-switch probability <= 10%. NA = not reached by 30.

| sample | drug | mode | n_detect_95 | n_dominant_90 | n_switch_plur_10 | n_switch_worst_10 |
| --- | --- | --- | --- | --- | --- | --- |
| P1 (cohort id 7) | amikacin | mean | 4 | NA | NA | 30 |
| P1 (cohort id 7) | amikacin | single | 4 | NA | NA | 15 |
| P1 (cohort id 7) | ciprofloxacin | mean | 1 | 1 | 3 | 2 |
| P1 (cohort id 7) | ciprofloxacin | single | 1 | 1 | 1 | 1 |
| P3 (cohort id 38) | ceftazidime | mean | NA | 1 | 3 | NA |
| P3 (cohort id 38) | ceftazidime | single | NA | 1 | 3 | NA |
| P3 (cohort id 38) | ciprofloxacin | mean | 1 | 5 | 8 | 3 |
| P3 (cohort id 38) | ciprofloxacin | single | 1 | 3 | 5 | 2 |
| P3 (cohort id 38) | piperacillin-tazobactam | mean | NA | 1 | 1 | 1 |
| P3 (cohort id 38) | piperacillin-tazobactam | single | NA | 1 | 1 | 1 |
