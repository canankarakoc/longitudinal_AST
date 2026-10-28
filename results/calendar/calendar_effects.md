# Calendar effects

## Sampling design by half-year

| half | patients | cultures | cultures_per_patient | colonies_tested_mean | pct_single_tested | pct_cultures_untested | median_gap_days | pct_APE_IV | pct_APE_out | pct_stable | pct_sputum | pct_mucoid_tested | pct_levo_tested | pct_aztreonam_tested |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2017-H1 | 64.00 | 141.00 | 2.20 | 1.65 | 45.39 | 0.00 | 49.00 | 39.01 | 24.82 | 32.62 | 97.87 | 53.45 | 100.00 | 100.00 |
| 2017-H2 | 64.00 | 127.00 | 1.98 | 1.71 | 42.06 | 0.79 | 63.00 | 42.52 | 21.26 | 31.50 | 98.43 | 48.15 | 100.00 | 93.52 |
| 2018-H1 | 62.00 | 146.00 | 2.35 | 1.49 | 60.00 | 0.68 | 57.50 | 45.89 | 17.81 | 34.25 | 97.95 | 51.39 | 100.00 | 81.02 |
| 2018-H2 | 65.00 | 136.00 | 2.09 | 1.47 | 61.07 | 3.68 | 83.50 | 20.59 | 34.56 | 41.91 | 96.32 | 47.15 | 99.48 | 88.08 |
| 2019-H1 | 61.00 | 141.00 | 2.31 | 1.55 | 57.58 | 6.38 | 55.00 | 24.82 | 37.59 | 37.59 | 97.87 | 50.73 | 100.00 | 85.85 |
| 2019-H2 | 66.00 | 127.00 | 1.92 | 1.50 | 61.42 | 0.00 | 77.00 | 37.80 | 21.26 | 37.01 | 100.00 | 49.21 | 99.48 | 79.06 |
| 2020-H1 | 46.00 | 63.00 | 1.37 | 1.60 | 53.97 | 0.00 | 81.00 | 22.22 | 22.22 | 53.97 | 100.00 | 46.53 | 100.00 | 77.23 |
| 2020-H2 | 25.00 | 30.00 | 1.20 | 1.37 | 66.67 | 0.00 | 229.50 | 13.33 | 3.33 | 76.67 | 93.33 | 53.66 | 97.56 | 75.61 |
| 2021-H1 | 29.00 | 38.00 | 1.31 | 1.34 | 65.79 | 0.00 | 172.50 | 23.68 | 10.53 | 63.16 | 100.00 | 58.82 | 100.00 | 76.47 |
| 2021-H2 | 33.00 | 48.00 | 1.45 | 1.21 | 80.85 | 2.08 | 168.00 | 18.75 | 8.33 | 68.75 | 100.00 | 59.65 | 100.00 | 84.21 |
| 2022-H1 | 28.00 | 41.00 | 1.46 | 1.24 | 75.61 | 0.00 | 166.50 | 17.07 | 14.63 | 63.41 | 100.00 | 47.06 | 100.00 | 66.67 |
| 2022-H2 | 30.00 | 40.00 | 1.33 | 1.25 | 80.00 | 0.00 | 168.00 | 25.00 | 17.50 | 57.50 | 100.00 | 42.00 | 100.00 | 56.00 |

## Change-point scan: best step date per drug

Max-LR null: 100 parametric simulations from the no-step model, quarterly grid.

| drug | best_change | LR | null 95% max-LR | p |
| --- | --- | --- | --- | --- |
| amikacin | 2017-07 | 14.9 | 7.7 | 0.010 |
| aztreonam | 2018-06 | 11.5 | 7.3 | 0.010 |
| cefepime | 2020-06 | 15.2 | 6.2 | 0.010 |
| ceftazidime | 2017-09 | 15.0 | 6.6 | 0.010 |
| gentamicin | 2021-12 | 6.9 | 6.3 | 0.050 |
| levofloxacin | 2018-12 | 40.1 | 8.3 | 0.010 |
| meropenem | 2017-10 | 32.7 | 7.9 | 0.010 |
| piperacillin_tazobactam | 2017-11 | 19.0 | 7.0 | 0.010 |
| tobramycin | 2020-03 | 20.2 | 7.6 | 0.010 |

## Step or drift: AIC relative to patient-only, per drug

| drug | best_change | dAIC trend | dAIC step (best) | dAIC step (CLSI) | dAIC trend+step | trend OR/yr given step | step OR given trend |
| --- | --- | --- | --- | --- | --- | --- | --- |
| levofloxacin | 2018-12 | -24.4 | -38.1 | -35.4 | -36.3 | 1.04 (p 0.687) | 2.87 (p 0.000) |
| tobramycin | 2020-03 | -13.9 | -18.3 | -1.7 | -17.0 | 0.93 (p 0.392) | 0.48 (p 0.023) |
| amikacin | 2017-07 | 1.4 | -12.9 | 1.9 | -21.4 | 0.83 (p 0.001) | 3.26 (p 0.000) |
| gentamicin | 2021-12 | 0.3 | -4.9 | 1.9 | -3.2 | 1.03 (p 0.615) | 0.42 (p 0.018) |
| meropenem | 2017-10 | -0.2 | -30.8 | -0.2 | -35.5 | 0.83 (p 0.009) | 4.04 (p 0.000) |
| cefepime | 2020-06 | 0.5 | -13.3 | 1.2 | -19.7 | 1.24 (p 0.003) | 0.25 (p 0.000) |
| ceftazidime | 2017-09 | 1.5 | -13.0 | -0.4 | -14.5 | 0.89 (p 0.061) | 2.42 (p 0.000) |
| piperacillin_tazobactam | 2017-11 | -1.9 | -17.0 | -5.2 | -16.2 | 0.93 (p 0.274) | 2.26 (p 0.000) |
| aztreonam | 2018-06 | 1.7 | -9.5 | -1.5 | -17.3 | 0.78 (p 0.002) | 2.74 (p 0.000) |

## Step odds ratios (patient random intercept)

| at | drug | change | outcome | OR |
| --- | --- | --- | --- | --- |
| best change point | amikacin | 2017-07 | ns | 2.18 (1.48–3.21) |
| best change point | amikacin | 2017-07 | s | 0.46 (0.31–0.68) |
| best change point | amikacin | 2017-07 | r | 1.83 (1.26–2.67) |
| best change point | aztreonam | 2018-06 | ns | 1.62 (1.23–2.13) |
| best change point | aztreonam | 2018-06 | s | 0.62 (0.47–0.81) |
| best change point | aztreonam | 2018-06 | r | 1.72 (1.30–2.28) |
| best change point | cefepime | 2020-06 | ns | 0.49 (0.34–0.70) |
| best change point | cefepime | 2020-06 | s | 2.05 (1.44–2.94) |
| best change point | cefepime | 2020-06 | r | 0.71 (0.49–1.03) |
| best change point | ceftazidime | 2017-09 | ns | 1.93 (1.38–2.69) |
| best change point | ceftazidime | 2017-09 | s | 0.52 (0.37–0.72) |
| best change point | ceftazidime | 2017-09 | r | 1.79 (1.27–2.52) |
| best change point | gentamicin | 2021-12 | ns | 0.47 (0.27–0.82) |
| best change point | gentamicin | 2021-12 | s | 2.11 (1.22–3.64) |
| best change point | gentamicin | 2021-12 | r | 0.90 (0.51–1.57) |
| best change point | levofloxacin | 2018-12 | ns | 3.13 (2.17–4.50) |
| best change point | levofloxacin | 2018-12 | s | 0.32 (0.22–0.46) |
| best change point | levofloxacin | 2018-12 | r | 2.81 (2.10–3.76) |
| best change point | meropenem | 2017-10 | ns | 2.80 (1.97–3.97) |
| best change point | meropenem | 2017-10 | s | 0.36 (0.25–0.51) |
| best change point | meropenem | 2017-10 | r | 4.03 (2.77–5.88) |
| best change point | piperacillin_tazobactam | 2017-11 | ns | 1.97 (1.45–2.67) |
| best change point | piperacillin_tazobactam | 2017-11 | s | 0.51 (0.37–0.69) |
| best change point | piperacillin_tazobactam | 2017-11 | r | 1.37 (1.00–1.89) |
| best change point | tobramycin | 2020-03 | ns | 0.39 (0.25–0.59) |
| best change point | tobramycin | 2020-03 | s | 2.58 (1.70–3.93) |
| best change point | tobramycin | 2020-03 | r | 0.30 (0.19–0.48) |
| CLSI 2019-01-01 | levofloxacin | 2019-01 | ns | 2.99 (2.08–4.29) |
| CLSI 2019-01-01 | levofloxacin | 2019-01 | s | 0.33 (0.23–0.48) |
| CLSI 2019-01-01 | levofloxacin | 2019-01 | r | 2.66 (1.99–3.56) |
| CLSI 2019-01-01 | tobramycin | 2019-01 | ns | 0.75 (0.56–1.00) |
| CLSI 2019-01-01 | tobramycin | 2019-01 | s | 1.33 (1.00–1.79) |
| CLSI 2019-01-01 | tobramycin | 2019-01 | r | 0.75 (0.56–1.02) |
| CLSI 2019-01-01 | amikacin | 2019-01 | ns | 1.06 (0.80–1.40) |
| CLSI 2019-01-01 | amikacin | 2019-01 | s | 0.95 (0.71–1.26) |
| CLSI 2019-01-01 | amikacin | 2019-01 | r | 1.04 (0.79–1.37) |
| CLSI 2019-01-01 | gentamicin | 2019-01 | ns | 1.05 (0.78–1.40) |
| CLSI 2019-01-01 | gentamicin | 2019-01 | s | 0.95 (0.71–1.28) |
| CLSI 2019-01-01 | gentamicin | 2019-01 | r | 1.07 (0.81–1.42) |
| CLSI 2019-01-01 | meropenem | 2019-01 | ns | 1.25 (0.93–1.68) |
| CLSI 2019-01-01 | meropenem | 2019-01 | s | 0.80 (0.59–1.07) |
| CLSI 2019-01-01 | meropenem | 2019-01 | r | 1.43 (1.06–1.91) |
| CLSI 2019-01-01 | cefepime | 2019-01 | ns | 1.13 (0.87–1.47) |
| CLSI 2019-01-01 | cefepime | 2019-01 | s | 0.88 (0.68–1.15) |
| CLSI 2019-01-01 | cefepime | 2019-01 | r | 1.29 (0.99–1.67) |
| CLSI 2019-01-01 | ceftazidime | 2019-01 | ns | 1.24 (0.95–1.61) |
| CLSI 2019-01-01 | ceftazidime | 2019-01 | s | 0.81 (0.62–1.05) |
| CLSI 2019-01-01 | ceftazidime | 2019-01 | r | 1.37 (1.04–1.80) |
| CLSI 2019-01-01 | piperacillin_tazobactam | 2019-01 | ns | 1.44 (1.11–1.87) |
| CLSI 2019-01-01 | piperacillin_tazobactam | 2019-01 | s | 0.70 (0.54–0.90) |
| CLSI 2019-01-01 | piperacillin_tazobactam | 2019-01 | r | 1.25 (0.95–1.64) |
| CLSI 2019-01-01 | aztreonam | 2019-01 | ns | 1.30 (0.99–1.71) |
| CLSI 2019-01-01 | aztreonam | 2019-01 | s | 0.77 (0.59–1.01) |
| CLSI 2019-01-01 | aztreonam | 2019-01 | r | 1.32 (1.01–1.75) |

## Levofloxacin: trend vs step

| model | AIC | OR per year | OR step | sd_patient | dAIC |
| --- | --- | --- | --- | --- | --- |
| patient only | 1214.77 |  |  | 2.96 | 58.80 |
| linear calendar trend | 1190.33 | 1.37 (1.21–1.55) |  | 3.11 | 34.35 |
| step at best date | 1176.65 |  | 3.13 (2.17–4.50) | 3.15 | 20.68 |
| step at CLSI date | 1179.40 |  | 2.99 (2.08–4.29) | 3.14 | 23.43 |
| trend + step (best) | 1178.50 | 1.04 (0.86–1.25) | 2.87 (1.64–5.00) | 3.16 | 22.52 |
| trend + step (best), adjusted | 1155.98 | 1.05 (0.86–1.27) | 3.02 (1.71–5.34) | 3.06 | 0.00 |
| trend, adjusted | 1168.69 | 1.39 (1.23–1.58) |  | 3.02 | 12.72 |

## Levofloxacin by year: cohort mean vs within-patient deviation

| year | patients | cohort_mean_NS | within_patient_deviation |
| --- | --- | --- | --- |
| 2017 | 76.00 | 0.61 | -0.06 |
| 2018 | 72.00 | 0.58 | -0.04 |
| 2019 | 73.00 | 0.66 | 0.04 |
| 2020 | 52.00 | 0.78 | 0.05 |
| 2021 | 45.00 | 0.77 | 0.08 |
| 2022 | 42.00 | 0.73 | 0.08 |

## Levofloxacin NS by year of a patient's last culture (attrition)

| last | patients | mean_NS |
| --- | --- | --- |
| 2017 | 6.00 | 0.67 |
| 2018 | 7.00 | 0.26 |
| 2019 | 12.00 | 0.29 |
| 2020 | 8.00 | 0.58 |
| 2021 | 15.00 | 0.67 |
| 2022 | 42.00 | 0.65 |

## Within-culture discordance by period (colony frame, culture-weighted)

| class | period | cultures | disc |
| --- | --- | --- | --- |
| Aminoglycoside | 2017-2018 | 259 | 0.18 (0.13–0.25) |
| Aminoglycoside | 2019 | 105 | 0.22 (0.12–0.32) |
| Aminoglycoside | 2020-2022 |  79 | 0.26 (0.14–0.37) |
| Carbapenem | 2017-2018 | 256 | 0.26 (0.18–0.34) |
| Carbapenem | 2019 | 104 | 0.32 (0.20–0.43) |
| Carbapenem | 2020-2022 |  76 | 0.26 (0.17–0.39) |
| Cephalosporin | 2017-2018 | 259 | 0.35 (0.27–0.44) |
| Cephalosporin | 2019 | 105 | 0.29 (0.19–0.39) |
| Cephalosporin | 2020-2022 |  79 | 0.31 (0.19–0.41) |
| Fluoroquinolone | 2017-2018 | 258 | 0.19 (0.13–0.26) |
| Fluoroquinolone | 2019 | 104 | 0.22 (0.14–0.32) |
| Fluoroquinolone | 2020-2022 |  78 | 0.16 (0.05–0.29) |

## Analysis-1 within-patient share, fluoroquinolone, by centring

| centring | V_tot | V_w | share | R_wp at 1 colony |
| --- | --- | --- | --- | --- |
| patient (analysis 1) | 0.101 | 0.004 | 0.04 (-0.17–0.19) | 0.04 |
| patient x side of 2018-12-01 | 0.101 | -0.001 | -0.01 (-0.21–0.12) | 0.00 |
