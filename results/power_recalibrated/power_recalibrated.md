# Power, recalibrated

2000 reps per cell (MC SE <= 0.011). Aminoglycoside calibration unless stated.
Effect sizes: step = 1 logit; coupling beta = 0.2 at 84 d (kappa = 4.52e-03 /day, tau = 131 d).
False-positive rate at baseline (2000 reps): step 0.056, coupling 0.053.

## Matched 2x levers (Table 2, recalibrated)

| budget | lever | N | L | dt | m | followup_days | step | gain step | coupling | gain coupling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| free | baseline (as run) |  90 | 12 | 84 | 1.5 |  924 | 0.44 | +0.00 | 0.42 | +0.00 |
| free | cohort 90 -> 180 | 180 | 12 | 84 | 1.5 |  924 | 0.73 | +0.30 | 0.70 | +0.28 |
| free | visits 12 -> 24, same follow-up |  90 | 24 | 40 | 1.5 |  920 | 0.52 | +0.08 | 0.57 | +0.15 |
| free | follow-up doubled (24 visits, 84 d) |  90 | 24 | 84 | 1.5 | 1932 | 0.66 | +0.23 | 0.81 | +0.39 |
| free | depth 1.5 -> 3 colonies |  90 | 12 | 84 | 3.0 |  924 | 0.49 | +0.05 | 0.59 | +0.17 |
| free | interval 84 -> 42 d, same visits |  90 | 12 | 42 | 1.5 |  462 | 0.40 | -0.04 | 0.23 | -0.19 |
| fixed cost (base) | baseline (as run) |  90 | 12 | 84 | 1.5 |  924 | 0.43 | +0.00 | 0.41 | +0.00 |
| fixed cost (base) | visits 12 -> 24, same follow-up |  53 | 24 | 40 | 1.5 |  920 | 0.32 | -0.11 | 0.37 | -0.04 |
| fixed cost (base) | follow-up doubled (24 visits, 84 d) |  53 | 24 | 84 | 1.5 | 1932 | 0.44 | +0.02 | 0.56 | +0.15 |
| fixed cost (base) | depth 1.5 -> 3 colonies |  69 | 12 | 84 | 3.0 |  924 | 0.40 | -0.02 | 0.46 | +0.05 |
| fixed cost (base) | interval 84 -> 42 d, same visits |  90 | 12 | 42 | 1.5 |  462 | 0.38 | -0.04 | 0.24 | -0.17 |

## All lever cells

| params | budget | lever | N | L | dt | m | followup_days | step | coupling |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Aminoglycoside (fitted tau) | N fixed at 90 | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.46 | 0.08 |
| Aminoglycoside (fitted tau) | N fixed at 90 | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.40 | 0.13 |
| Aminoglycoside (fitted tau) | N fixed at 90 | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.40 | 0.25 |
| Aminoglycoside (fitted tau) | N fixed at 90 | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.42 | 0.32 |
| Aminoglycoside (fitted tau) | N fixed at 90 | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.45 | 0.42 |
| Aminoglycoside (fitted tau) | N fixed at 90 | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.49 | 0.46 |
| Aminoglycoside (fitted tau) | N fixed at 90 | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.51 | 0.41 |
| Aminoglycoside (fitted tau) | N fixed at 90 | visits (follow-up fixed) |  90 |  4 | 308 | 1.5 |  924 | 0.27 | 0.08 |
| Aminoglycoside (fitted tau) | N fixed at 90 | visits (follow-up fixed) |  90 |  6 | 185 | 1.5 |  924 | 0.34 | 0.17 |
| Aminoglycoside (fitted tau) | N fixed at 90 | visits (follow-up fixed) |  90 |  8 | 132 | 1.5 |  924 | 0.38 | 0.29 |
| Aminoglycoside (fitted tau) | N fixed at 90 | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.43 | 0.45 |
| Aminoglycoside (fitted tau) | N fixed at 90 | visits (follow-up fixed) |  90 | 16 |  62 | 1.5 |  924 | 0.47 | 0.49 |
| Aminoglycoside (fitted tau) | N fixed at 90 | visits (follow-up fixed) |  90 | 24 |  40 | 1.5 |  924 | 0.51 | 0.57 |
| Aminoglycoside (fitted tau) | N fixed at 90 | visits (follow-up fixed) |  90 | 32 |  30 | 1.5 |  924 | 0.54 | 0.62 |
| Aminoglycoside (fitted tau) | N fixed at 90 | visits (follow-up fixed) |  90 | 48 |  20 | 1.5 |  924 | 0.53 | 0.65 |
| Aminoglycoside (fitted tau) | N fixed at 90 | follow-up (interval fixed) |  90 |  4 |  84 | 1.5 |  252 | 0.26 | 0.08 |
| Aminoglycoside (fitted tau) | N fixed at 90 | follow-up (interval fixed) |  90 |  6 |  84 | 1.5 |  420 | 0.31 | 0.14 |
| Aminoglycoside (fitted tau) | N fixed at 90 | follow-up (interval fixed) |  90 |  8 |  84 | 1.5 |  588 | 0.35 | 0.21 |
| Aminoglycoside (fitted tau) | N fixed at 90 | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.46 | 0.42 |
| Aminoglycoside (fitted tau) | N fixed at 90 | follow-up (interval fixed) |  90 | 16 |  84 | 1.5 | 1260 | 0.52 | 0.56 |
| Aminoglycoside (fitted tau) | N fixed at 90 | follow-up (interval fixed) |  90 | 24 |  84 | 1.5 | 1932 | 0.68 | 0.81 |
| Aminoglycoside (fitted tau) | N fixed at 90 | follow-up (interval fixed) |  90 | 32 |  84 | 1.5 | 2604 | 0.78 | 0.92 |
| Aminoglycoside (fitted tau) | N fixed at 90 | follow-up (interval fixed) |  90 | 48 |  84 | 1.5 | 3948 | 0.89 | 0.98 |
| Aminoglycoside (fitted tau) | N fixed at 90 | depth |  90 | 12 |  84 | 1.0 |  924 | 0.42 | 0.35 |
| Aminoglycoside (fitted tau) | N fixed at 90 | depth |  90 | 12 |  84 | 1.5 |  924 | 0.44 | 0.44 |
| Aminoglycoside (fitted tau) | N fixed at 90 | depth |  90 | 12 |  84 | 2.0 |  924 | 0.45 | 0.50 |
| Aminoglycoside (fitted tau) | N fixed at 90 | depth |  90 | 12 |  84 | 3.0 |  924 | 0.47 | 0.56 |
| Aminoglycoside (fitted tau) | N fixed at 90 | depth |  90 | 12 |  84 | 4.0 |  924 | 0.51 | 0.58 |
| Aminoglycoside (fitted tau) | N fixed at 90 | depth |  90 | 12 |  84 | 6.0 |  924 | 0.52 | 0.64 |
| Aminoglycoside (fitted tau) | N fixed at 90 | depth |  90 | 12 |  84 | 8.0 |  924 | 0.51 | 0.64 |
| Aminoglycoside (fitted tau) | fixed cost: base | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.45 | 0.07 |
| Aminoglycoside (fitted tau) | fixed cost: base | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.39 | 0.14 |
| Aminoglycoside (fitted tau) | fixed cost: base | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.40 | 0.24 |
| Aminoglycoside (fitted tau) | fixed cost: base | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.41 | 0.30 |
| Aminoglycoside (fitted tau) | fixed cost: base | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.45 | 0.42 |
| Aminoglycoside (fitted tau) | fixed cost: base | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.48 | 0.47 |
| Aminoglycoside (fitted tau) | fixed cost: base | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.53 | 0.44 |
| Aminoglycoside (fitted tau) | fixed cost: base | visits (follow-up fixed) | 164 |  4 | 308 | 1.5 |  924 | 0.43 | 0.09 |
| Aminoglycoside (fitted tau) | fixed cost: base | visits (follow-up fixed) | 136 |  6 | 185 | 1.5 |  924 | 0.49 | 0.24 |
| Aminoglycoside (fitted tau) | fixed cost: base | visits (follow-up fixed) | 116 |  8 | 132 | 1.5 |  924 | 0.48 | 0.33 |
| Aminoglycoside (fitted tau) | fixed cost: base | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.45 | 0.42 |
| Aminoglycoside (fitted tau) | fixed cost: base | visits (follow-up fixed) |  73 | 16 |  62 | 1.5 |  924 | 0.40 | 0.44 |
| Aminoglycoside (fitted tau) | fixed cost: base | visits (follow-up fixed) |  53 | 24 |  40 | 1.5 |  924 | 0.32 | 0.37 |
| Aminoglycoside (fitted tau) | fixed cost: base | visits (follow-up fixed) |  42 | 32 |  30 | 1.5 |  924 | 0.26 | 0.33 |
| Aminoglycoside (fitted tau) | fixed cost: base | visits (follow-up fixed) |  29 | 48 |  20 | 1.5 |  924 | 0.21 | 0.24 |
| Aminoglycoside (fitted tau) | fixed cost: base | follow-up (interval fixed) | 164 |  4 |  84 | 1.5 |  252 | 0.43 | 0.12 |
| Aminoglycoside (fitted tau) | fixed cost: base | follow-up (interval fixed) | 136 |  6 |  84 | 1.5 |  420 | 0.45 | 0.20 |
| Aminoglycoside (fitted tau) | fixed cost: base | follow-up (interval fixed) | 116 |  8 |  84 | 1.5 |  588 | 0.44 | 0.29 |
| Aminoglycoside (fitted tau) | fixed cost: base | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.44 | 0.41 |
| Aminoglycoside (fitted tau) | fixed cost: base | follow-up (interval fixed) |  73 | 16 |  84 | 1.5 | 1260 | 0.42 | 0.49 |
| Aminoglycoside (fitted tau) | fixed cost: base | follow-up (interval fixed) |  53 | 24 |  84 | 1.5 | 1932 | 0.43 | 0.56 |
| Aminoglycoside (fitted tau) | fixed cost: base | follow-up (interval fixed) |  42 | 32 |  84 | 1.5 | 2604 | 0.43 | 0.63 |
| Aminoglycoside (fitted tau) | fixed cost: base | follow-up (interval fixed) |  29 | 48 |  84 | 1.5 | 3948 | 0.43 | 0.61 |
| Aminoglycoside (fitted tau) | fixed cost: base | depth |  99 | 12 |  84 | 1.0 |  924 | 0.45 | 0.40 |
| Aminoglycoside (fitted tau) | fixed cost: base | depth |  90 | 12 |  84 | 1.5 |  924 | 0.42 | 0.40 |
| Aminoglycoside (fitted tau) | fixed cost: base | depth |  82 | 12 |  84 | 2.0 |  924 | 0.44 | 0.45 |
| Aminoglycoside (fitted tau) | fixed cost: base | depth |  69 | 12 |  84 | 3.0 |  924 | 0.39 | 0.47 |
| Aminoglycoside (fitted tau) | fixed cost: base | depth |  60 | 12 |  84 | 4.0 |  924 | 0.35 | 0.44 |
| Aminoglycoside (fitted tau) | fixed cost: base | depth |  48 | 12 |  84 | 6.0 |  924 | 0.31 | 0.39 |
| Aminoglycoside (fitted tau) | fixed cost: base | depth |  39 | 12 |  84 | 8.0 |  924 | 0.26 | 0.33 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.45 | 0.07 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.40 | 0.14 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.39 | 0.23 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.40 | 0.32 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.44 | 0.42 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.49 | 0.44 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.54 | 0.41 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | visits (follow-up fixed) | 215 |  4 | 308 | 1.5 |  924 | 0.54 | 0.10 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | visits (follow-up fixed) | 159 |  6 | 185 | 1.5 |  924 | 0.51 | 0.25 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | visits (follow-up fixed) | 126 |  8 | 132 | 1.5 |  924 | 0.51 | 0.36 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.45 | 0.41 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | visits (follow-up fixed) |  69 | 16 |  62 | 1.5 |  924 | 0.39 | 0.40 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | visits (follow-up fixed) |  48 | 24 |  40 | 1.5 |  924 | 0.30 | 0.32 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | visits (follow-up fixed) |  36 | 32 |  30 | 1.5 |  924 | 0.21 | 0.30 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | visits (follow-up fixed) |  24 | 48 |  20 | 1.5 |  924 | 0.17 | 0.21 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | follow-up (interval fixed) | 215 |  4 |  84 | 1.5 |  252 | 0.56 | 0.11 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | follow-up (interval fixed) | 159 |  6 |  84 | 1.5 |  420 | 0.51 | 0.23 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | follow-up (interval fixed) | 126 |  8 |  84 | 1.5 |  588 | 0.49 | 0.29 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.45 | 0.42 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | follow-up (interval fixed) |  69 | 16 |  84 | 1.5 | 1260 | 0.41 | 0.46 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | follow-up (interval fixed) |  48 | 24 |  84 | 1.5 | 1932 | 0.37 | 0.54 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | follow-up (interval fixed) |  36 | 32 |  84 | 1.5 | 2604 | 0.39 | 0.54 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | follow-up (interval fixed) |  24 | 48 |  84 | 1.5 | 3948 | 0.34 | 0.52 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | depth |  93 | 12 |  84 | 1.0 |  924 | 0.43 | 0.37 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | depth |  90 | 12 |  84 | 1.5 |  924 | 0.45 | 0.42 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | depth |  86 | 12 |  84 | 2.0 |  924 | 0.46 | 0.47 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | depth |  80 | 12 |  84 | 3.0 |  924 | 0.45 | 0.51 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | depth |  75 | 12 |  84 | 4.0 |  924 | 0.44 | 0.52 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | depth |  67 | 12 |  84 | 6.0 |  924 | 0.42 | 0.52 |
| Aminoglycoside (fitted tau) | fixed cost: visits dear | depth |  60 | 12 |  84 | 8.0 |  924 | 0.37 | 0.50 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.38 | 0.08 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.41 | 0.12 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.46 | 0.15 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.48 | 0.16 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.55 | 0.12 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.59 | 0.10 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.60 | 0.07 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | visits (follow-up fixed) |  90 |  4 | 308 | 1.5 |  924 | 0.28 | 0.06 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | visits (follow-up fixed) |  90 |  6 | 185 | 1.5 |  924 | 0.37 | 0.05 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | visits (follow-up fixed) |  90 |  8 | 132 | 1.5 |  924 | 0.42 | 0.07 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.56 | 0.13 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | visits (follow-up fixed) |  90 | 16 |  62 | 1.5 |  924 | 0.60 | 0.20 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | visits (follow-up fixed) |  90 | 24 |  40 | 1.5 |  924 | 0.68 | 0.28 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | visits (follow-up fixed) |  90 | 32 |  30 | 1.5 |  924 | 0.69 | 0.35 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | visits (follow-up fixed) |  90 | 48 |  20 | 1.5 |  924 | 0.74 | 0.39 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | follow-up (interval fixed) |  90 |  4 |  84 | 1.5 |  252 | 0.25 | 0.07 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | follow-up (interval fixed) |  90 |  6 |  84 | 1.5 |  420 | 0.33 | 0.08 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | follow-up (interval fixed) |  90 |  8 |  84 | 1.5 |  588 | 0.40 | 0.09 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.54 | 0.14 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | follow-up (interval fixed) |  90 | 16 |  84 | 1.5 | 1260 | 0.65 | 0.16 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | follow-up (interval fixed) |  90 | 24 |  84 | 1.5 | 1932 | 0.81 | 0.25 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | follow-up (interval fixed) |  90 | 32 |  84 | 1.5 | 2604 | 0.92 | 0.32 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | follow-up (interval fixed) |  90 | 48 |  84 | 1.5 | 3948 | 0.97 | 0.44 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | depth |  90 | 12 |  84 | 1.0 |  924 | 0.52 | 0.12 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | depth |  90 | 12 |  84 | 1.5 |  924 | 0.55 | 0.14 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | depth |  90 | 12 |  84 | 2.0 |  924 | 0.61 | 0.15 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | depth |  90 | 12 |  84 | 3.0 |  924 | 0.62 | 0.17 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | depth |  90 | 12 |  84 | 4.0 |  924 | 0.64 | 0.18 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | depth |  90 | 12 |  84 | 6.0 |  924 | 0.66 | 0.20 |
| Aminoglycoside, tau = 60 d | N fixed at 90 | depth |  90 | 12 |  84 | 8.0 |  924 | 0.68 | 0.21 |
| Aminoglycoside, tau = 60 d | fixed cost: base | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.41 | 0.07 |
| Aminoglycoside, tau = 60 d | fixed cost: base | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.43 | 0.11 |
| Aminoglycoside, tau = 60 d | fixed cost: base | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.46 | 0.14 |
| Aminoglycoside, tau = 60 d | fixed cost: base | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.50 | 0.14 |
| Aminoglycoside, tau = 60 d | fixed cost: base | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.55 | 0.12 |
| Aminoglycoside, tau = 60 d | fixed cost: base | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.59 | 0.10 |
| Aminoglycoside, tau = 60 d | fixed cost: base | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.61 | 0.08 |
| Aminoglycoside, tau = 60 d | fixed cost: base | visits (follow-up fixed) | 164 |  4 | 308 | 1.5 |  924 | 0.43 | 0.06 |
| Aminoglycoside, tau = 60 d | fixed cost: base | visits (follow-up fixed) | 136 |  6 | 185 | 1.5 |  924 | 0.51 | 0.06 |
| Aminoglycoside, tau = 60 d | fixed cost: base | visits (follow-up fixed) | 116 |  8 | 132 | 1.5 |  924 | 0.52 | 0.09 |
| Aminoglycoside, tau = 60 d | fixed cost: base | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.54 | 0.13 |
| Aminoglycoside, tau = 60 d | fixed cost: base | visits (follow-up fixed) |  73 | 16 |  62 | 1.5 |  924 | 0.52 | 0.16 |
| Aminoglycoside, tau = 60 d | fixed cost: base | visits (follow-up fixed) |  53 | 24 |  40 | 1.5 |  924 | 0.44 | 0.19 |
| Aminoglycoside, tau = 60 d | fixed cost: base | visits (follow-up fixed) |  42 | 32 |  30 | 1.5 |  924 | 0.40 | 0.18 |
| Aminoglycoside, tau = 60 d | fixed cost: base | visits (follow-up fixed) |  29 | 48 |  20 | 1.5 |  924 | 0.28 | 0.15 |
| Aminoglycoside, tau = 60 d | fixed cost: base | follow-up (interval fixed) | 164 |  4 |  84 | 1.5 |  252 | 0.40 | 0.07 |
| Aminoglycoside, tau = 60 d | fixed cost: base | follow-up (interval fixed) | 136 |  6 |  84 | 1.5 |  420 | 0.48 | 0.08 |
| Aminoglycoside, tau = 60 d | fixed cost: base | follow-up (interval fixed) | 116 |  8 |  84 | 1.5 |  588 | 0.51 | 0.10 |
| Aminoglycoside, tau = 60 d | fixed cost: base | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.55 | 0.13 |
| Aminoglycoside, tau = 60 d | fixed cost: base | follow-up (interval fixed) |  73 | 16 |  84 | 1.5 | 1260 | 0.56 | 0.15 |
| Aminoglycoside, tau = 60 d | fixed cost: base | follow-up (interval fixed) |  53 | 24 |  84 | 1.5 | 1932 | 0.58 | 0.16 |
| Aminoglycoside, tau = 60 d | fixed cost: base | follow-up (interval fixed) |  42 | 32 |  84 | 1.5 | 2604 | 0.59 | 0.17 |
| Aminoglycoside, tau = 60 d | fixed cost: base | follow-up (interval fixed) |  29 | 48 |  84 | 1.5 | 3948 | 0.59 | 0.19 |
| Aminoglycoside, tau = 60 d | fixed cost: base | depth |  99 | 12 |  84 | 1.0 |  924 | 0.53 | 0.14 |
| Aminoglycoside, tau = 60 d | fixed cost: base | depth |  90 | 12 |  84 | 1.5 |  924 | 0.51 | 0.15 |
| Aminoglycoside, tau = 60 d | fixed cost: base | depth |  82 | 12 |  84 | 2.0 |  924 | 0.55 | 0.15 |
| Aminoglycoside, tau = 60 d | fixed cost: base | depth |  69 | 12 |  84 | 3.0 |  924 | 0.51 | 0.15 |
| Aminoglycoside, tau = 60 d | fixed cost: base | depth |  60 | 12 |  84 | 4.0 |  924 | 0.47 | 0.13 |
| Aminoglycoside, tau = 60 d | fixed cost: base | depth |  48 | 12 |  84 | 6.0 |  924 | 0.39 | 0.14 |
| Aminoglycoside, tau = 60 d | fixed cost: base | depth |  39 | 12 |  84 | 8.0 |  924 | 0.33 | 0.12 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.67 | 0.06 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.51 | 0.08 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.47 | 0.15 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.46 | 0.27 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.40 | 0.55 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.40 | 0.85 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.40 | 0.94 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | visits (follow-up fixed) |  90 |  4 | 308 | 1.5 |  924 | 0.25 | 0.31 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | visits (follow-up fixed) |  90 |  6 | 185 | 1.5 |  924 | 0.32 | 0.46 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | visits (follow-up fixed) |  90 |  8 | 132 | 1.5 |  924 | 0.36 | 0.49 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.42 | 0.54 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | visits (follow-up fixed) |  90 | 16 |  62 | 1.5 |  924 | 0.41 | 0.57 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | visits (follow-up fixed) |  90 | 24 |  40 | 1.5 |  924 | 0.47 | 0.58 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | visits (follow-up fixed) |  90 | 32 |  30 | 1.5 |  924 | 0.49 | 0.59 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | visits (follow-up fixed) |  90 | 48 |  20 | 1.5 |  924 | 0.49 | 0.59 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | follow-up (interval fixed) |  90 |  4 |  84 | 1.5 |  252 | 0.33 | 0.05 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | follow-up (interval fixed) |  90 |  6 |  84 | 1.5 |  420 | 0.38 | 0.12 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | follow-up (interval fixed) |  90 |  8 |  84 | 1.5 |  588 | 0.39 | 0.24 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.40 | 0.54 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | follow-up (interval fixed) |  90 | 16 |  84 | 1.5 | 1260 | 0.43 | 0.82 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | follow-up (interval fixed) |  90 | 24 |  84 | 1.5 | 1932 | 0.46 | 0.98 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | follow-up (interval fixed) |  90 | 32 |  84 | 1.5 | 2604 | 0.50 | 1.00 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | follow-up (interval fixed) |  90 | 48 |  84 | 1.5 | 3948 | 0.59 | 1.00 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | depth |  90 | 12 |  84 | 1.0 |  924 | 0.37 | 0.43 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | depth |  90 | 12 |  84 | 1.5 |  924 | 0.41 | 0.57 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | depth |  90 | 12 |  84 | 2.0 |  924 | 0.44 | 0.66 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | depth |  90 | 12 |  84 | 3.0 |  924 | 0.45 | 0.75 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | depth |  90 | 12 |  84 | 4.0 |  924 | 0.45 | 0.81 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | depth |  90 | 12 |  84 | 6.0 |  924 | 0.47 | 0.85 |
| Aminoglycoside, tau = 480 d | N fixed at 90 | depth |  90 | 12 |  84 | 8.0 |  924 | 0.48 | 0.86 |
| Aminoglycoside, tau = 480 d | fixed cost: base | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.66 | 0.05 |
| Aminoglycoside, tau = 480 d | fixed cost: base | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.53 | 0.08 |
| Aminoglycoside, tau = 480 d | fixed cost: base | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.46 | 0.15 |
| Aminoglycoside, tau = 480 d | fixed cost: base | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.45 | 0.27 |
| Aminoglycoside, tau = 480 d | fixed cost: base | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.42 | 0.53 |
| Aminoglycoside, tau = 480 d | fixed cost: base | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.40 | 0.83 |
| Aminoglycoside, tau = 480 d | fixed cost: base | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.40 | 0.94 |
| Aminoglycoside, tau = 480 d | fixed cost: base | visits (follow-up fixed) | 164 |  4 | 308 | 1.5 |  924 | 0.43 | 0.52 |
| Aminoglycoside, tau = 480 d | fixed cost: base | visits (follow-up fixed) | 136 |  6 | 185 | 1.5 |  924 | 0.45 | 0.61 |
| Aminoglycoside, tau = 480 d | fixed cost: base | visits (follow-up fixed) | 116 |  8 | 132 | 1.5 |  924 | 0.44 | 0.62 |
| Aminoglycoside, tau = 480 d | fixed cost: base | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.40 | 0.56 |
| Aminoglycoside, tau = 480 d | fixed cost: base | visits (follow-up fixed) |  73 | 16 |  62 | 1.5 |  924 | 0.38 | 0.48 |
| Aminoglycoside, tau = 480 d | fixed cost: base | visits (follow-up fixed) |  53 | 24 |  40 | 1.5 |  924 | 0.30 | 0.37 |
| Aminoglycoside, tau = 480 d | fixed cost: base | visits (follow-up fixed) |  42 | 32 |  30 | 1.5 |  924 | 0.24 | 0.29 |
| Aminoglycoside, tau = 480 d | fixed cost: base | visits (follow-up fixed) |  29 | 48 |  20 | 1.5 |  924 | 0.19 | 0.21 |
| Aminoglycoside, tau = 480 d | fixed cost: base | follow-up (interval fixed) | 164 |  4 |  84 | 1.5 |  252 | 0.58 | 0.08 |
| Aminoglycoside, tau = 480 d | fixed cost: base | follow-up (interval fixed) | 136 |  6 |  84 | 1.5 |  420 | 0.53 | 0.17 |
| Aminoglycoside, tau = 480 d | fixed cost: base | follow-up (interval fixed) | 116 |  8 |  84 | 1.5 |  588 | 0.48 | 0.30 |
| Aminoglycoside, tau = 480 d | fixed cost: base | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.40 | 0.57 |
| Aminoglycoside, tau = 480 d | fixed cost: base | follow-up (interval fixed) |  73 | 16 |  84 | 1.5 | 1260 | 0.35 | 0.71 |
| Aminoglycoside, tau = 480 d | fixed cost: base | follow-up (interval fixed) |  53 | 24 |  84 | 1.5 | 1932 | 0.28 | 0.88 |
| Aminoglycoside, tau = 480 d | fixed cost: base | follow-up (interval fixed) |  42 | 32 |  84 | 1.5 | 2604 | 0.24 | 0.92 |
| Aminoglycoside, tau = 480 d | fixed cost: base | follow-up (interval fixed) |  29 | 48 |  84 | 1.5 | 3948 | 0.21 | 0.94 |
| Aminoglycoside, tau = 480 d | fixed cost: base | depth |  99 | 12 |  84 | 1.0 |  924 | 0.41 | 0.51 |
| Aminoglycoside, tau = 480 d | fixed cost: base | depth |  90 | 12 |  84 | 1.5 |  924 | 0.41 | 0.53 |
| Aminoglycoside, tau = 480 d | fixed cost: base | depth |  82 | 12 |  84 | 2.0 |  924 | 0.40 | 0.62 |
| Aminoglycoside, tau = 480 d | fixed cost: base | depth |  69 | 12 |  84 | 3.0 |  924 | 0.35 | 0.64 |
| Aminoglycoside, tau = 480 d | fixed cost: base | depth |  60 | 12 |  84 | 4.0 |  924 | 0.32 | 0.64 |
| Aminoglycoside, tau = 480 d | fixed cost: base | depth |  48 | 12 |  84 | 6.0 |  924 | 0.25 | 0.55 |
| Aminoglycoside, tau = 480 d | fixed cost: base | depth |  39 | 12 |  84 | 8.0 |  924 | 0.23 | 0.47 |
| Cephalosporin (fitted) | N fixed at 90 | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.98 | 0.07 |
| Cephalosporin (fitted) | N fixed at 90 | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.96 | 0.14 |
| Cephalosporin (fitted) | N fixed at 90 | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.95 | 0.27 |
| Cephalosporin (fitted) | N fixed at 90 | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.94 | 0.43 |
| Cephalosporin (fitted) | N fixed at 90 | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.95 | 0.66 |
| Cephalosporin (fitted) | N fixed at 90 | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.96 | 0.82 |
| Cephalosporin (fitted) | N fixed at 90 | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.97 | 0.87 |
| Cephalosporin (fitted) | N fixed at 90 | visits (follow-up fixed) |  90 |  4 | 308 | 1.5 |  924 | 0.71 | 0.18 |
| Cephalosporin (fitted) | N fixed at 90 | visits (follow-up fixed) |  90 |  6 | 185 | 1.5 |  924 | 0.83 | 0.39 |
| Cephalosporin (fitted) | N fixed at 90 | visits (follow-up fixed) |  90 |  8 | 132 | 1.5 |  924 | 0.89 | 0.54 |
| Cephalosporin (fitted) | N fixed at 90 | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.94 | 0.68 |
| Cephalosporin (fitted) | N fixed at 90 | visits (follow-up fixed) |  90 | 16 |  62 | 1.5 |  924 | 0.97 | 0.74 |
| Cephalosporin (fitted) | N fixed at 90 | visits (follow-up fixed) |  90 | 24 |  40 | 1.5 |  924 | 0.98 | 0.81 |
| Cephalosporin (fitted) | N fixed at 90 | visits (follow-up fixed) |  90 | 32 |  30 | 1.5 |  924 | 0.98 | 0.83 |
| Cephalosporin (fitted) | N fixed at 90 | visits (follow-up fixed) |  90 | 48 |  20 | 1.5 |  924 | 0.98 | 0.88 |
| Cephalosporin (fitted) | N fixed at 90 | follow-up (interval fixed) |  90 |  4 |  84 | 1.5 |  252 | 0.77 | 0.09 |
| Cephalosporin (fitted) | N fixed at 90 | follow-up (interval fixed) |  90 |  6 |  84 | 1.5 |  420 | 0.85 | 0.19 |
| Cephalosporin (fitted) | N fixed at 90 | follow-up (interval fixed) |  90 |  8 |  84 | 1.5 |  588 | 0.88 | 0.35 |
| Cephalosporin (fitted) | N fixed at 90 | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.95 | 0.70 |
| Cephalosporin (fitted) | N fixed at 90 | follow-up (interval fixed) |  90 | 16 |  84 | 1.5 | 1260 | 0.98 | 0.87 |
| Cephalosporin (fitted) | N fixed at 90 | follow-up (interval fixed) |  90 | 24 |  84 | 1.5 | 1932 | 0.99 | 0.99 |
| Cephalosporin (fitted) | N fixed at 90 | follow-up (interval fixed) |  90 | 32 |  84 | 1.5 | 2604 | 1.00 | 1.00 |
| Cephalosporin (fitted) | N fixed at 90 | follow-up (interval fixed) |  90 | 48 |  84 | 1.5 | 3948 | 1.00 | 1.00 |
| Cephalosporin (fitted) | N fixed at 90 | depth |  90 | 12 |  84 | 1.0 |  924 | 0.92 | 0.54 |
| Cephalosporin (fitted) | N fixed at 90 | depth |  90 | 12 |  84 | 1.5 |  924 | 0.94 | 0.67 |
| Cephalosporin (fitted) | N fixed at 90 | depth |  90 | 12 |  84 | 2.0 |  924 | 0.97 | 0.83 |
| Cephalosporin (fitted) | N fixed at 90 | depth |  90 | 12 |  84 | 3.0 |  924 | 0.98 | 0.93 |
| Cephalosporin (fitted) | N fixed at 90 | depth |  90 | 12 |  84 | 4.0 |  924 | 0.99 | 0.96 |
| Cephalosporin (fitted) | N fixed at 90 | depth |  90 | 12 |  84 | 6.0 |  924 | 0.99 | 0.99 |
| Cephalosporin (fitted) | N fixed at 90 | depth |  90 | 12 |  84 | 8.0 |  924 | 0.99 | 0.99 |
| Cephalosporin (fitted) | fixed cost: base | interval (L fixed) |  90 | 12 |  14 | 1.5 |  154 | 0.99 | 0.06 |
| Cephalosporin (fitted) | fixed cost: base | interval (L fixed) |  90 | 12 |  28 | 1.5 |  308 | 0.96 | 0.15 |
| Cephalosporin (fitted) | fixed cost: base | interval (L fixed) |  90 | 12 |  42 | 1.5 |  462 | 0.95 | 0.26 |
| Cephalosporin (fitted) | fixed cost: base | interval (L fixed) |  90 | 12 |  56 | 1.5 |  616 | 0.95 | 0.43 |
| Cephalosporin (fitted) | fixed cost: base | interval (L fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.94 | 0.67 |
| Cephalosporin (fitted) | fixed cost: base | interval (L fixed) |  90 | 12 | 126 | 1.5 | 1386 | 0.95 | 0.80 |
| Cephalosporin (fitted) | fixed cost: base | interval (L fixed) |  90 | 12 | 168 | 1.5 | 1848 | 0.98 | 0.86 |
| Cephalosporin (fitted) | fixed cost: base | visits (follow-up fixed) | 164 |  4 | 308 | 1.5 |  924 | 0.93 | 0.27 |
| Cephalosporin (fitted) | fixed cost: base | visits (follow-up fixed) | 136 |  6 | 185 | 1.5 |  924 | 0.96 | 0.54 |
| Cephalosporin (fitted) | fixed cost: base | visits (follow-up fixed) | 116 |  8 | 132 | 1.5 |  924 | 0.96 | 0.65 |
| Cephalosporin (fitted) | fixed cost: base | visits (follow-up fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.95 | 0.66 |
| Cephalosporin (fitted) | fixed cost: base | visits (follow-up fixed) |  73 | 16 |  62 | 1.5 |  924 | 0.92 | 0.65 |
| Cephalosporin (fitted) | fixed cost: base | visits (follow-up fixed) |  53 | 24 |  40 | 1.5 |  924 | 0.86 | 0.55 |
| Cephalosporin (fitted) | fixed cost: base | visits (follow-up fixed) |  42 | 32 |  30 | 1.5 |  924 | 0.81 | 0.50 |
| Cephalosporin (fitted) | fixed cost: base | visits (follow-up fixed) |  29 | 48 |  20 | 1.5 |  924 | 0.67 | 0.41 |
| Cephalosporin (fitted) | fixed cost: base | follow-up (interval fixed) | 164 |  4 |  84 | 1.5 |  252 | 0.95 | 0.10 |
| Cephalosporin (fitted) | fixed cost: base | follow-up (interval fixed) | 136 |  6 |  84 | 1.5 |  420 | 0.96 | 0.24 |
| Cephalosporin (fitted) | fixed cost: base | follow-up (interval fixed) | 116 |  8 |  84 | 1.5 |  588 | 0.96 | 0.40 |
| Cephalosporin (fitted) | fixed cost: base | follow-up (interval fixed) |  90 | 12 |  84 | 1.5 |  924 | 0.95 | 0.67 |
| Cephalosporin (fitted) | fixed cost: base | follow-up (interval fixed) |  73 | 16 |  84 | 1.5 | 1260 | 0.94 | 0.79 |
| Cephalosporin (fitted) | fixed cost: base | follow-up (interval fixed) |  53 | 24 |  84 | 1.5 | 1932 | 0.93 | 0.91 |
| Cephalosporin (fitted) | fixed cost: base | follow-up (interval fixed) |  42 | 32 |  84 | 1.5 | 2604 | 0.93 | 0.93 |
| Cephalosporin (fitted) | fixed cost: base | follow-up (interval fixed) |  29 | 48 |  84 | 1.5 | 3948 | 0.90 | 0.95 |
| Cephalosporin (fitted) | fixed cost: base | depth |  99 | 12 |  84 | 1.0 |  924 | 0.94 | 0.59 |
| Cephalosporin (fitted) | fixed cost: base | depth |  90 | 12 |  84 | 1.5 |  924 | 0.95 | 0.64 |
| Cephalosporin (fitted) | fixed cost: base | depth |  82 | 12 |  84 | 2.0 |  924 | 0.95 | 0.80 |
| Cephalosporin (fitted) | fixed cost: base | depth |  69 | 12 |  84 | 3.0 |  924 | 0.94 | 0.85 |
| Cephalosporin (fitted) | fixed cost: base | depth |  60 | 12 |  84 | 4.0 |  924 | 0.92 | 0.85 |
| Cephalosporin (fitted) | fixed cost: base | depth |  48 | 12 |  84 | 6.0 |  924 | 0.87 | 0.86 |
| Cephalosporin (fitted) | fixed cost: base | depth |  39 | 12 |  84 | 8.0 |  924 | 0.81 | 0.81 |
