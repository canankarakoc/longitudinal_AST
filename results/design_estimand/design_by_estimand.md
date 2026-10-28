# Sampling design by estimand

D = 20000 simulated visits per cell; plate of 50 colonies for morphotype visibility.
Pooled / whole-sputum readouts use ASSUMED observation rules (NS if NS share >= 1% or >= 50%).

## Cohort anchor: morphotype alignment inside cultures (colony frame, patient bootstrap)

| class | cultures_between | cultures_within | P(differ), different morphotype | P(differ), same morphotype | alignment |
| --- | --- | --- | --- | --- | --- |
| Aminoglycoside | 286 | 221 | 0.24 (0.16–0.34) | 0.17 (0.11–0.24) | 0.31 (-0.04–0.54) |
| Carbapenem | 282 | 217 | 0.30 (0.21–0.38) | 0.21 (0.14–0.29) | 0.29 (-0.08–0.54) |
| Cephalosporin | 286 | 221 | 0.40 (0.31–0.48) | 0.22 (0.15–0.30) | 0.45 (0.23–0.62) |
| Fluoroquinolone | 284 | 220 | 0.22 (0.15–0.29) | 0.15 (0.09–0.22) | 0.32 (0.04–0.55) |

alignment = 1 - within/between: 1 = resistance follows morphotype exactly, 0 = morphotype uninformative.

## Snapshot: theta = 0.025 (rare minority), no assay error

| structure | design | cost | P(report NS) | freq MAE | naive freq MAE | dominant correct | false change (detect) | false change (dominant) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aligned | random-1 | 1.0 | 3% | 0.05 | – | 97% | 5% | 5% |
| aligned | random-4 | 4.0 | 10% | 0.05 | – | 100% | 18% | 0% |
| aligned | random-8 | 8.0 | 18% | 0.04 | – | 100% | 30% | 0% |
| aligned | stratified-1 | 1.7 | 72% | 0.01 | 0.35 | 100% | 40% | 0% |
| aligned | stratified-2 | 3.4 | 72% | 0.01 | 0.35 | 100% | 40% | 0% |
| aligned | pooled-10, NS if >= 1% | 1.0 | 22% | – | – | 78% | 35% | 35% |
| aligned | pooled-10, NS if >= 50% | 1.0 | 0% | – | – | 100% | 0% | 0% |
| aligned | whole sputum, NS if >= 1% | 1.0 | 100% | – | – | 0% | 0% | 0% |
| aligned | whole sputum, NS if >= 50% | 1.0 | 0% | – | – | 100% | 0% | 0% |
| unaligned | random-1 | 1.0 | 3% | 0.05 | – | 97% | 5% | 5% |
| unaligned | random-4 | 4.0 | 10% | 0.05 | – | 100% | 18% | 0% |
| unaligned | random-8 | 8.0 | 18% | 0.04 | – | 100% | 30% | 0% |
| unaligned | stratified-1 | 2.0 | 5% | 0.05 | 0.05 | 98% | 9% | 5% |
| unaligned | stratified-2 | 4.0 | 9% | 0.04 | 0.04 | 100% | 17% | 0% |
| unaligned | pooled-10, NS if >= 1% | 1.0 | 23% | – | – | 77% | 35% | 35% |
| unaligned | pooled-10, NS if >= 50% | 1.0 | 0% | – | – | 100% | 0% | 0% |
| unaligned | whole sputum, NS if >= 1% | 1.0 | 100% | – | – | 0% | 0% | 0% |
| unaligned | whole sputum, NS if >= 50% | 1.0 | 0% | – | – | 100% | 0% | 0% |

## Snapshot: theta = 0.225, no assay error

| structure | design | cost | P(report NS) | freq MAE | naive freq MAE | dominant correct | false change (detect) | false change (dominant) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aligned | random-1 | 1.0 | 23% | 0.35 | – | 77% | 35% | 35% |
| aligned | random-4 | 4.0 | 65% | 0.16 | – | 87% | 46% | 23% |
| aligned | random-8 | 8.0 | 87% | 0.12 | – | 95% | 23% | 10% |
| aligned | stratified-1 | 2.0 | 100% | 0.00 | 0.28 | 100% | 0% | 0% |
| aligned | stratified-2 | 4.0 | 100% | 0.00 | 0.28 | 100% | 0% | 0% |
| aligned | pooled-10, NS if >= 1% | 1.0 | 92% | – | – | 8% | 15% | 15% |
| aligned | pooled-10, NS if >= 50% | 1.0 | 5% | – | – | 95% | 10% | 10% |
| aligned | whole sputum, NS if >= 1% | 1.0 | 100% | – | – | 0% | 0% | 0% |
| aligned | whole sputum, NS if >= 50% | 1.0 | 0% | – | – | 100% | 0% | 0% |
| unaligned | random-1 | 1.0 | 23% | 0.35 | – | 77% | 35% | 35% |
| unaligned | random-4 | 4.0 | 64% | 0.16 | – | 88% | 46% | 22% |
| unaligned | random-8 | 8.0 | 87% | 0.12 | – | 95% | 23% | 9% |
| unaligned | stratified-1 | 2.0 | 40% | 0.27 | 0.27 | 78% | 48% | 35% |
| unaligned | stratified-2 | 4.0 | 64% | 0.19 | 0.16 | 87% | 46% | 23% |
| unaligned | pooled-10, NS if >= 1% | 1.0 | 92% | – | – | 8% | 14% | 14% |
| unaligned | pooled-10, NS if >= 50% | 1.0 | 6% | – | – | 94% | 10% | 10% |
| unaligned | whole sputum, NS if >= 1% | 1.0 | 100% | – | – | 0% | 0% | 0% |
| unaligned | whole sputum, NS if >= 50% | 1.0 | 0% | – | – | 100% | 0% | 0% |

## Snapshot: theta = 0.225, 3% per-test error

| structure | design | cost | P(report NS) | freq MAE | naive freq MAE | dominant correct | false change (detect) | false change (dominant) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aligned | random-1 | 1.0 | 24% | 0.36 | – | 76% | 37% | 37% |
| aligned | random-4 | 4.0 | 67% | 0.16 | – | 85% | 44% | 25% |
| aligned | random-8 | 8.0 | 89% | 0.12 | – | 94% | 20% | 12% |
| aligned | stratified-1 | 2.0 | 97% | 0.03 | 0.29 | 97% | 6% | 6% |
| aligned | stratified-2 | 4.0 | 100% | 0.03 | 0.28 | 94% | 0% | 11% |
| aligned | pooled-10, NS if >= 1% | 1.0 | 90% | – | – | 10% | 19% | 19% |
| aligned | pooled-10, NS if >= 50% | 1.0 | 8% | – | – | 92% | 15% | 15% |
| aligned | whole sputum, NS if >= 1% | 1.0 | 97% | – | – | 3% | 6% | 6% |
| aligned | whole sputum, NS if >= 50% | 1.0 | 3% | – | – | 97% | 6% | 6% |
| unaligned | random-1 | 1.0 | 25% | 0.36 | – | 75% | 37% | 37% |
| unaligned | random-4 | 4.0 | 67% | 0.17 | – | 85% | 44% | 25% |
| unaligned | random-8 | 8.0 | 89% | 0.12 | – | 94% | 20% | 12% |
| unaligned | stratified-1 | 2.0 | 43% | 0.27 | 0.27 | 76% | 49% | 36% |
| unaligned | stratified-2 | 4.0 | 67% | 0.20 | 0.17 | 85% | 45% | 25% |
| unaligned | pooled-10, NS if >= 1% | 1.0 | 89% | – | – | 11% | 19% | 19% |
| unaligned | pooled-10, NS if >= 50% | 1.0 | 8% | – | – | 92% | 14% | 14% |
| unaligned | whole sputum, NS if >= 1% | 1.0 | 97% | – | – | 3% | 6% | 6% |
| unaligned | whole sputum, NS if >= 50% | 1.0 | 3% | – | – | 97% | 6% | 6% |

## Snapshot: theta = 0 (false positives), 3% per-test error

| structure | design | cost | P(report NS) | freq MAE | naive freq MAE | dominant correct | false change (detect) | false change (dominant) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aligned | random-1 | 1.0 | 3% | 0.03 | – | 97% | 6% | 6% |
| aligned | random-4 | 4.0 | 11% | 0.03 | – | 100% | 20% | 1% |
| aligned | random-8 | 8.0 | 22% | 0.03 | – | 100% | 34% | 0% |
| aligned | stratified-1 | 1.0 | 3% | 0.03 | 0.03 | 97% | 6% | 6% |
| aligned | stratified-2 | 2.0 | 6% | 0.03 | 0.03 | 97% | 11% | 6% |
| aligned | pooled-10, NS if >= 1% | 1.0 | 3% | – | – | 97% | 6% | 6% |
| aligned | pooled-10, NS if >= 50% | 1.0 | 3% | – | – | 97% | 6% | 6% |
| aligned | whole sputum, NS if >= 1% | 1.0 | 3% | – | – | 97% | 6% | 6% |
| aligned | whole sputum, NS if >= 50% | 1.0 | 3% | – | – | 97% | 6% | 6% |
| unaligned | random-1 | 1.0 | 3% | 0.03 | – | 97% | 6% | 6% |
| unaligned | random-4 | 4.0 | 11% | 0.03 | – | 100% | 20% | 1% |
| unaligned | random-8 | 8.0 | 22% | 0.03 | – | 100% | 34% | 0% |
| unaligned | stratified-1 | 2.0 | 6% | 0.03 | 0.03 | 97% | 11% | 6% |
| unaligned | stratified-2 | 4.0 | 12% | 0.03 | 0.03 | 100% | 21% | 1% |
| unaligned | pooled-10, NS if >= 1% | 1.0 | 3% | – | – | 97% | 6% | 6% |
| unaligned | pooled-10, NS if >= 50% | 1.0 | 3% | – | – | 97% | 6% | 6% |
| unaligned | whole sputum, NS if >= 1% | 1.0 | 3% | – | – | 97% | 6% | 6% |
| unaligned | whole sputum, NS if >= 50% | 1.0 | 3% | – | – | 97% | 6% | 6% |

## True change vs false change (no assay error)

| structure | change | design | P(differ), detect call | false, detect | P(differ), dominant call | false, dominant | P(|d freq| >= 0.25) | false, freq |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| aligned | 0.000 -> 0.225 | random-1 | 22% | 0% | 22% | 0% | 22% | 0% |
| aligned | 0.000 -> 0.225 | random-4 | 64% | 0% | 13% | 0% | 64% | 0% |
| aligned | 0.000 -> 0.225 | random-8 | 87% | 0% | 5% | 0% | 57% | 0% |
| aligned | 0.000 -> 0.225 | stratified-1 | 100% | 0% | 0% | 0% | 0% | 0% |
| aligned | 0.000 -> 0.225 | stratified-2 | 100% | 0% | 0% | 0% | 0% | 0% |
| aligned | 0.000 -> 0.225 | pooled-10, NS if >= 1% | 92% | 0% | 92% | 0% | – | – |
| aligned | 0.000 -> 0.225 | pooled-10, NS if >= 50% | 5% | 0% | 5% | 0% | – | – |
| aligned | 0.000 -> 0.225 | whole sputum, NS if >= 1% | 100% | 0% | 100% | 0% | – | – |
| aligned | 0.000 -> 0.225 | whole sputum, NS if >= 50% | 0% | 0% | 0% | 0% | – | – |
| aligned | 0.075 -> 0.375 | random-1 | 40% | 14% | 40% | 14% | 40% | 14% |
| aligned | 0.075 -> 0.375 | random-4 | 66% | 40% | 32% | 3% | 79% | 41% |
| aligned | 0.075 -> 0.375 | random-8 | 53% | 49% | 24% | 0% | 73% | 13% |
| aligned | 0.075 -> 0.375 | stratified-1 | 2% | 4% | 0% | 0% | 100% | 0% |
| aligned | 0.075 -> 0.375 | stratified-2 | 2% | 4% | 0% | 0% | 100% | 0% |
| aligned | 0.075 -> 0.375 | pooled-10, NS if >= 1% | 46% | 50% | 46% | 50% | – | – |
| aligned | 0.075 -> 0.375 | pooled-10, NS if >= 50% | 30% | 0% | 30% | 0% | – | – |
| aligned | 0.075 -> 0.375 | whole sputum, NS if >= 1% | 0% | 0% | 0% | 0% | – | – |
| aligned | 0.075 -> 0.375 | whole sputum, NS if >= 50% | 0% | 0% | 0% | 0% | – | – |
| aligned | 0.225 -> 0.725 | random-1 | 62% | 35% | 62% | 35% | 62% | 35% |
| aligned | 0.225 -> 0.725 | random-4 | 36% | 47% | 72% | 23% | 91% | 66% |
| aligned | 0.225 -> 0.725 | random-8 | 13% | 23% | 87% | 9% | 92% | 36% |
| aligned | 0.225 -> 0.725 | stratified-1 | 0% | 0% | 100% | 0% | 100% | 0% |
| aligned | 0.225 -> 0.725 | stratified-2 | 0% | 0% | 100% | 0% | 100% | 0% |
| aligned | 0.225 -> 0.725 | pooled-10, NS if >= 1% | 8% | 15% | 8% | 15% | – | – |
| aligned | 0.225 -> 0.725 | pooled-10, NS if >= 50% | 92% | 10% | 92% | 10% | – | – |
| aligned | 0.225 -> 0.725 | whole sputum, NS if >= 1% | 0% | 0% | 0% | 0% | – | – |
| aligned | 0.225 -> 0.725 | whole sputum, NS if >= 50% | 100% | 0% | 100% | 0% | – | – |
| aligned | 0.725 -> 0.975 | random-1 | 29% | 40% | 29% | 40% | 29% | 40% |
| aligned | 0.725 -> 0.975 | random-4 | 1% | 1% | 19% | 30% | 71% | 69% |
| aligned | 0.725 -> 0.975 | random-8 | 0% | 0% | 10% | 17% | 63% | 40% |
| aligned | 0.725 -> 0.975 | stratified-1 | 0% | 0% | 0% | 0% | 100% | 0% |
| aligned | 0.725 -> 0.975 | stratified-2 | 0% | 0% | 0% | 0% | 100% | 0% |
| aligned | 0.725 -> 0.975 | pooled-10, NS if >= 1% | 0% | 0% | 0% | 0% | – | – |
| aligned | 0.725 -> 0.975 | pooled-10, NS if >= 50% | 3% | 6% | 3% | 6% | – | – |
| aligned | 0.725 -> 0.975 | whole sputum, NS if >= 1% | 0% | 0% | 0% | 0% | – | – |
| aligned | 0.725 -> 0.975 | whole sputum, NS if >= 50% | 0% | 0% | 0% | 0% | – | – |
| unaligned | 0.000 -> 0.225 | random-1 | 23% | 0% | 23% | 0% | 23% | 0% |
| unaligned | 0.000 -> 0.225 | random-4 | 64% | 0% | 13% | 0% | 64% | 0% |
| unaligned | 0.000 -> 0.225 | random-8 | 87% | 0% | 5% | 0% | 56% | 0% |
| unaligned | 0.000 -> 0.225 | stratified-1 | 40% | 0% | 22% | 0% | 40% | 0% |
| unaligned | 0.000 -> 0.225 | stratified-2 | 64% | 0% | 13% | 0% | 43% | 0% |
| unaligned | 0.000 -> 0.225 | pooled-10, NS if >= 1% | 92% | 0% | 92% | 0% | – | – |
| unaligned | 0.000 -> 0.225 | pooled-10, NS if >= 50% | 5% | 0% | 5% | 0% | – | – |
| unaligned | 0.000 -> 0.225 | whole sputum, NS if >= 1% | 100% | 0% | 100% | 0% | – | – |
| unaligned | 0.000 -> 0.225 | whole sputum, NS if >= 50% | 0% | 0% | 0% | 0% | – | – |
| unaligned | 0.075 -> 0.375 | random-1 | 40% | 14% | 40% | 14% | 40% | 14% |
| unaligned | 0.075 -> 0.375 | random-4 | 66% | 40% | 33% | 3% | 80% | 41% |
| unaligned | 0.075 -> 0.375 | random-8 | 53% | 50% | 24% | 0% | 73% | 14% |
| unaligned | 0.075 -> 0.375 | stratified-1 | 57% | 25% | 39% | 14% | 63% | 26% |
| unaligned | 0.075 -> 0.375 | stratified-2 | 66% | 39% | 32% | 3% | 60% | 22% |
| unaligned | 0.075 -> 0.375 | pooled-10, NS if >= 1% | 46% | 49% | 46% | 49% | – | – |
| unaligned | 0.075 -> 0.375 | pooled-10, NS if >= 50% | 31% | 0% | 31% | 0% | – | – |
| unaligned | 0.075 -> 0.375 | whole sputum, NS if >= 1% | 0% | 0% | 0% | 0% | – | – |
| unaligned | 0.075 -> 0.375 | whole sputum, NS if >= 50% | 0% | 0% | 0% | 0% | – | – |
| unaligned | 0.225 -> 0.725 | random-1 | 62% | 35% | 62% | 35% | 62% | 35% |
| unaligned | 0.225 -> 0.725 | random-4 | 37% | 46% | 73% | 22% | 92% | 66% |
| unaligned | 0.225 -> 0.725 | random-8 | 13% | 23% | 86% | 10% | 92% | 36% |
| unaligned | 0.225 -> 0.725 | stratified-1 | 58% | 48% | 62% | 35% | 86% | 58% |
| unaligned | 0.225 -> 0.725 | stratified-2 | 36% | 46% | 74% | 22% | 80% | 42% |
| unaligned | 0.225 -> 0.725 | pooled-10, NS if >= 1% | 8% | 14% | 8% | 14% | – | – |
| unaligned | 0.225 -> 0.725 | pooled-10, NS if >= 50% | 92% | 10% | 92% | 10% | – | – |
| unaligned | 0.225 -> 0.725 | whole sputum, NS if >= 1% | 0% | 0% | 0% | 0% | – | – |
| unaligned | 0.225 -> 0.725 | whole sputum, NS if >= 50% | 100% | 0% | 100% | 0% | – | – |
| unaligned | 0.725 -> 0.975 | random-1 | 28% | 40% | 28% | 40% | 28% | 40% |
| unaligned | 0.725 -> 0.975 | random-4 | 1% | 1% | 19% | 31% | 71% | 69% |
| unaligned | 0.725 -> 0.975 | random-8 | 0% | 0% | 10% | 18% | 63% | 40% |
| unaligned | 0.725 -> 0.975 | stratified-1 | 7% | 14% | 28% | 39% | 49% | 64% |
| unaligned | 0.725 -> 0.975 | stratified-2 | 1% | 1% | 19% | 31% | 50% | 46% |
| unaligned | 0.725 -> 0.975 | pooled-10, NS if >= 1% | 0% | 0% | 0% | 0% | – | – |
| unaligned | 0.725 -> 0.975 | pooled-10, NS if >= 50% | 3% | 6% | 3% | 6% | – | – |
| unaligned | 0.725 -> 0.975 | whole sputum, NS if >= 1% | 0% | 0% | 0% | 0% | – | – |
| unaligned | 0.725 -> 0.975 | whole sputum, NS if >= 50% | 0% | 0% | 0% | 0% | – | – |

## Empirical: heterogeneous Vanderwoude cells

| sample | drug | NS | design | P(report NS) | freq MAE | dominant correct | false change (detect) | false change (dominant) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 (cohort id 7) | amikacin | 57% | random-1 | 57% | 0.49 | 57% | 50% | 50% |
| P1 (cohort id 7) | amikacin | 57% | random-4 | 97% | 0.20 | 61% | 6% | 47% |
| P1 (cohort id 7) | amikacin | 57% | random-8 | 100% | 0.14 | 67% | 0% | 43% |
| P1 (cohort id 7) | amikacin | 57% | pooled-10, NS if >= 1% | 100% | – | 100% | 0% | 0% |
| P1 (cohort id 7) | amikacin | 57% | pooled-10, NS if >= 50% | 80% | – | 80% | 31% | 31% |
| P1 (cohort id 7) | amikacin | 57% | whole sputum, NS if >= 1% | 100% | – | 100% | 0% | 0% |
| P1 (cohort id 7) | amikacin | 57% | whole sputum, NS if >= 50% | 100% | – | 100% | 0% | 0% |
| P1 (cohort id 7) | ciprofloxacin | 100% | random-1 | 100% | 0.00 | 100% | 0% | 0% |
| P1 (cohort id 7) | ciprofloxacin | 100% | random-4 | 100% | 0.00 | 100% | 0% | 0% |
| P1 (cohort id 7) | ciprofloxacin | 100% | random-8 | 100% | 0.00 | 100% | 0% | 0% |
| P1 (cohort id 7) | ciprofloxacin | 100% | pooled-10, NS if >= 1% | 100% | – | 100% | 0% | 0% |
| P1 (cohort id 7) | ciprofloxacin | 100% | pooled-10, NS if >= 50% | 100% | – | 100% | 0% | 0% |
| P1 (cohort id 7) | ciprofloxacin | 100% | whole sputum, NS if >= 1% | 100% | – | 100% | 0% | 0% |
| P1 (cohort id 7) | ciprofloxacin | 100% | whole sputum, NS if >= 50% | 100% | – | 100% | 0% | 0% |
| P3 (cohort id 38) | ceftazidime | 7% | random-1 | 7% | 0.12 | 93% | 13% | 13% |
| P3 (cohort id 38) | ceftazidime | 7% | random-4 | 23% | 0.10 | 99% | 36% | 3% |
| P3 (cohort id 38) | ceftazidime | 7% | random-8 | 45% | 0.08 | 100% | 51% | 0% |
| P3 (cohort id 38) | ceftazidime | 7% | pooled-10, NS if >= 1% | 52% | – | 48% | 49% | 49% |
| P3 (cohort id 38) | ceftazidime | 7% | pooled-10, NS if >= 50% | 0% | – | 100% | 0% | 0% |
| P3 (cohort id 38) | ceftazidime | 7% | whole sputum, NS if >= 1% | 100% | – | 0% | 0% | 0% |
| P3 (cohort id 38) | ceftazidime | 7% | whole sputum, NS if >= 50% | 0% | – | 100% | 0% | 0% |
| P3 (cohort id 38) | ciprofloxacin | 100% | random-1 | 100% | 0.00 | 100% | 0% | 0% |
| P3 (cohort id 38) | ciprofloxacin | 100% | random-4 | 100% | 0.00 | 100% | 0% | 0% |
| P3 (cohort id 38) | ciprofloxacin | 100% | random-8 | 100% | 0.00 | 100% | 0% | 0% |
| P3 (cohort id 38) | ciprofloxacin | 100% | pooled-10, NS if >= 1% | 100% | – | 100% | 0% | 0% |
| P3 (cohort id 38) | ciprofloxacin | 100% | pooled-10, NS if >= 50% | 100% | – | 100% | 0% | 0% |
| P3 (cohort id 38) | ciprofloxacin | 100% | whole sputum, NS if >= 1% | 100% | – | 100% | 0% | 0% |
| P3 (cohort id 38) | ciprofloxacin | 100% | whole sputum, NS if >= 50% | 100% | – | 100% | 0% | 0% |
| P3 (cohort id 38) | piperacillin-tazobactam | 4% | random-1 | 4% | 0.07 | 96% | 8% | 8% |
| P3 (cohort id 38) | piperacillin-tazobactam | 4% | random-4 | 15% | 0.07 | 100% | 27% | 1% |
| P3 (cohort id 38) | piperacillin-tazobactam | 4% | random-8 | 28% | 0.06 | 100% | 39% | 0% |
| P3 (cohort id 38) | piperacillin-tazobactam | 4% | pooled-10, NS if >= 1% | 36% | – | 64% | 46% | 46% |
| P3 (cohort id 38) | piperacillin-tazobactam | 4% | pooled-10, NS if >= 50% | 0% | – | 100% | 0% | 0% |
| P3 (cohort id 38) | piperacillin-tazobactam | 4% | whole sputum, NS if >= 1% | 100% | – | 0% | 0% | 0% |
| P3 (cohort id 38) | piperacillin-tazobactam | 4% | whole sputum, NS if >= 50% | 0% | – | 100% | 0% | 0% |
