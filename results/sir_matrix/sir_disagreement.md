# Within-culture S/I/R disagreement

Nine drugs (amikacin, gentamicin, tobramycin, meropenem, cefepime, ceftazidime, levofloxacin, piperacillin_tazobactam, aztreonam). Patient-cluster bootstrap B = 1000; 95% percentile intervals.
'culture' weighting: each culture x drug counts once (primary). 'pairs': each ordered colony pair counts once (as in AMR_sir_ordinal.R).

## Matrix, all drugs, same culture, culture-weighted

| first \ other | S | I | R |
| --- | --- | --- | --- |
| S | 73% (66–79) | 8% (7–10) | 19% (14–25) |
| I | 36% (30–42) | 25% (21–30) | 39% (33–45) |
| R | 25% (20–30) | 11% (10–13) | 64% (58–70) |

## Matrix, all drugs, same culture, pair-weighted (reproduces 29 / 74 / 34%)

| first \ other | S | I | R |
| --- | --- | --- | --- |
| S | 71% (64–77) | 9% (7–10) | 21% (15–27) |
| I | 34% (29–41) | 26% (21–32) | 40% (34–46) |
| R | 23% (18–30) | 11% (9–13) | 66% (58–72) |

## Matrix, same culture, same morphotype (culture-weighted)

| first \ other | S | I | R |
| --- | --- | --- | --- |
| S | 79% (71–85) | 7% (5–10) | 14% (10–20) |
| I | 30% (24–36) | 33% (24–42) | 37% (29–46) |
| R | 17% (13–23) | 10% (8–13) | 73% (66–78) |

## Matrix, same culture, different morphotype (culture-weighted)

| first \ other | S | I | R |
| --- | --- | --- | --- |
| S | 69% (60–76) | 9% (7–11) | 22% (15–31) |
| I | 40% (33–48) | 20% (16–24) | 40% (33–47) |
| R | 30% (25–37) | 12% (10–14) | 58% (50–63) |

## Matrix, all drugs, next culture, multi-colony index cultures (culture-weighted)

| first \ other | S | I | R |
| --- | --- | --- | --- |
| S | 73% (67–79) | 8% (7–10) | 19% (14–24) |
| I | 35% (31–40) | 24% (21–28) | 40% (34–46) |
| R | 22% (18–27) | 11% (10–13) | 66% (61–71) |

## Matrix, all drugs, next culture, all index cultures (culture-weighted)

| first \ other | S | I | R |
| --- | --- | --- | --- |
| S | 77% (73–81) | 8% (7–9) | 15% (12–19) |
| I | 34% (30–39) | 27% (24–30) | 39% (34–44) |
| R | 18% (15–22) | 10% (8–12) | 72% (67–77) |

## Overall discordance, all drugs and by class (culture-weighted)

| comparison | stratum | patients | cultures | any S/I/R difference | S<->R | S vs R+I flip | R vs S+I flip |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Next culture | All drugs | 81 | 8587 | 30% (27–33) | 15% (13–17) | 22% (20–24) | 23% (20–25) |
| Next culture, multi-colony index | All drugs | 61 | 3703 | 35% (31–39) | 18% (15–20) | 26% (23–29) | 27% (23–31) |
| Same culture | All drugs | 65 | 3849 | 36% (31–41) | 19% (15–23) | 27% (23–31) | 28% (23–32) |
| Next culture | Aminoglycoside | 81 | 2892 | 29% (25–33) | 11% (9–14) | 20% (17–22) | 21% (17–24) |
| Next culture | Carbapenem | 81 |  965 | 25% (19–30) | 14% (11–17) | 19% (15–23) | 20% (15–24) |
| Next culture | Cephalosporin | 81 | 1934 | 32% (28–36) | 18% (15–21) | 25% (22–28) | 25% (21–29) |
| Next culture | Fluoroquinolone | 81 |  971 | 29% (24–33) | 10% (7–12) | 17% (12–21) | 22% (18–26) |
| Next culture | Monobactam | 73 |  873 | 32% (27–37) | 22% (19–26) | 28% (24–32) | 27% (22–31) |
| Next culture | Penicillin | 80 |  952 | 35% (29–41) | 18% (14–21) | 28% (23–33) | 25% (20–30) |
| Next culture, multi-colony index | Aminoglycoside | 61 | 1263 | 33% (28–38) | 14% (11–18) | 24% (20–27) | 24% (19–28) |
| Next culture, multi-colony index | Carbapenem | 60 |  420 | 30% (23–37) | 18% (13–22) | 24% (18–29) | 24% (18–29) |
| Next culture, multi-colony index | Cephalosporin | 61 |  843 | 38% (33–42) | 22% (18–27) | 30% (26–33) | 31% (26–36) |
| Next culture, multi-colony index | Fluoroquinolone | 60 |  422 | 36% (30–40) | 12% (9–16) | 18% (14–23) | 30% (24–34) |
| Next culture, multi-colony index | Monobactam | 57 |  334 | 36% (29–41) | 26% (21–31) | 31% (26–37) | 30% (24–35) |
| Next culture, multi-colony index | Penicillin | 61 |  421 | 39% (31–46) | 19% (15–24) | 31% (25–36) | 28% (21–34) |
| Same culture | Aminoglycoside | 65 | 1315 | 34% (29–39) | 16% (11–20) | 24% (20–29) | 25% (20–31) |
| Same culture | Carbapenem | 63 |  436 | 33% (24–40) | 20% (14–27) | 27% (20–34) | 26% (19–33) |
| Same culture | Cephalosporin | 65 |  879 | 39% (33–44) | 24% (18–29) | 32% (27–37) | 31% (24–37) |
| Same culture | Fluoroquinolone | 64 |  440 | 38% (32–44) | 14% (9–20) | 19% (13–26) | 33% (27–39) |
| Same culture | Monobactam | 61 |  344 | 38% (30–46) | 27% (19–35) | 33% (25–42) | 31% (23–39) |
| Same culture | Penicillin | 64 |  435 | 37% (28–45) | 19% (13–25) | 31% (23–38) | 25% (18–31) |

## By drug (culture-weighted)

| comparison | stratum | patients | cultures | any S/I/R difference | S<->R | S vs R+I flip | R vs S+I flip |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Next culture | amikacin | 81 | 966 | 29% (24–35) | 13% (10–16) | 20% (15–24) | 23% (18–27) |
| Next culture | aztreonam | 73 | 873 | 32% (27–37) | 22% (19–26) | 28% (24–32) | 27% (22–31) |
| Next culture | cefepime | 81 | 970 | 36% (31–41) | 15% (12–18) | 25% (21–29) | 26% (22–30) |
| Next culture | ceftazidime | 81 | 964 | 29% (24–34) | 20% (16–25) | 25% (21–29) | 24% (19–29) |
| Next culture | gentamicin | 81 | 965 | 32% (26–36) | 8% (6–11) | 18% (14–22) | 22% (18–26) |
| Next culture | levofloxacin | 81 | 971 | 29% (24–33) | 10% (7–12) | 17% (12–21) | 22% (18–26) |
| Next culture | meropenem | 81 | 965 | 25% (19–30) | 14% (11–17) | 19% (15–23) | 20% (15–24) |
| Next culture | piperacillin_tazobactam | 80 | 952 | 35% (29–41) | 18% (14–21) | 28% (23–33) | 25% (20–30) |
| Next culture | tobramycin | 81 | 961 | 26% (20–32) | 12% (9–16) | 21% (17–26) | 18% (13–22) |
| Next culture, multi-colony index | amikacin | 61 | 420 | 35% (28–41) | 15% (12–19) | 24% (18–29) | 26% (20–32) |
| Next culture, multi-colony index | aztreonam | 57 | 334 | 36% (29–41) | 26% (21–31) | 31% (26–37) | 30% (24–35) |
| Next culture, multi-colony index | cefepime | 61 | 424 | 41% (35–47) | 19% (15–23) | 29% (24–34) | 31% (26–36) |
| Next culture, multi-colony index | ceftazidime | 60 | 419 | 35% (28–41) | 26% (20–31) | 31% (25–36) | 30% (23–37) |
| Next culture, multi-colony index | gentamicin | 61 | 423 | 37% (30–42) | 11% (8–15) | 23% (17–28) | 25% (19–30) |
| Next culture, multi-colony index | levofloxacin | 60 | 422 | 36% (30–40) | 12% (9–16) | 18% (14–23) | 30% (24–34) |
| Next culture, multi-colony index | meropenem | 60 | 420 | 30% (23–37) | 18% (13–22) | 24% (18–29) | 24% (18–29) |
| Next culture, multi-colony index | piperacillin_tazobactam | 61 | 421 | 39% (31–46) | 19% (15–24) | 31% (25–36) | 28% (21–34) |
| Next culture, multi-colony index | tobramycin | 61 | 420 | 29% (21–37) | 16% (11–21) | 25% (18–32) | 20% (14–26) |
| Same culture | amikacin | 65 | 437 | 35% (28–41) | 15% (11–19) | 24% (19–31) | 25% (20–31) |
| Same culture | aztreonam | 61 | 344 | 38% (30–46) | 27% (19–35) | 33% (25–42) | 31% (23–39) |
| Same culture | cefepime | 65 | 442 | 42% (34–50) | 22% (16–29) | 32% (25–39) | 32% (25–39) |
| Same culture | ceftazidime | 64 | 437 | 36% (27–43) | 25% (18–32) | 31% (24–39) | 29% (21–37) |
| Same culture | gentamicin | 65 | 440 | 38% (31–45) | 13% (9–18) | 24% (18–30) | 28% (21–34) |
| Same culture | levofloxacin | 64 | 440 | 38% (32–44) | 14% (9–20) | 19% (13–26) | 33% (27–39) |
| Same culture | meropenem | 63 | 436 | 33% (24–40) | 20% (14–27) | 27% (20–34) | 26% (19–33) |
| Same culture | piperacillin_tazobactam | 64 | 435 | 37% (28–45) | 19% (13–25) | 31% (23–38) | 25% (18–31) |
| Same culture | tobramycin | 65 | 438 | 29% (20–38) | 19% (11–26) | 25% (17–33) | 23% (15–31) |

## By morphotype pairing (culture-weighted)

| comparison | stratum | patients | cultures | any S/I/R difference | S<->R | S vs R+I flip | R vs S+I flip |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Same culture | different morphotype | 50 | 2288 | 41% (35–47) | 23% (17–28) | 32% (27–37) | 32% (26–38) |
| Same culture | same morphotype | 53 | 1561 | 29% (24–34) | 14% (11–17) | 21% (16–25) | 22% (18–26) |
| Same culture | Aminoglycoside / different morphotype | 50 |  781 | 40% (33–46) | 19% (13–25) | 29% (23–35) | 30% (22–38) |
| Same culture | Aminoglycoside / same morphotype | 53 |  534 | 26% (20–32) | 11% (7–15) | 17% (14–22) | 19% (14–25) |
| Same culture | Carbapenem / different morphotype | 49 |  261 | 39% (28–48) | 23% (15–33) | 31% (20–40) | 32% (22–41) |
| Same culture | Carbapenem / same morphotype | 52 |  175 | 24% (15–32) | 14% (9–21) | 21% (13–29) | 17% (11–23) |
| Same culture | Cephalosporin / different morphotype | 50 |  525 | 44% (38–51) | 28% (21–35) | 38% (32–44) | 35% (27–42) |
| Same culture | Cephalosporin / same morphotype | 53 |  354 | 31% (22–39) | 17% (11–23) | 23% (16–30) | 25% (17–33) |
| Same culture | Fluoroquinolone / different morphotype | 50 |  262 | 43% (34–49) | 16% (10–23) | 22% (15–31) | 36% (29–42) |
| Same culture | Fluoroquinolone / same morphotype | 52 |  178 | 32% (24–40) | 10% (5–17) | 15% (9–22) | 28% (20–35) |
| Same culture | Monobactam / different morphotype | 45 |  201 | 44% (33–54) | 32% (22–41) | 39% (28–50) | 37% (27–46) |
| Same culture | Monobactam / same morphotype | 51 |  143 | 29% (20–38) | 19% (12–29) | 25% (17–35) | 23% (15–32) |
| Same culture | Penicillin / different morphotype | 49 |  259 | 36% (25–46) | 21% (14–30) | 34% (23–43) | 24% (16–33) |
| Same culture | Penicillin / same morphotype | 53 |  176 | 38% (26–48) | 15% (10–22) | 27% (19–35) | 26% (17–34) |

## Different minus same morphotype

| weighting | stat | same | different | diff | P(diff <= 0) |
| --- | --- | --- | --- | --- | --- |
| culture | discord | 0.29 | 0.41 | 0.12 (0.05–0.19) | 0.000 |
| culture | flip_S_vs_NS | 0.21 | 0.32 | 0.11 (0.06–0.17) | 0.000 |
| pairs | discord | 0.27 | 0.44 | 0.16 (0.07–0.24) | 0.000 |
| pairs | flip_S_vs_NS | 0.19 | 0.34 | 0.15 (0.07–0.21) | 0.000 |

## Same culture / next culture discordance ratio

| stat | stratum | same | next, all | ratio, all | next, multi-colony index | ratio, matched |
| --- | --- | --- | --- | --- | --- | --- |
| discord | All drugs | 0.36 | 0.30 | 1.19 (1.07–1.31) | 0.35 | 1.03 (0.95–1.10) |
| discord | Aminoglycoside | 0.34 | 0.29 | 1.17 (1.05–1.31) | 0.33 | 1.02 (0.94–1.11) |
| discord | Carbapenem | 0.33 | 0.25 | 1.30 (1.07–1.55) | 0.30 | 1.09 (0.96–1.23) |
| discord | Cephalosporin | 0.39 | 0.32 | 1.20 (1.04–1.36) | 0.38 | 1.02 (0.92–1.12) |
| discord | Fluoroquinolone | 0.38 | 0.29 | 1.34 (1.12–1.60) | 0.36 | 1.08 (0.96–1.19) |
| discord | Monobactam | 0.38 | 0.32 | 1.17 (0.94–1.42) | 0.36 | 1.06 (0.90–1.21) |
| discord | Penicillin | 0.37 | 0.35 | 1.05 (0.87–1.23) | 0.39 | 0.95 (0.83–1.06) |
| flip_S_vs_NS | All drugs | 0.27 | 0.22 | 1.23 (1.09–1.38) | 0.26 | 1.05 (0.96–1.14) |
| flip_S_vs_NS | Aminoglycoside | 0.24 | 0.20 | 1.24 (1.09–1.44) | 0.24 | 1.02 (0.90–1.13) |
| flip_S_vs_NS | Carbapenem | 0.27 | 0.19 | 1.38 (1.06–1.73) | 0.24 | 1.12 (0.94–1.32) |
| flip_S_vs_NS | Cephalosporin | 0.32 | 0.25 | 1.27 (1.10–1.44) | 0.30 | 1.08 (0.95–1.19) |
| flip_S_vs_NS | Fluoroquinolone | 0.19 | 0.17 | 1.17 (0.79–1.60) | 0.18 | 1.05 (0.82–1.25) |
| flip_S_vs_NS | Monobactam | 0.33 | 0.28 | 1.20 (0.95–1.48) | 0.31 | 1.06 (0.87–1.24) |
| flip_S_vs_NS | Penicillin | 0.31 | 0.28 | 1.12 (0.89–1.34) | 0.31 | 1.01 (0.86–1.16) |
