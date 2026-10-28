# Pseudo-trajectories: one colony per culture

2000 random picks per patient x class; patient bootstrap B = 2000.
Colony frame, class-level NS (any agent R or I). Patients with >= 2 cultures in the class.
History types: stable S, stable NS, acquired (one S->NS), lost (one NS->S), fluctuating (>= 2 switches).
Lower bound: single-colony cultures enter as fixed.

| variant | class | patients | cultures | patients with a heterogeneous culture | switches/patient, one colony | switches/patient, all colonies any-NS | share of switches set by the pick | patients whose switch count varies | SD of switch count across picks | P(two records differ by >= 2 switches) | patients with >= 2 histories at >= 10% | P(one-colony history = all-colony history) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| All cultures | Aminoglycoside | 81 | 1052 | 43% (33–54) | 1.92 (1.41–2.47) | 1.58 (1.14–2.04) | 48% (39–57) | 40% (30–51) | 0.43 (0.31–0.57) | 16% (11–21) | 16% (9–25) | 89% (84–94) |
| All cultures | Carbapenem | 81 | 1046 | 44% (33–56) | 2.35 (1.74–2.96) | 1.91 (1.40–2.48) | 51% (40–61) | 42% (31–53) | 0.51 (0.37–0.66) | 21% (15–27) | 11% (5–19) | 89% (83–95) |
| All cultures | Cephalosporin | 81 | 1052 | 49% (38–60) | 3.00 (2.37–3.68) | 2.36 (1.78–2.95) | 47% (37–56) | 49% (38–60) | 0.63 (0.47–0.78) | 23% (17–28) | 15% (7–23) | 88% (82–94) |
| All cultures | Fluoroquinolone | 81 | 1052 | 42% (32–53) | 1.98 (1.42–2.67) | 1.46 (0.94–2.11) | 43% (30–59) | 41% (31–52) | 0.45 (0.33–0.59) | 18% (13–24) | 21% (12–31) | 84% (77–91) |
| Multi-colony cultures only | Aminoglycoside | 53 |  431 | 64% (51–77) | 1.52 (1.07–2.07) | 1.04 (0.60–1.57) | 83% (71–92) | 58% (45–72) | 0.59 (0.44–0.77) | 18% (12–25) | 40% (26–53) | 73% (64–82) |
| Multi-colony cultures only | Carbapenem | 53 |  426 | 68% (55–81) | 1.83 (1.31–2.45) | 1.15 (0.75–1.60) | 89% (81–97) | 66% (53–79) | 0.75 (0.58–0.92) | 27% (20–34) | 40% (26–53) | 66% (56–76) |
| Multi-colony cultures only | Cephalosporin | 53 |  431 | 68% (55–79) | 2.04 (1.46–2.64) | 1.08 (0.64–1.57) | 90% (83–96) | 68% (55–79) | 0.81 (0.63–0.99) | 26% (20–33) | 42% (28–55) | 64% (53–74) |
| Multi-colony cultures only | Fluoroquinolone | 53 |  429 | 62% (49–75) | 1.30 (0.91–1.72) | 0.55 (0.28–0.83) | 88% (79–96) | 62% (49–75) | 0.62 (0.47–0.77) | 21% (14–28) | 43% (30–57) | 67% (58–77) |
