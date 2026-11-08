# Final report notes

> Running summary of the results that matter for the final report, kept very short. Details are in [experiments/](experiments/).

## Instance set

**26 of the 216 verified instances are provably infeasible** and are now in `data/verified_instances/unfeasible/`. The 190 remaining instances pass the two checks below, which does not prove they are solvable.

| Class | Instances | Infeasible | Remaining |
|---|---|---|---|
| Petit | 54 | 0 | 54 |
| Moyen | 54 | 4 | 50 |
| Grand | 54 | 5 | 49 |
| TrèsGrand | 54 | 17 | 37 |

- **Power test (7 fail, all winter):** |sum of demands| ≤ sum over sources of min(source capacity, max feeders × feeder power), which is min(70, 8 × 8) = 64 MW per source.
- **Distance test (19 fail):** every Conso and Prod node is within the 50 km feeder length of its nearest source along the graph. The farthest load is at 57.91 km.

Full list: [Experiment 3](experiments/experiment_3_unfeasible_instances.md).

## Lower bounds

**Fast bounds are weak.** The H0 bound (under 20 ms) reaches 53-55% of the optimum on two instances with a proven optimum, and ½ (MST + longest MST edge) 57-60%. [Experiment 4](experiments/experiment_4_lower_bounds.md).

**Feeder-based bounds are nearly exact on small instances.** On the 18 Petit_N60 instances, they are computed in 29 s in total, against 24 minutes for Gurobi. The dynamic program (DP) equals the optimum on 16 of 18, and the dual bound on 14 of 18. [Experiment 5](experiments/experiment_5_gurobi_lower_bounds_petit_n60.md).

## How far exact methods go

**Optimality is proven up to about 20 loads, and rarely beyond.** Gurobi stops as soon as it reaches the DP bound, which is valid proof of optimality. Limit: 20 minutes per run.

| Instances | Loads | Runs | Proven optimal |
|---|---|---|---|
| Petit_N60 | 8-16 | 18 | 18 |
| Petit_N75 and N90 (random) | 14-22 | 10 | 8 |
| Moyen_N150 | 22-30 | 8 | 1 |

The other Moyen runs end with a gap of 1.9 to 24%, or no network at all (1 run). The DP costs about 3^(loads), so it stops at 20 loads (584 s at 19 loads). [Experiment 6](experiments/experiment_6_gurobi_campaign_petit_moyen.md).

## TDA techniques

**Merge-tree feeders are the best fast heuristic:** valid on 50 of 54 Petit instances in 0.2 s (Voronoi: 23), and about 4% above the optimum after local search. **Persistence cuts and a TDA warm start help Gurobi:** on Petit_N60, 14 of 18 proven instead of 10 at equal work, with 45% less work. [Experiment 7](experiments/experiment_7_tda_heuristics_and_plne.md).

## Scaling to TrèsGrand (proposal, not implemented)

Exact methods cannot reach 300 to 1,000 loads, but a feeder stays small at any size (about 5 to 20 loads). The proposal re-optimizes areas of 2 or 3 feeders exactly, recombines the feeders found, and certifies a gap with column generation. [Scaling](docs/scaling.md).

## Known gap

`verify.py` and the PLNE do not check per-line power or P×L. One Petit_N60 optimum breaks both limits (`S2_0.25_Ete`).

Methods: [exact](docs/exact_methods_explained.md), [TDA](docs/tda_algorithms_explained.md).
