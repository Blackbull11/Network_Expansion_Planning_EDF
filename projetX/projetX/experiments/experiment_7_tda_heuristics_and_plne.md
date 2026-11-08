# Experiment 7: TDA heuristics, and TDA techniques inside the exact PLNE

> Run on 8 November 2026 with [run_tda_plne.py](../run_tda_plne.py). Raw values: `output/tda_plne/heuristics.csv` and `output/tda_plne/exact.csv`. Methods are explained in [tda_algorithms_explained.md](../docs/tda_algorithms_explained.md).

## What was tested

**Heuristics:** every TDA heuristic on the 54 Petit instances, each result judged by `verify.py`.

**Exact:** the tightened Gurobi model on the 18 Petit_N60 instances, in five variants, with 1 thread and the same limit of 120 Gurobi work units. Work units do not depend on the load of the machine, unlike seconds.

| Variant | Description |
|---|---|
| base | Tightened model, no starting solution |
| warm | Starts from the best valid TDA heuristic |
| cuts | Adds rounded capacity cuts on persistence clusters (`TDA_algos/plne_cuts.py`), 8.5 per instance on average |
| warm+cuts | Both |
| pipeline | Solves the PLNE on a restricted set of edges (27% kept on average), then runs warm+cuts from that solution. Its work includes the restricted solve |

## Results: heuristics on 54 Petit instances

"Gap" is the mean excess length over the proven optimum, measured only on the Petit_N60 instances where the heuristic is valid, so the instances differ between rows.

| Heuristic | Valid (of 54) | Mean time (s) | Valid on N60 (of 18) | Mean gap on N60 |
|---|---|---|---|---|
| Portfolio (best of all below) | **51** | 3.13 | 18 | 4.2% |
| Merge tree + local search (all candidates) | 50 | 0.20 | 17 | 3.9% |
| Merge tree + local search (TDA candidates) | 50 | 0.19 | 17 | 4.1% |
| Merge tree feeders | 50 | 0.18 | 17 | 9.6% |
| Savings feeders | 49 | 0.21 | 16 | 9.7% |
| Kruskal feeders | 44 | 0.39 | 16 | 13.2% |
| Voronoi variable distance | 23 | 0.44 | 9 | 10.3% |
| ToMATo initial weights | 23 | 0.61 | 9 | 8.2% |
| Spectral edge costs | 23 | 0.46 | 9 | 5.6% |
| Shrink + grow | 23 | 0.45 | 9 | 11.0% |
| Anchored spectral cells | 11 | 0.01 | 6 | 13.9% |

No heuristic finds a valid network on 3 instances: `Petit_N75_S1_0.20_Hiver`, `Petit_N90_S1_0.20_Hiver` and `Petit_N90_S1_0.25_Hiver`. None is proven infeasible, and Gurobi did find a valid network on the first one ([Experiment 6](experiment_6_gurobi_campaign_petit_moyen.md)).

## Results: exact variants on 18 Petit_N60 instances

Runs that hit the work limit count as 120 in the mean work.

| Variant | Proven optimal (of 18) | Mean work | Median work | Mean final gap | Work vs base (geometric mean) | Better / worse than base |
|---|---|---|---|---|---|---|
| base | 10 | 57.5 | 19.7 | 4.5% | 1 | |
| warm | 11 | 53.6 | 16.8 | 3.7% | 0.82 | 9 / 1 |
| cuts | 13 | 44.7 | 10.0 | 1.8% | 0.78 | 8 / 2 |
| warm+cuts | 13 | 41.2 | 7.8 | 1.7% | 0.60 | 12 / 0 |
| pipeline | **14** | 38.8 | 6.2 | 1.6% | **0.55** | 11 / 1 |

Better and worse mean a work difference of more than 5%. The restricted solve alone found the optimum on 9 of 18 instances.

Checks:

- All 90 runs ended with a network accepted by `verify.py`.
- Every proven optimum equals the one in [Experiment 5](experiment_5_gurobi_lower_bounds_petit_n60.md).
- No persistence cut was violated by a base solution.

## Conclusions

- **Designing feeders beats cutting Voronoi cells.** Merge-tree feeders are valid on 50 of 54 Petit instances in 0.2 s, against 23 of 54 for Voronoi variable distance. The local search brings them from 9.6% to about 4% above the optimum.
- **The TDA techniques help the exact solver.** Persistence cuts raise the proven instances from 10 to 13 of 18, and the full pipeline to 14 of 18, with 45% less work in geometric mean.
- **The warm start alone gives little** (11 of 18), but it combines well with the cuts.

## Limits

- 18 small instances with a single thread and one work limit. The comparison shows a trend.
- The four instances the pipeline leaves unproven (`S1_0.25_Ete`, `S1_0.25_Hiver`, `S1_0.15_Hiver`, `S2_0.25_Hiver`) include 3 single-source instances, and these variants were not tested on Moyen.
- The cuts are checked against the base solutions only. They were not checked against every optimal network.

## Reproduce

From `projetX/projetX`:

```powershell
python run_tda_plne.py heuristics
python run_tda_plne.py exact --work-limit 120 --workers 2
python run_tda_plne.py report
```
