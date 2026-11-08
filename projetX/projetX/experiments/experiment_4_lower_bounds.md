# Experiment 4: lower bounds on J

> Run on 8 November 2026 on the 190 instances of `data/verified_instances/`, with `lower_bound` from [bounds.py](../TDA_algos/bounds.py). The moat bound below was computed by a temporary script and is not in the repository. The subset dynamic program is now [subset_dp.py](../subset_dp.py), which also enforces the source capacity and the number of feeders.

## Bounds compared

| Bound | Formula | What it ignores |
|---|---|---|
| `lb_nn` (repo) | ½ Σ distance to nearest other point | L, P, number of feeders, node disjointness |
| `lb_mst` (repo) | MST · t / (2(t − 1)) | same |
| Moat (new) | ½ (MST + longest MST edge) | same |
| Steiner* (DP) | Dreyfus-Wagner on loads + source | same, but solved exactly |
| Feeder partition (DP) | best partition of the loads into feeders within L and P | node disjointness between feeders, source capacity |

**Moat bound proof.** An Euler tour of the optimal Steiner tree, shortcut to the terminals, is a cycle of length at most 2 Steiner*. It crosses the cut of the longest MST edge at least twice, each crossing being at least that edge. Deleting one crossing leaves a spanning path, of length at least MST. So 2 Steiner* ≥ MST + longest edge.

## Results on all 190 instances

- `lb_mst` is the larger repo bound on all 190 instances (`lb_nn` never beats it).
- The moat bound is at least `max(lb_mst, lb_nn)` on all 190, by 3.5% on average and 20.2% at most.

| Class | Instances | Nodes | Points (loads + root) | Mean time | Max time |
|---|---|---|---|---|---|
| Petit | 54 | 60-90 | 9-23 | 0.6 ms | 0.9 ms |
| Moyen | 50 | 150-250 | 23-63 | 1.1 ms | 2.1 ms |
| Grand | 49 | 500-700 | 76-177 | 2.6 ms | 3.6 ms |
| TrèsGrand | 37 | 2000-4000 | 301-1001 | 12.3 ms | 18.3 ms |

## Quality against the proven optimum

| | Petit_N60_S1_0.25_Hiver | % of J* | Time | Petit_N60_S1_0.20_Ete | % of J* |
|---|---|---|---|---|---|
| `lb_nn` | 37.20 | 44% | under 1 ms | 32.42 | 47% |
| `lb_mst` | 45.13 | 53% | under 1 ms | 38.14 | 55% |
| Moat | 48.20 | 57% | under 1 ms | 41.23 | 60% |
| HiGHS bound at 600 s | 65.85 | 77% | 600 s | | |
| Steiner* (DP) | 79.38 | 93% | 1.8 s | | |
| Feeder partition (DP) | 85.22 | 100% | 5.2 s | | |
| Gurobi optimum J* | 85.22 | 100% | 295 s | 68.72 | 100% |

The dynamic programs were run on Petit_N60_S1_0.25_Hiver only (15 loads).

## Conclusions

- The H0 and moat bounds are extremely cheap, and only reach 53-60% of J* on the two instances with a proven optimum.
- On Petit_N60_S1_0.25_Hiver, the feeder limits account for 5.84 km of the gap (79.38 to 85.22). The rest comes from the weakness of the fast bound, which is 31.18 km below Steiner*.
- The partition DP matches the optimum on this instance. Its cost is about 3^k in the number of loads k, so it applies to the smallest instances only: Petit goes up to 22 loads (3^22 ≈ 3·10^10).

## Limits

- Only two instances had a proven optimum when this was run, so the percentages show a trend. [Experiment 5](experiment_5_gurobi_lower_bounds_petit_n60.md) compares the partition DP with the optimum on all 18 Petit_N60 instances.
- The partition DP is a lower bound in general. It equals J* here because Gurobi proves the same value.
