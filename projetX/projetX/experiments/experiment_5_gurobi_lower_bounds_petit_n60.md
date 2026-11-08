# Experiment 5: Gurobi against the dynamic-program and dual lower bounds on Petit_N60

> Run on 8 November 2026 on the 18 instances `Petit_N60_*`. Raw values: [results.csv](../output/gurobi_lower_bounds_n60/results.csv).

## What was tested

Each instance has 8 to 16 loads and 1 to 3 sources. Limits: 8 MW and 50 km per feeder, 70 MW and at most 8 feeders per source.

| Method | What it computes | Setup |
|---|---|---|
| Dynamic program (DP) | Best split of the loads into feeders, within the power, length, source-capacity and feeder-count limits. Ignores node disjointness between feeders, so it is a lower bound. | Dreyfus-Wagner on the graph without sources, then a partition over load subsets |
| Dual bound | Dual of the LP with one variable per allowed feeder and source, each reduced cost rechecked | HiGHS through SciPy |
| Gurobi | PLNE of `plne_reseau.py` aligned on `verify.py`, tightened with the valid inequalities below | 4 threads, 1,200 s limit, MIPGap 10⁻⁶, no starting solution |

The Gurobi model is independent from the two bounds. Its tightening keeps every optimal solution:

- at most min(`n_departs_max`, degree of the source) feeders per source;
- `f ≥ x` on every arc;
- a transit node used by a feeder has a child;
- one unit of flow from the sources to each load on the arcs of the network, summed over the feeders.

## Results

"DP status" is **exact** when the feeders found are node-disjoint and `verify.py` accepts the network, which proves optimality on its own. Bold values are strictly below the optimum.

| Instance | Loads | Sources | DP bound | DP status | Dual bound | Gurobi optimum | Gurobi time (s) |
|---|---|---|---|---|---|---|---|
| S1_0.15_Ete | 8 | 1 | 76.31 | exact | 76.31 | 76.31 | 82 |
| S1_0.15_Hiver | 9 | 1 | **68.22** | bound only | **66.71** | 69.60 | 365 |
| S1_0.20_Ete | 12 | 1 | 68.72 | exact | 68.72 | 68.72 | 1 |
| S1_0.20_Hiver | 12 | 1 | 73.62 | exact | **72.12** | 73.62 | 855 |
| S1_0.25_Ete | 16 | 1 | 91.94 | exact | **90.26** | 91.94 | 136 |
| S1_0.25_Hiver | 15 | 1 | 85.22 | exact | 85.22 | 85.22 | 267 |
| S2_0.15_Ete | 8 | 2 | 91.19 | exact | 91.19 | 91.19 | 2 |
| S2_0.15_Hiver | 9 | 2 | **91.10** | bound only | **91.10** | 95.79 | 83 |
| S2_0.20_Ete | 12 | 2 | 88.13 | exact | 88.13 | 88.13 | 7 |
| S2_0.20_Hiver | 12 | 2 | 101.42 | exact | 101.42 | 101.42 | 8 |
| S2_0.25_Ete | 16 | 2 | 108.65 | exact | 108.65 | 108.65 | 6 |
| S2_0.25_Hiver | 15 | 2 | 136.67 | exact | 136.67 | 136.67 (not proven) | 1,200 (limit) |
| S3_0.15_Ete | 8 | 3 | 70.51 | exact | 70.51 | 70.51 | 1 |
| S3_0.15_Hiver | 9 | 3 | 100.29 | exact | 100.29 | 100.29 | 3 |
| S3_0.20_Ete | 12 | 3 | 138.84 | exact | 138.84 | 138.84 | 6 |
| S3_0.20_Hiver | 12 | 3 | 128.27 | exact | 128.27 | 128.27 | 11 |
| S3_0.25_Ete | 16 | 3 | 136.88 | exact | 136.88 | 136.88 | 7 |
| S3_0.25_Hiver | 15 | 3 | 152.14 | exact | 152.14 | 152.14 | 22 |

| | DP bound | Dual bound | Gurobi |
|---|---|---|---|
| Equal to the optimum | 16 / 18 | 14 / 18 | 17 / 18 proven |
| Largest shortfall | 4.69 km | 4.69 km | none on proven instances |
| Time for all 18 | 29 s together, 8 processes | included in the 29 s | 24 min, 4 instances at a time (51 min of solve time) |

Both bounds are never above the Gurobi optimum, as expected of valid lower bounds.

## Why the bounds differ

- **DP below the optimum (2 instances).** In `S1_0.15_Hiver` and `S2_0.15_Hiver`, the best split has two feeders starting at the same node (24 and 40). A radial network forbids this, and the DP does not see it.
- **Dual below the optimum (4 instances).** On `S1_0.15_Hiver`, `S1_0.20_Hiver` and `S1_0.25_Ete` the feeder LP solution is fractional, so the dual stays 1.5 to 2.9 km under the optimum. On `S2_0.15_Hiver` the LP solution is integer and equals the DP bound, so the dual inherits the DP's feeder overlap and is 4.69 km low.

## Why S2_0.25_Hiver reached the time limit

Gurobi found the optimal 136.67 km after about 40 s. The remaining time went to raising the lower bound, which reached only 129.62 km.

- On all 18 instances the root LP bound of the PLNE equals the optimum with no feeder limit (121.72 km here). The LP spreads each load across the feeders, so the 8 MW and 50 km limits never bind.
- The feeder limits add 14.95 km on this instance, the largest of the 18, and Gurobi must recover them by branching. It explored about 11 nodes per second, and no cutting plane moved the root bound.
- The dual over complete feeders reaches 136.67 km at once with an integer solution, so the difficulty comes from the edge-by-edge formulation.

## Conclusions

- **The DP and dual bounds are exact on most small instances**, and all 18 are done in 29 s instead of 24 minutes.
- **The DP proves optimality on its own on 16 of 18 instances**, including the one Gurobi could not close.
- **Gurobi remains the reference when feeders overlap:** the DP alone cannot settle `S1_0.15_Hiver` or `S2_0.15_Hiver`.
- **Stopping Gurobi once it reaches the DP bound is valid**, since a feasible network at a lower bound is optimal. Tested after this campaign, the rule proves `S2_0.25_Hiver` optimal in 29 s of Gurobi time instead of timing out at 1,200 s.

## Limits

- Only 8 to 16 loads were tested. The DP costs about 3^(loads) and does not scale to Moyen instances (22 to 62 loads).
- `verify.py` and the PLNE check the power only at the head of each feeder and do not check P×L. On `S2_0.25_Ete` the optimum carries 9 MW on one line (limit 8) and reaches P×L = 231 (limit 200). So 108.65 km is the optimum of the repository's model, not of the file's full limits. The other 17 networks respect both limits.
- The timing comparison used 4 Gurobi threads. Other thread counts change the times.
- This table was produced by temporary scripts. `subset_dp.py` was written afterwards and gives the same 18 DP bounds to the fourth decimal. The dual bound is not in the repository.

## Reproduce

The tightened Gurobi model and the DP bound are in the repository since this experiment. From `projetX/projetX`:

```powershell
python plne_gurobi.py verified --prefix Petit_N60 --tighten --dp-stop --no-warm-start --threads 4 --time-limit 1200
```

`--dp-stop` is the improved stop rule, so this command ends earlier than the run reported here on the instances the DP solves. Without it, the command repeats the setup reported above, one instance at a time.
