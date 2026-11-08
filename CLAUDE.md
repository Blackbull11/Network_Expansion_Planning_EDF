# Network Expansion Planning (EDF): project log

Project memory: how to work here, what has been done, what is known and what is open. Run instructions are in [README.md](README.md).

Last updated: 2026-11-08, after the Gurobi campaign and the TDA study.

## Working rules

- **Run from `projetX/projetX`.** Interpreter: `C:\Users\siand\dev\X_3A\RO_EDF\.venv\Scripts\python.exe`. `python` on PATH is the Microsoft Store stub. Set `PYTHONUTF8=1`, because instance names contain `è`.
- **`verify.py` is the only judge of feasibility.** Never report a J without its `verify.py` result.
- **After each result, update three places:**
  1. `experiments/experiment_N_<topic>.md`: method, full tables, limits, reproduction command.
  2. `projetX/projetX/report.md`: conclusions only, very short. Replace outdated statements.
  3. The tables in this file: one line per experiment, and move finished work out of "In progress".
- **Separate measured values from estimates.** State the instance count, the time limit and the thread count with each comparison.
- **Do not write into an output folder while a job is using it.** Check for running jobs before relaunching a campaign.
- **Exclude `unfeasible/` from statistics.** Batches read only the top level of `data/verified_instances/`.

## Problem

Choose lines in a candidate graph to build a forest with one radial network per source, connecting every load (`Conso`, `Prod`) at minimum total length J. A feeder (*départ*) is the subtree hanging from one neighbour of the source.

Limits on the verified instances: 50 km and 8 MW per feeder, 70 MW and 8 feeders per source. Feeders are node-disjoint.

Vocabulary: *poste source* = source, *charge* = load, *départ* = feeder, *PLNE* = MILP.

## Code map

| Area | Files |
|---|---|
| Data | `nx_graph.py` (`GraphData`), `src/instances.py`, `data/generate_instances.py`, `plot_demands.py` |
| Judge | `verify.py`: `check_solution`, `necessary_conditions`, `solution_key` |
| Heuristics | `voronoi.py`, `voronoi_variable_distance.py`, `steiner.py`, `TDA_algos/` |
| Exact | `plne_reseau.py` (PuLP/HiGHS), `plne_gurobi.py` (model, `--tighten`, `--dp-stop`), `subset_dp.py` |
| Runners | `main.py`, `run_voronoi_verified.py`, `run_plne_verified.py`, `run_gurobi_campaign.py`, `run_tda_plne.py`, `run_tomato_decomposition.py`, `TDA.py` |
| Explanations | `docs/exact_methods_explained.md`, `docs/tda_algorithms_explained.md`, `docs/scaling.md`, `docs/voronoi_algorithms.tex`, `docs/reflexion.md` (whiteboard notes) |

## Project timeline

| Date | Change |
|---|---|
| 2026-09-27 | Initial network expansion model and 4-substation example |
| 2026-10-04 | Voronoi + Steiner with variable distance, tested on 40 generated graphs (Experiment 1) |
| 2026-10-06 | Jean-Côme's graphs added as JSON |
| 2026-11-02 | Verified instance dataset added (216 instances) |
| 2026-11-03 | Repository cleaned, TDA algorithms added, Voronoi run on the verified instances, first HiGHS PLNE runs |
| 2026-11-08 | Full Gurobi licence. Experiments 2 to 7, `subset_dp.py`, tightened model, DP stop rule. Instances renamed (no more ` (1)`), 26 infeasible instances moved out. ToMATo decomposition, persistence cuts, Gurobi campaign. Scaling proposal written. **Not committed yet.** |

## Experiments

| # | Date | Question | Main result | Status | Details |
|---|---|---|---|---|---|
| 1 | 10-04 | Which Voronoi distance on g01 to g40? | Variable distance is valid on 31/40, against 20/40 for the best fixed distance. Fixed-distance methods are not valid on the 9 other instances either. | Done | [experiment_1](projetX/projetX/experiments/experiment_1.md) |
| V | 11-03 | Voronoi + Steiner on verified instances | Valid on 28/213 at best, nearly none beyond Petit. | Stale, see below | `output/voronoi_verified/report.md` |
| H | 11-03 | HiGHS PLNE on Petit_N60 | No proof in 600 s, and often no integer solution. | Superseded by 2 | `output/plne_verified/` |
| 2 | 11-08 | Gurobi or HiGHS? | One instance: Gurobi proves 85.22 km in 295 s. HiGHS stops at 91.71 km with a 28% gap at 600 s. | Done | [experiment_2](projetX/projetX/experiments/experiment_2_plne_gurobi_vs_highs.md) |
| 3 | 11-08 | Are all instances feasible? | 26 of 216 are provably infeasible (7 on power, 19 on distance). The other 190 are only unproven. | Done | [experiment_3](projetX/projetX/experiments/experiment_3_unfeasible_instances.md) |
| 4 | 11-08 | How good are cheap lower bounds? | Under 20 ms but weak: 53 to 60% of J* on two instances. The moat bound beats the H0 bounds on all 190. | Done | [experiment_4](projetX/projetX/experiments/experiment_4_lower_bounds.md) |
| 5 | 11-08 | Gurobi against DP and dual on the 18 Petit_N60 | DP equals J* on 16/18, dual on 14/18, 29 s in total against 24 min for Gurobi. Gurobi proves 17/18. | Done | [experiment_5](projetX/projetX/experiments/experiment_5_gurobi_lower_bounds_petit_n60.md) |
| 6 | 11-08 | How far does Gurobi prove optimality with the DP stop rule? | 1,200 s per run, 5 threads. Petit N75/N90: 8/10 proven (7 at the DP bound). Moyen N150: 1/8 proven, gaps of 1.9 to 24%, no network on 1 run. **Moyen is the current limit of the exact approach.** | Done | [experiment_6](projetX/projetX/experiments/experiment_6_gurobi_campaign_petit_moyen.md) |
| 7 | 11-08 | Do TDA techniques help heuristics and the PLNE? | 54 Petit: merge tree valid on 50, Voronoi variable on 23. 18 Petit_N60, 1 thread, work limit 120: 10 proven for base, 14 for the TDA pipeline, with 45% less work. | Done | [experiment_7](projetX/projetX/experiments/experiment_7_tda_heuristics_and_plne.md) |
| T | 11-08 | ToMATo decomposition on one TrèsGrand | Only 1 significant density peak, best configuration still has 10 violations. Run on an instance later proven infeasible. | Done, negative | `output/tomato_decomposition/` |

## In progress

Nothing is running.

## What we know

- **Feeder limits make the problem hard.** Pure geometric methods give short networks that violate them. Every network below the optimum in Experiment 2 is invalid.
- **Designing feeders beats cutting Voronoi cells.** Merge-tree feeders are valid on far more Petit instances than Voronoi methods.
- **Up to about 20 loads, the DP is the best tool.** It often proves the optimum alone. A feasible network at the DP bound is optimal, which is the `--dp-stop` rule.
- **The weak point of the PLNE is its LP.** The root bound ignores feeder limits, since the LP spreads loads over feeders. The feeder LP (dual) closes this immediately on small instances.
- **Topology helps the exact solver.** Persistence cuts and a TDA warm start prove more Petit_N60 instances for less work (Experiment 7).
- **Density clustering does not help as a decomposition.** Instances have no natural clusters and the regions stay too large.
- **Local structure should allow scaling.** A feeder covers about 5 to 20 loads at every size, which motivates the proposal in `docs/scaling.md` (LNS, feeder pool, column generation). Nothing from it is implemented.

## Known issues

- **`src/instances.py` has a wrong path.** `VERIFIED_DIR` is `verified_instances/` instead of `data/verified_instances/`, so `verified_files()` returns nothing. This breaks `main.py --verified`, `TDA.py --set verified` and the batch mode of `run_plne_verified.py`.
- **Model gap:** `verify.py` and the PLNE check power only at the head of each feeder and do not check P×L. The optimum of `Petit_N60_S2_0.25_Ete` carries 9 MW on one line (limit 8) and reaches P×L = 231 (limit 200).
- **Stale results:** `output/voronoi_verified/` (213 instances, ` (1)` names, 25 infeasible counted against 26 today, some class labels wrong) and `output/plne_verified/`. The Voronoi counts in `docs/scaling.md` come from that run.
- **Lost results:** the Experiment 1 outputs (`recap_*`, charts) were removed in commit 1ec1776. `data/instances/` must be regenerated.
- **Scripts not in the repository:** the moat bound, the dual bound, the Experiment 3 sweep and the all-methods figure.
- **`run_tomato_decomposition.py` default instance** is now in `unfeasible/`, so a motif must be passed.
- **No `.gitignore`:** `.pyc` files and `output/normalized_instances/` are tracked.

## Next steps

1. Fix `VERIFIED_DIR`, then rerun Voronoi on the 190 clean instances.
2. Decide whether per-line power and P×L belong in `verify.py` and the PLNE. This can change reported optima.
3. Try the first steps of `docs/scaling.md` on Moyen: merge-tree start, then LNS with `subset_dp` and Gurobi.
4. Bring the dual bound and the moat bound into the repository.
