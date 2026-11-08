# Experiment 6: Gurobi with early stop on larger Petit and Moyen instances

> Run on 8 November 2026 with [run_gurobi_campaign.py](../run_gurobi_campaign.py). Raw values: [summary.csv](../output/gurobi_campaign/summary.csv), run by run in [runs/](../output/gurobi_campaign/runs/), timeline in [progress.log](../output/gurobi_campaign/progress.log).

## What was tested

**Question:** how far does Gurobi prove optimality when it can stop at the dynamic-program (DP) bound?

**Instances:** 5 Petit_N75 and 5 Petit_N90 drawn at random (seed 2026), then Moyen instances from the fewest loads upwards, until 5 of them reach a limit with no proven optimum.

**Setup:** same tightened model as [Experiment 5](experiment_5_gurobi_lower_bounds_petit_n60.md), no starting solution, 1,200 s of Gurobi time per run, MIPGap 10⁻⁶. 3 runs at a time with 5 threads each, 3 GB memory limit per run.

**Stop rule:** Gurobi stops as soon as it finds a network at the bound from [subset_dp.py](../subset_dp.py), computed up to 20 loads. A valid network at a lower bound is optimal, and `verify.py` must accept it (status `OPTIMAL_BY_DP`). Above 20 loads there is no DP bound, so Gurobi must close its own gap. That covers one Petit_N90 instance and all Moyen instances, so a Moyen run counts as a timeout when it hits a limit unproven.

## Results

The optimality gap uses the best lower bound available (DP or Gurobi). "Time" is the Gurobi time. No run hit the memory limit.

| Instance | Loads | Sources | Status | J (km) | DP bound | Gurobi bound | Gap | DP time (s) | Time (s) |
|---|---|---|---|---|---|---|---|---|---|
| Petit_N75_S1_0.20_Hiver | 15 | 1 | Time limit | 89.11 | 84.59 | 81.06 | 5.1% | 6 | 1,200 |
| Petit_N75_S2_0.25_Ete | 18 | 2 | Optimal (Gurobi) | 114.70 | 114.63 | 114.70 | 0% | 162 | 11 |
| Petit_N75_S3_0.20_Ete | 16 | 3 | Optimal at DP bound | 190.23 | 190.23 | | 0% | 28 | 387 |
| Petit_N75_S3_0.20_Hiver | 15 | 3 | Optimal at DP bound | 149.12 | 149.12 | | 0% | 14 | 100 |
| Petit_N75_S3_0.25_Hiver | 19 | 3 | Optimal at DP bound | 186.23 | 186.23 | | 0% | 584 | 46 |
| Petit_N90_S3_0.15_Ete | 14 | 3 | Optimal at DP bound | 125.01 | 125.01 | | 0% | 5 | 752 |
| Petit_N90_S3_0.15_Hiver | 14 | 3 | Optimal at DP bound | 122.21 | 122.21 | | 0% | 6 | 11 |
| Petit_N90_S3_0.20_Ete | 18 | 3 | Optimal at DP bound | 167.96 | 167.96 | | 0% | 246 | 302 |
| Petit_N90_S3_0.20_Hiver | 18 | 3 | Optimal at DP bound | 158.40 | 158.40 | | 0% | 239 | 30 |
| Petit_N90_S3_0.25_Hiver | 22 | 3 | Time limit | 180.22 | none | 171.10 | 5.1% | | 1,200 |
| Moyen_N150_S6_0.15_Ete | 22 | 6 | Optimal (Gurobi) | 173.10 | none | 173.10 | 0% | | 108 |
| Moyen_N150_S4_0.15_Hiver | 22 | 4 | Time limit | 186.52 | none | 183.06 | 1.9% | | 1,200 |
| Moyen_N150_S4_0.15_Ete | 22 | 4 | Time limit | 198.96 | none | 193.07 | 3.0% | | 1,200 |
| Moyen_N150_S2_0.15_Ete | 22 | 2 | Time limit | 164.31 | none | 150.17 | 8.6% | | 1,200 |
| Moyen_N150_S6_0.15_Hiver | 22 | 6 | Time limit | 225.04 | none | 204.44 | 9.2% | | 1,200 |
| Moyen_N150_S2_0.15_Hiver | 22 | 2 | Time limit | 167.81 | none | 140.95 | 16.0% | | 1,200 |
| Moyen_N150_S2_0.20_Hiver | 30 | 2 | Time limit | 217.90 | none | 164.88 | 24.3% | | 1,200 |
| Moyen_N150_S2_0.20_Ete | 30 | 2 | Time limit, **no network found** | | none | 175.72 | | | 1,200 |

Every network found was accepted by `verify.py`.

| Set | Runs | Proven optimal | Unproven at the limit |
|---|---|---|---|
| Petit_N75 and N90 | 10 | 8 (7 stopped at the DP bound, 1 by Gurobi) | 2 |
| Moyen_N150 | 8 | 1 | 7 |

The rule triggered at the 6th Moyen run. The two runs still in progress then also timed out, giving 7 timeouts. The other 42 Moyen instances were not run. The campaign took 1 h 24 min.

## Conclusions

- **The early stop works.** 7 of the 8 Petit runs proven optimal stopped at the DP bound, without Gurobi closing its own gap.
- **The DP bound equals the optimum on 7 of the 9 Petit runs that have one.** On `Petit_N75_S2_0.25_Ete` it is 0.07 km below the proven optimum. On `Petit_N75_S1_0.20_Hiver` the optimum is unknown, and the DP bound (84.59 km) is tighter than Gurobi's (81.06 km).
- **The limit of the exact approach is reached at the smallest Moyen instances.** Within 20 minutes, Gurobi proves 1 of 6 instances at 22 loads and leaves gaps of 1.9 to 16%. At 30 loads, the gap is 24% or no network is found at all.
- **The DP is the bottleneck near its cap,** at 584 s for 19 loads, longer than Gurobi needed after it (46 s).

## Limits

- This is a small sample. The Petit instances are random, and 3 of the 5 Petit_N90 draws share a size class.
- Above 20 loads, a timeout says nothing about how far the network is from the optimum beyond Gurobi's own gap.
- Without a start from a heuristic, Gurobi can end with no network, as on `Moyen_N150_S2_0.20_Ete`.
- As in Experiment 5, `verify.py` and the PLNE do not check the per-line power or P×L limits.

## Reproduce

From `projetX/projetX`, with the interpreter of the project's virtual environment and `PYTHONUTF8=1`:

```powershell
python run_gurobi_campaign.py --workers 3 --threads 5 --time-limit 1200 --seed 2026
```

The script skips runs already saved in `output/gurobi_campaign/runs/`, so an interrupted campaign resumes.
