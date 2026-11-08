# Network Expansion Planning (EDF)

Design a radial electricity distribution network at minimum total line length. Each source feeds a set of feeders (*départs*), and every load must be connected without breaking the feeder and source limits.

| Limit, on the verified instances | Value |
|---|---|
| Feeder length | 50 km |
| Feeder power | 8 MW |
| Source capacity | 70 MW |
| Feeders per source | 8 |

The objective is **J = total length of the network (km)**. `verify.py` is the only judge of feasibility: every algorithm is checked by it, whatever built the solution.

For the project history, results and open problems, see [CLAUDE.md](CLAUDE.md). For the short results summary, see [report.md](projetX/projetX/report.md).

## Setup

All commands run from `projetX/projetX`.

```powershell
cd projetX\projetX
..\..\..\.venv\Scripts\Activate.ps1   # the existing environment, next to the repository
$env:PYTHONUTF8 = "1"                 # required: instance names contain "è"
```

To build a new environment instead (tested on Python 3.14):

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install networkx geopandas pulp highspy gurobipy matplotlib scipy
```

**Gurobi needs a full licence for the verified instances.** The licence bundled with `pip install gurobipy` stops at 2,000 variables and 2,000 constraints, and the smallest verified instance has about 7,000 variables. The bundled licence is enough for `plne_gurobi.py demo --n 10`.

## Instances

| Set | Location | Content |
|---|---|---|
| Verified, usable | `data/verified_instances/*.json` | 190 instances that pass the necessary feasibility conditions (not proven solvable) |
| Verified, infeasible | `data/verified_instances/unfeasible/` | 26 instances proven infeasible |
| Generated | `data/instances/` | g01 to g40. Not in the repository: run `python data/generate_instances.py` first |
| Worked example | `data/exemple_ete_4_PS/` | Empty folder at the moment |

Verified names follow `<Class>_N<nodes>_S<sources>_<ratio>_<season>`, for example `Petit_N60_S1_0.25_Hiver`.

| Class | Nodes | Loads | Usable instances |
|---|---|---|---|
| Petit | 60 to 90 | 8 to 22 | 54 |
| Moyen | 150 to 250 | 22 to 62 | 50 |
| Grand | 500 to 700 | 75 to 176 | 49 |
| TrèsGrand | 2,000 to 4,000 | 300 to 1,000 | 37 |

## Which method to use

| Goal | Method | Reach |
|---|---|---|
| Prove an instance infeasible | `verify.necessary_conditions` | All sizes. An empty result proves nothing |
| Fast valid network | Merge-tree feeders (`run_tda_plne.py heuristics`) | Valid on 50 of 54 Petit, about 0.2 s on average. Not tested on larger classes in this run |
| Exact optimum or strong lower bound | `subset_dp.py`, alone or inside Gurobi with `--dp-stop` | Up to about 20 loads |
| Proven optimum | `plne_gurobi.py verified --tighten` | Most of Petit. Only 1 of 8 Moyen_N150 instances was proven within the 20-minute limit |
| Cheap lower bound | `TDA_algos/bounds.py` | All sizes, but weak (53 to 60% of the optimum) |
| TrèsGrand | None gives valid networks yet | Proposal in [docs/scaling.md](projetX/projetX/docs/scaling.md) |

## Running the algorithms

### Exact methods

```powershell
# Gurobi, one instance or a whole class
python plne_gurobi.py verified --prefix Petit_N60_S1_0.25_Hiver --log
python plne_gurobi.py verified --prefix Petit_N60 --tighten --dp-stop --no-warm-start --threads 4 --time-limit 1200

# Gurobi on a small random instance (works with the bundled licence)
python plne_gurobi.py demo --n 10 --seed 5 --lmax-mode shortest --disjoint

# HiGHS through PuLP
python plne_reseau.py --n 30 --seed 3 --disjoint
python run_plne_verified.py --single data/verified_instances/Petit_N60_S1_0.25_Hiver.json --time-limit 600
```

Main `plne_gurobi.py verified` options:

| Option | Effect |
|---|---|
| `--prefix` | Start of the instance name: a class or one instance |
| `--max-nodes`, `--n` | Keep only small instances, or the n smallest |
| `--tighten` | Stronger model with valid inequalities |
| `--dp-stop` | Compute the `subset_dp` bound (up to 20 loads) and stop as soon as Gurobi reaches it |
| `--no-warm-start` | Do not start from the heuristic solution |
| `--time-limit`, `--threads`, `--max-feeders`, `--log` | Solver controls |

The dynamic program has no command line. Call it from Python:

```powershell
python -c "from src.instances import load_instance; from subset_dp import lower_bound; print(lower_bound(load_instance('data/verified_instances/Petit_N60_S1_0.25_Hiver.json')))"
```

It returns `bound`, `n_loads` and `time`. The cost is about 3^(loads), so stay at or below about 20 loads.

### Heuristics

```powershell
# Voronoi + Steiner, 4 distances, on the verified instances
python run_voronoi_verified.py Petit_ Moyen_ --no-figures --out output/voronoi_verified_clean
python run_voronoi_verified.py --report      # rebuild the report only

# Same on the generated instances (generate them first)
python main.py --petites --label petites_variable
python main.py g05 --show

# Every TDA heuristic on the Petit instances
python run_tda_plne.py heuristics
```

Use `--out` with `run_voronoi_verified.py`: the results already in `output/voronoi_verified/` predate the instance cleanup, and `--restart` would erase them.

### TDA experiments

```powershell
python TDA.py e7 g05 g10 --label test                  # generated instances
python TDA.py grid --petites
python TDA.py e10 --labels output/tda/<run>_grid/results.csv
```

| Experiment | Compares |
|---|---|
| `e1` | Persistence diagrams as instance descriptors |
| `e2` | Necessary conditions: `verify.py` against `TDA_algos/bounds.py` |
| `e3` | H0 lower bound against heuristics and PLNE |
| `e4` | Exact graph reductions (`--plne T` to time the PLNE with and without) |
| `e5` | Variable distance from fixed weights against ToMATo initialisation |
| `e6` | Voronoi on lengths against spectral costs and spectral cells |
| `e7` | Steiner tree against designed feeders (Kruskal, savings, merge tree) |
| `e8` | Local search after the heuristic |
| `e9` | Shrinking only against shrinking and growing |
| `e10` | Predict the best configuration from features (`--labels` required) |
| `e11` | Time of each TDA component |
| `grid` | Full ablation of init, distance, feeders and local search |
| `all` | `e1` to `e9`, `grid`, then `e10` |

Select instances with `--petites`, `--grandes` or name motifs. Results go to `output/tda/<date>_<experiment>[_label]/`. `--set verified` does not work yet, see [Known problems](#known-problems).

### Other tools

```powershell
python run_tda_plne.py exact --work-limit 120 --workers 2   # 5 Gurobi variants on Petit_N60
python run_tda_plne.py report                               # write the CSV files
python run_gurobi_campaign.py --workers 3 --threads 5 --time-limit 1200 --seed 2026
python run_tomato_decomposition.py Grand_N500_S3_0.25_Ete --workers 12 --ls-time 20
python plot_demands.py Petit_ Moyen_                        # demand maps, in data/verified_instances/demand_maps/
```

`run_tomato_decomposition.py` needs a name motif that matches exactly one instance. Its default instance is now in `unfeasible/`, so the command fails without an argument.

## Reproducing the experiments

| Experiment | Command | Notes |
|---|---|---|
| [1. Voronoi on 40 graphs](projetX/projetX/experiments/experiment_1.md) | `python data/generate_instances.py`, then `python main.py --petites --label petites_variable` and `python main.py --grandes --label grandes_variable` | |
| [2. Gurobi vs HiGHS](projetX/projetX/experiments/experiment_2_plne_gurobi_vs_highs.md) | `python plne_gurobi.py verified --prefix Petit_N60_S1_0.25_Hiver`, then the `run_plne_verified.py --single` command above | |
| [3. Infeasible instances](projetX/projetX/experiments/experiment_3_unfeasible_instances.md) | Call `verify.necessary_conditions(load_instance(path))` on each file | The sweep itself was a temporary script |
| [4. Lower bounds](projetX/projetX/experiments/experiment_4_lower_bounds.md) | `TDA_algos.bounds.lower_bound(graph)` and `subset_dp.lower_bound(graph)` | The moat bound script is not in the repository |
| [5. Gurobi vs bounds, Petit_N60](projetX/projetX/experiments/experiment_5_gurobi_lower_bounds_petit_n60.md) | `python plne_gurobi.py verified --prefix Petit_N60 --tighten --dp-stop --no-warm-start --threads 4 --time-limit 1200` | `--dp-stop` ends earlier than the reported run. The dual bound is not in the repository |
| [6. Gurobi campaign, Petit and Moyen](projetX/projetX/experiments/experiment_6_gurobi_campaign_petit_moyen.md) | `python run_gurobi_campaign.py --workers 3 --threads 5 --time-limit 1200 --seed 2026` | Skips the runs already saved |
| [7. TDA heuristics and PLNE](projetX/projetX/experiments/experiment_7_tda_heuristics_and_plne.md) | `python run_tda_plne.py heuristics`, `exact --work-limit 120 --workers 2`, then `report` | |

## Results on disk

| Script | Output in `output/` | On a rerun |
|---|---|---|
| `plne_gurobi.py verified` | `plne_gurobi/<instance>.json`, `.png`, `.log`, `summary.csv` | Replaces |
| `run_plne_verified.py` | `plne_verified/<instance>.json`, `summary.csv` | Replaces |
| `run_voronoi_verified.py` | `voronoi_verified/` or `--out`: `results.jsonl`, `results.csv`, `report.md`, `charts/`, `figures/` | Resumes |
| `run_gurobi_campaign.py` | `gurobi_campaign/progress.log`, `summary.csv`, `runs/` | Resumes |
| `run_tda_plne.py` | `tda_plne/heuristics.jsonl`, `exact.jsonl`, CSV after `report` | Resumes |
| `run_tomato_decomposition.py` | `tomato_decomposition/<instance>/` with `report.md`, `results.csv`, `solution.json`, figures | Replaces |
| `main.py` | `recap_<stamp>.md`, `.csv`, `details_<stamp>.txt`, `charts/`, `figures/` | New stamp, figures replaced |
| `TDA.py` | `tda/<date>_<experiment>[_label]/` | New folder |

## Known problems

- **Wrong path for the verified instances.** `src/instances.py` looks in `verified_instances/` and the files are in `data/verified_instances/`. Until `VERIFIED_DIR` is fixed, `main.py --verified`, `TDA.py --set verified` and the batch mode of `run_plne_verified.py` find no instance. `run_plne_verified.py --single <path>` is not affected, and the other runners read the correct folder.
- **`verify.py` and the PLNE do not enforce everything in the instance files.** Power is checked only at the head of each feeder, and P×L is not checked. So the reported optima are optimal for the repository's model, which can differ from the file's full limits (see Experiment 5).
