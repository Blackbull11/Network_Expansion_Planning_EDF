# Experiment 3: instances that fail the naive feasibility conditions

> Run on 8 November 2026 on all 216 instances of `data/verified_instances/`. **26 instances are provably infeasible** and were moved to `data/verified_instances/unfeasible/`. The other 190 are not proven solvable: they only pass these two necessary conditions.

## What was tested

All 216 instances share the same limits: 70 MW and at most 8 feeders per source, 8 MW and 50 km per feeder.

**Condition 1, power.** |sum of the demands| must not exceed the sum over the sources of min(capacity, number of feeders × feeder power). Here that is min(70, 8 × 8) = **64 MW per source**, so the feeders are the binding limit, not the 70 MW source capacity.

**Condition 2, distance.** Every Conso and Prod node must be within 50 km of a source, along the shortest path in the graph. A feeder is a tree whose total length is at most 50 km, so the path from its source to any of its loads cannot be longer than that.

Failing either condition proves that no solution exists. The checks were computed directly from the JSON files and agree with `necessary_conditions` in `verify.py` on all 216 instances. Conso and Prod totals and feeder capacities match `ratios_charge_nette.csv`.

## Results

| Class | Instances | Fail power | Fail distance | Infeasible | Remaining |
|---|---|---|---|---|---|
| Petit | 54 | 0 | 0 | 0 | 54 |
| Moyen | 54 | 0 | 4 | 4 | 50 |
| Grand | 54 | 1 | 4 | 5 | 49 |
| TrèsGrand | 54 | 6 | 11 | 17 | 37 |
| **Total** | **216** | **7** | **19** | **26** | **190** |

No instance fails both conditions.

## Power failures (7 instances)

All are winter instances, which have no production to offset consumption. "Source limit only" shows whether the instance would have passed if only the 70 MW source capacity were checked.

| Instance | Net demand (MW) | Feeder limit (MW) | Excess (MW) | Source limit (MW) | Source limit only |
|---|---|---|---|---|---|
| `Grand_N700_S3_0.25_Hiver` | 207 | 192 | 15 | 210 | Passes |
| `TrèsGrand_N3000_S12_0.25_Hiver` | 888 | 768 | 120 | 840 | Fails |
| `TrèsGrand_N4000_S12_0.20_Hiver` | 915 | 768 | 147 | 840 | Fails |
| `TrèsGrand_N4000_S12_0.25_Hiver` | 1197 | 768 | 429 | 840 | Fails |
| `TrèsGrand_N4000_S15_0.20_Hiver` | 991 | 960 | 31 | 1050 | Passes |
| `TrèsGrand_N4000_S15_0.25_Hiver` | 1188 | 960 | 228 | 1050 | Fails |
| `TrèsGrand_N4000_S18_0.25_Hiver` | 1185 | 1152 | 33 | 1260 | Passes |

## Distance failures (19 instances)

| Instance | Loads beyond 50 km | Loads in total | Farthest load (km) |
|---|---|---|---|
| `Moyen_N150_S6_0.25_Hiver` | 1 | 38 | 50.14 |
| `Moyen_N200_S4_0.15_Ete` | 2 | 30 | 53.27 |
| `Moyen_N200_S4_0.15_Hiver` | 2 | 30 | 53.51 |
| `Moyen_N200_S4_0.25_Hiver` | 1 | 50 | 52.65 |
| `Grand_N500_S9_0.15_Ete` | 1 | 76 | 53.44 |
| `Grand_N500_S9_0.20_Ete` | 2 | 100 | 52.72 |
| `Grand_N500_S9_0.25_Ete` | 1 | 124 | 50.25 |
| `Grand_N500_S9_0.25_Hiver` | 1 | 125 | 56.66 |
| `TrèsGrand_N2000_S15_0.15_Ete` | 1 | 300 | 52.29 |
| `TrèsGrand_N2000_S15_0.20_Ete` | 1 | 400 | 51.01 |
| `TrèsGrand_N2000_S15_0.20_Hiver` | 1 | 400 | 52.14 |
| `TrèsGrand_N2000_S15_0.25_Hiver` | 1 | 500 | 52.37 |
| `TrèsGrand_N2000_S18_0.20_Ete` | 1 | 400 | 51.76 |
| `TrèsGrand_N3000_S15_0.15_Ete` | 1 | 450 | 54.10 |
| `TrèsGrand_N3000_S15_0.15_Hiver` | 5 | 450 | 55.57 |
| `TrèsGrand_N3000_S15_0.20_Ete` | 2 | 600 | 57.91 |
| `TrèsGrand_N3000_S15_0.20_Hiver` | 7 | 600 | 57.01 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 8 | 750 | 56.65 |
| `TrèsGrand_N3000_S15_0.25_Hiver` | 6 | 750 | 55.41 |

### Every load beyond 50 km (45 nodes)

"Node" is `id_noeud` and "Nearest source" is the `id_noeud` of the closest source.

| Instance | Node | Type | Demand (MW) | Nearest source | Distance (km) |
|---|---|---|---|---|---|
| `Moyen_N150_S6_0.25_Hiver` | 14 | Conso | 2 | 3 | 50.14 |
| `Moyen_N200_S4_0.15_Ete` | 25 | Prod | −1 | 3 | 53.27 |
| `Moyen_N200_S4_0.15_Ete` | 31 | Prod | −1 | 3 | 51.00 |
| `Moyen_N200_S4_0.15_Hiver` | 33 | Conso | 1 | 1 | 53.51 |
| `Moyen_N200_S4_0.15_Hiver` | 28 | Conso | 1 | 3 | 51.74 |
| `Moyen_N200_S4_0.25_Hiver` | 29 | Conso | 2 | 3 | 52.65 |
| `Grand_N500_S9_0.15_Ete` | 54 | Prod | −3 | 4 | 53.44 |
| `Grand_N500_S9_0.20_Ete` | 21 | Conso | 2 | 6 | 52.72 |
| `Grand_N500_S9_0.20_Ete` | 58 | Conso | 1 | 6 | 50.96 |
| `Grand_N500_S9_0.25_Ete` | 70 | Conso | 2 | 6 | 50.25 |
| `Grand_N500_S9_0.25_Hiver` | 61 | Conso | 1 | 4 | 56.66 |
| `TrèsGrand_N2000_S15_0.15_Ete` | 100 | Conso | 1 | 6 | 52.29 |
| `TrèsGrand_N2000_S15_0.20_Ete` | 154 | Conso | 2 | 6 | 51.01 |
| `TrèsGrand_N2000_S15_0.20_Hiver` | 340 | Conso | 0 | 6 | 52.14 |
| `TrèsGrand_N2000_S15_0.25_Hiver` | 474 | Conso | 2 | 6 | 52.37 |
| `TrèsGrand_N2000_S18_0.20_Ete` | 230 | Prod | −1 | 14 | 51.76 |
| `TrèsGrand_N3000_S15_0.15_Ete` | 290 | Prod | −1 | 6 | 54.10 |
| `TrèsGrand_N3000_S15_0.15_Hiver` | 317 | Conso | 2 | 6 | 55.57 |
| `TrèsGrand_N3000_S15_0.15_Hiver` | 394 | Conso | 2 | 6 | 54.17 |
| `TrèsGrand_N3000_S15_0.15_Hiver` | 338 | Conso | 0 | 14 | 53.81 |
| `TrèsGrand_N3000_S15_0.15_Hiver` | 330 | Conso | 2 | 6 | 52.82 |
| `TrèsGrand_N3000_S15_0.15_Hiver` | 332 | Conso | 1 | 6 | 50.93 |
| `TrèsGrand_N3000_S15_0.20_Ete` | 337 | Prod | −3 | 6 | 57.91 |
| `TrèsGrand_N3000_S15_0.20_Ete` | 205 | Conso | 2 | 14 | 52.46 |
| `TrèsGrand_N3000_S15_0.20_Hiver` | 312 | Conso | 2 | 14 | 57.01 |
| `TrèsGrand_N3000_S15_0.20_Hiver` | 119 | Conso | 0 | 14 | 55.47 |
| `TrèsGrand_N3000_S15_0.20_Hiver` | 110 | Conso | 1 | 14 | 53.99 |
| `TrèsGrand_N3000_S15_0.20_Hiver` | 452 | Conso | 2 | 6 | 53.86 |
| `TrèsGrand_N3000_S15_0.20_Hiver` | 456 | Conso | 0 | 6 | 51.48 |
| `TrèsGrand_N3000_S15_0.20_Hiver` | 75 | Conso | 1 | 6 | 51.41 |
| `TrèsGrand_N3000_S15_0.20_Hiver` | 281 | Conso | 2 | 14 | 50.30 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 402 | Prod | −2 | 6 | 56.65 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 745 | Prod | −2 | 14 | 54.06 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 376 | Conso | 2 | 6 | 52.65 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 664 | Prod | −2 | 14 | 51.58 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 514 | Prod | −2 | 14 | 51.19 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 660 | Prod | −3 | 14 | 50.98 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 45 | Conso | 0 | 14 | 50.79 |
| `TrèsGrand_N3000_S15_0.25_Ete` | 251 | Conso | 2 | 14 | 50.55 |
| `TrèsGrand_N3000_S15_0.25_Hiver` | 407 | Conso | 2 | 6 | 55.41 |
| `TrèsGrand_N3000_S15_0.25_Hiver` | 253 | Conso | 2 | 14 | 55.19 |
| `TrèsGrand_N3000_S15_0.25_Hiver` | 596 | Conso | 1 | 6 | 54.99 |
| `TrèsGrand_N3000_S15_0.25_Hiver` | 98 | Conso | 1 | 6 | 53.93 |
| `TrèsGrand_N3000_S15_0.25_Hiver` | 459 | Conso | 1 | 14 | 52.85 |
| `TrèsGrand_N3000_S15_0.25_Hiver` | 635 | Conso | 1 | 14 | 50.89 |

## Limits of this test

- These conditions are only necessary, so the 190 remaining instances may still be infeasible.
- `TrèsGrand_N2000_S15_0.20_Hiver` fails only because of a Conso node with 0 MW demand. `verify.py` requires every Conso to be connected, so it counts as infeasible. It would pass if 0 MW nodes did not need connecting.
- Transit nodes do not have to be served and were not tested. 14 instances that pass have Transit nodes beyond 50 km.
- The P×L limit is not part of these conditions.
- The demand maps in `data/verified_instances/demand_maps/` and the rows of `ratios_charge_nette.csv` still cover all 216 instances.

## Reproduce

From `projetX/projetX`, this prints the violated conditions of each moved instance:

```powershell
$env:PYTHONUTF8 = 1
& C:\Users\siand\dev\X_3A\RO_EDF\.venv\Scripts\python.exe -c "from pathlib import Path; from src.instances import load_instance; from verify import necessary_conditions; [print(p.stem, necessary_conditions(load_instance(p))) for p in sorted(Path('data/verified_instances/unfeasible').glob('*.json'))]"
```
