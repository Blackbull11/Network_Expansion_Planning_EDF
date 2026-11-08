# Scaling to TrèsGrand instances

> Written on 8 November 2026. **This is a proposal: nothing below is implemented yet.** Measured facts come from the instance files and past experiments. Anything labelled *estimate* or *hypothesis* has not been tested.

## Recommendation

Stop treating TrèsGrand as one problem. **A feeder is small and local, whatever the size of the instance**, so the instance can be solved as many small exact problems that are combined and improved.

The proposed method has three parts:

1. **Large neighbourhood search (LNS):** start from a heuristic network, then repeatedly re-optimize a small area of 2 or 3 feeders exactly, with the tools we already have (`subset_dp.py` and Gurobi). This produces valid networks of good quality.
2. **Feeder pool:** keep every valid feeder met during the search, and let Gurobi pick the best combination.
3. **Column generation lower bound:** use feeder prices, computed source by source, to measure how far the network is from the optimum.

The realistic target is a **valid network on every feasible instance, with a certified optimality gap**, not a proven optimum. A proof of optimality on 1,000 loads is out of reach of any method considered here.

## 1. Why the current methods cannot scale

| Method | What blocks it on TrèsGrand (300 to 1,000 loads, 2,000 to 4,000 nodes) |
|---|---|
| Gurobi PLNE | The model has `feeders × (5·edges + nodes + loads)` variables: **3.1 to 9.3 million** on TrèsGrand, and **6.7 to 33 million** with the tightening. The smallest instance today has about 7,000. In 20 minutes, Gurobi already proves only 1 of 8 Moyen_N150 instances ([Experiment 6](../experiments/experiment_6_gurobi_campaign_petit_moyen.md)). |
| Dynamic program | About 3^(loads) operations. It is limited to about 20 loads, and TrèsGrand has up to 1,000. |
| Dual bound as implemented | One column per allowed group of loads, so up to 2^(loads) per source. |
| Voronoi + Steiner | Fast, but it rarely satisfies the feeder limits: valid on **1 of 53** TrèsGrand instances tested (variable distance), 0 of 52 Grand, 5 of 54 Moyen. These counts predate the removal of the infeasible instances. |
| ToMATo decomposition | On `TrèsGrand_N2000_S15_0.25_Hiver` (later proven infeasible, because of one unreachable load), the bootstrap finds only 1 significant density peak, and the best configuration still has 10 violations. Its regions hold up to 392 loads. |

The last line matters: splitting by **load density** does not produce small regions, and there is no natural clustering to exploit. Any split has to follow **feeders**, not density.

## 2. The observation that makes scaling possible

Figures are measured on the 37 TrèsGrand instances not proven infeasible, except the nearest-load distance (sample of 8 instances, 150 loads each) and the Petit_N60 row.

**A feeder is bounded by the limits, not by the instance.**

| Quantity | Value | Consequence |
|---|---|---|
| Feeder power limit / mean load | 8 MW / 1.5 to 1.8 MW | About 5 loads of the same sign per feeder |
| Feeder length limit / nearest-load distance | 50 km / 2.4 to 3.7 km | A feeder cannot reach more than about 15 to 20 loads |
| Loads per feeder if all feeders are full | 6.1 to 8.3 | Same order |
| Optimal feeders on Petit_N60 (64 feeders) | mean 3.1 loads, at most 8 | Direct observation on proven optima |

So a feeder stays around 5 to 20 loads whether the instance has 300 or 1,000 loads.

**Loads have few choices.** Each load has only **4.3 to 5.1 sources** within 50 km, and each source reaches **73 to 363 loads**. Long-range interactions between distant parts of the network come only through limits that are mostly not tight.

**Exact solvers already handle this size.** The dynamic program takes 5 to 28 s on 15 to 16 loads, and Gurobi proved optima on Petit instances of up to 19 loads, stopping at the bound within minutes.

**Not all instances are easy.** Two have net demand at 93% of the feeder budget (`N3000_S12_0.20_Hiver`, `N4000_S12_0.15_Hiver`), so a valid assignment is hard to find. Two have gross production above the feeder budget (`N4000_S12_0.20_Ete` at 1.06, `N4000_S12_0.25_Ete` at 1.33), so producers and consumers **must** be mixed in the same feeders.

## 3. Proposed solution

```text
                instance
                   |
     (0) reduce the graph, screen infeasible instances
                   |
     (1) assign loads to sources with an LP (capacity, reach)
                   |
     (2) build feeders with the existing merge-tree heuristic
                   |
     (3) LNS: re-optimize small areas exactly, in parallel   <--+
                   |                                            |
     (4) feeder pool: Gurobi picks the best combination   ------+
                   |
     (5) column generation: valid lower bound and new feeders
                   |
         network + certified gap, judged by verify.py
```

### Step 0: reduce and screen

Contract degree-1 and degree-2 transit nodes (`TDA_algos/reductions.py` exists), and remove loads and lines that cannot be used within 50 km of any source. The two necessary conditions of Experiment 3 and the transportation LP in `TDA_algos/bounds.py` then reject provably infeasible instances early.

### Step 1: assign loads to sources

A linear program over the pairs (load, source within 50 km) gives each source a load respecting its feeder budget and capacity. This is what makes the tight instances tractable: the heuristics grow Voronoi cells and only discover overloaded sources afterwards.

### Step 2: build feeders

Use `TDA_algos/departs.py` (merge tree) inside each source's area. It is valid on 50 of 54 Petit instances in 0.2 s, against 23 of 54 for Voronoi variable distance ([Experiment 7](../experiments/experiment_7_tda_heuristics_and_plne.md)). It remains to be tested at scale.

### Step 3: large neighbourhood search

This is the core of the method.

```text
LNS(network):
    while time remains:
        pick a seed feeder: violated first, then long or recently unimproved
        area <- the seed feeder + its 1 or 2 nearest feeders       # 15 to 20 loads
        subproblem:
            the loads of area,
            the nodes of area + free nodes nearby,
            excluding the nodes of all other feeders,
            with the remaining source capacity and feeder count
        solve exactly: subset_dp if its feeders are disjoint, else Gurobi with a time limit
        if result is better (violations, excess, length): replace the feeders of area
        add every valid feeder seen to the pool
```

**Why it should work.** Every step uses a solver that is proven on this size. The acceptance rule already exists in `TDA_algos/local_search.py`. Areas that do not touch can run at the same time on all cores.

**Cost estimate.** A pass over 1,000 loads is about 67 areas of 15 loads. At 20 s each on 16 cores, that is about 1.5 minutes per pass. *Estimate from DP times on 15 to 16 loads. Areas that need Gurobi will be slower.*

### Step 4: feeder pool

The search visits many valid feeders that it then discards. A set-partitioning model picks the best combination of feeders seen so far:

```text
minimize    sum of cost[feeder] * λ[feeder]
subject to  each load is in exactly one chosen feeder
            each node is in at most one chosen feeder     # node-disjointness, which the DP ignores
            each source: at most nmax feeders, |power| <= capacity
            λ binary
```

The feeder LP (tested without the node rows) was integer on 15 of the 18 Petit_N60 instances, so this model is expected to be easy for Gurobi *(hypothesis at large scale)*. The chosen network can restart the LNS.

### Step 5: lower bound and new feeders by column generation

Lower bounds are what turn a good network into a *certified* one. Prices come from the same LP, with feeders generated on demand.

```text
COLUMN_GENERATION:
    master <- LP over the feeder pool
    repeat:
        solve master, read prices π (loads) and σ (nodes)
        for each source r, in parallel:
            look for feeders of r with sum of prices above their cost
        add them to master
    lower_bound <- sum of π + sum over r of nmax[r] * (best reduced cost of r)   # valid at every iteration
```

**The decomposition by source is the main benefit.** Given prices, each source is an independent problem on its own 73 to 363 reachable loads. So the work is split into 12 to 18 parallel tasks, regardless of instance size.

**The hard part is pricing:** finding the best feeder of a source is a Steiner tree problem with power and length limits. A valid lower bound needs an exact pricing, or at least a lower bound on it, which can come from an LP relaxation or a relaxed route. Heuristic pricing still produces good new feeders for the pool.

### A side benefit: real constraints become easy

The per-line power and P×L limits are not in `verify.py` or the PLNE, and the optimum of `S2_0.25_Ete` breaks both (Experiment 5). In a feeder-based method, a feeder is an explicit tree, so these limits are a simple check when a feeder is created. No new variables or big-M constraints are needed.

## 4. Alternatives considered

| Alternative | Decision | Reason |
|---|---|---|
| Bigger Gurobi model, or longer time limit | Rejected | 3 to 33 million variables. Its LP is weak even on small instances |
| Dynamic program with sampled loads | Rejected | Bounds are weak, since dropped loads hide most of the power and length pressure |
| Dynamic program only on nearby groups | Kept as pricing tool | Enumerating groups of at most about 8 nearby loads is the natural exact pricing to try |
| Density decomposition (ToMATo) | Rejected as the main approach | No natural clusters, and regions stay too large |
| Fix-and-optimize by source | Included inside LNS | Areas are smaller, and they can cross source boundaries |
| Persistence-based cuts in the PLNE (`TDA_algos/plne_cuts.py`) | Separate line of work | They strengthen the LP but leave the model size unchanged. Could help the Gurobi subproblems of the LNS |
| Full branch-and-price | Later, if the gap is too large | The same model as step 5 with branching added. It is the natural route to exact solutions if the pricing is fast |

## 5. What this can and cannot promise

- **Can be expected:** valid networks on most feasible instances, within an hour, with a quality that improves over time.
- **Hypothesis to test:** feeder-based bounds give a small gap, because the feeder LP was tight on 14 of 18 small instances. That may not hold when many feeders interact.
- **Main risks:** exact pricing may be too slow, and instances near the feeder budget may leave the search stuck without a valid network.
- **Not achievable:** a proof of optimality on instances with hundreds of loads, in general.

## 6. Roadmap

Each stage has a measurable outcome, so a failure stops the work early.

| Stage | Work | Test | Success criterion |
|---|---|---|---|
| 1 | LNS with `subset_dp` and Gurobi subproblems | Petit_N60 and N75, where the optimum is known | Within 1% of the optimum |
| 2 | Feasibility phase (steps 0 to 2) | Moyen and Grand | Valid network on most instances, compared with 5 of 54 and 0 of 52 today |
| 3 | Feeder pool (step 4) | Same instances | Length reduced compared with stage 2 alone |
| 4 | Scaling to TrèsGrand with parallel areas | The 37 TrèsGrand instances | Valid networks within 1 hour on 16 cores |
| 5 | Column generation lower bound | Petit first (known optima), then TrèsGrand | Certified gap reported on each instance |

Stages 1 and 2 use existing components only and already answer whether the local approach works.
