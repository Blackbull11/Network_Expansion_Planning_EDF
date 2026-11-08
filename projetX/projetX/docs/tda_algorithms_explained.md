# Topological methods for distribution network planning

> **Summary of the results (§5).** Merge-tree feeders are the best fast heuristic on Petit, and persistence-based cuts with a TDA warm start help Gurobi prove optimality. Density clustering does not give a useful decomposition.

## 1. The problem

We are given a graph of candidate lines with lengths in km. Each node is a **source** (substation), a **load** (consumer `D > 0`, producer `D < 0`, or a consumer of 0 MW) or a **transit** node. We choose lines of minimal total length `J`, so that every load is connected to exactly one source and the network is radial: a forest of node-disjoint trees, one per source, none passing through a second source.

A **feeder** (*départ*) is the subtree hanging from one neighbour of a source.

| Constraint | Scope | Value on the test instances |
|---|---|---|
| Total length `L_départ` | each feeder | ≤ 50 km |
| Absolute net power `\|P_départ\|` | each feeder | ≤ 8 MW |
| Number of feeders | each source | ≤ 8 |
| Absolute net power `\|P_source\|` | each source | ≤ 70 MW |

Every solution in this document is judged by one independent checker, `verify.py`, whatever algorithm produced it. The exact formulation (PLNE), the dynamic program and the dual bound are explained in [exact_methods_explained.md](exact_methods_explained.md).

**Why the problem is hard.** Without the feeder limits it is a Steiner tree problem with the sources contracted, already NP-hard. The limits on length and power add a capacitated partition of the loads into feeders on top.

**Notation.** `G = (V, E)` is the graph, `n = |V|`, `m = |E|`. `d(u, v)` is the shortest-path distance. `S` is the set of sources, `C` the set of loads, `L_max = 50`, `P_max = 8`.

## 2. Topological data analysis in one page

Topological data analysis studies a shape **across all scales at once**, instead of picking one threshold.

- A **filtration** is a growing family of subgraphs indexed by a scale `t`. Example: keep the nodes where a density is at least `t`, and let `t` decrease.
- **H0** tracks connected components. A component is *born* at some `t` and *dies* when it merges with an older one (the **elder rule**).
- **H1** tracks independent cycles: a loop is born when it closes, and dies when it is filled in.
- The **persistence diagram** is the set of points (birth, death). The **persistence** or **prominence** of a feature is `|birth − death|`. Long-lived features are structure, short-lived ones are noise. Diagrams are stable: a small perturbation of the input moves the points by a small amount.
- A **merge tree** keeps the full history of which components merge with which, not only their lifetimes.

Four filtrations carry most of what we use:

| Filtration | What it sees in a network | Used by |
|---|---|---|
| Distance between loads (H0) | clusters of loads at every scale, length of the MST | bounds, cuts, capacitated clustering of feeders |
| Superlevel sets of a load density (H0) | towns and villages, valleys between them | ToMATo, decomposition, cuts |
| Superlevel sets of the distance to a source (H0) | corridors of the network as seen from the source | merge tree, feeder design |
| Edge length plus faces of the drawing (H1) | holes of the mesh, real alternative routes | filter of the local search |

## 3. Map of the algorithms

The code is in [TDA_algos/](../TDA_algos/). Each algorithm plays one of six roles.

```text
                   instance
                      │
     ┌────────────────┼──────────────────────────────┐
     ▼                ▼                              ▼
 DIAGNOSE          HEURISTIC                       EXACT (PLNE)
 conditions        partition into cells            warm start  ◄── heuristic
 lower bounds      (Voronoi, spectral, ToMATo)     cuts        ◄── H0 / ToMATo clusters
 descriptors          │                            restricted  ◄── heuristic edges
                      ▼                            edge set
                   design the feeders
                   (Steiner, merge tree, savings)
                      │
                      ▼
                   local search (H1, margin)
                      │
                      ▼
                   verify.py
```

| § | Algorithm | Module | Role | Objective |
|---|---|---|---|---|
| 4.1 | Union-find | `union_find` | tool | Maintain merging components |
| 4.2 | H0 of a vertex function | `persistence` | tool | Components of sub- or superlevel sets |
| 4.3 | H0 of the terminals | `persistence` | tool, bound | Single-linkage of the loads in near-linear time |
| 4.4 | H1 of the planarised graph | `persistence` | filter | Find the real holes of the mesh |
| 4.5 | Diagram summaries and distances | `vectorize` | descriptors | Turn diagrams into numbers, compare them |
| 4.6 | ToMATo | `clustering` | partition | Cluster the loads by density peaks |
| 4.7 | Split tree | `merge_tree` | feeder design | Corridor structure inside a cell |
| 4.8 | Bounds and necessary conditions | `bounds` | diagnose | Prove infeasibility, bound `J` from below |
| 4.9 | Spectral tools | `spectral` | partition | Cut the graph at its bottlenecks |
| 4.10 | Exact reductions | `reductions` | preprocessing | Remove useless transit nodes |
| 4.11 | Feeder design | `departs` | heuristic | Build feeders that respect L and P |
| 4.12 | Variable-distance Voronoi | `partition` | heuristic | Balance the cells of the sources |
| 4.13 | Local search | `local_search` | heuristic | Improve a solution with filtered moves |
| 4.14 | Instance features, k-NN | `features` | selection | Predict the best configuration |
| 4.15 | ToMATo decomposition | `decomposition` | decomposition | Split a large instance into subproblems |
| 4.16 | Persistence-based cuts | `plne_cuts` | exact | Strengthen the LP relaxation of the PLNE |

## 4. The algorithms

### 4.1 Union-find

- **Objective:** maintain a partition under merges. Every filtration below uses it.
- **Input:** items. **Operations:** `add`, `find`, `union`, `groups`.
- **Output:** the root of a set, or all the sets.
- **Cost:** almost constant per operation (path halving and union by size).

```text
FIND(x):      while parent[x] != x:  parent[x] <- parent[parent[x]];  x <- parent[x]
              return x
UNION(a, b):  ra, rb <- FIND(a), FIND(b);  if ra == rb: return None
              attach the smaller tree under the larger;  return the new root
```

### 4.2 H0 of a vertex function

- **Objective:** the persistence of the connected components of `{v : f(v) ≤ t}` (`h0_lower_star`), or of `{v : f(v) ≥ t}` (`h0_superlevel`, run on `−f`).
- **Input:** graph `G`, function `f` on the nodes.
- **Output:** diagram of the pairs (birth, death), and the node where each pair is born and dies. Essential classes have death `+∞`.
- **Cost:** `O(n log n + m α(n))`.

```text
sort the nodes by increasing f
for v in that order:
    new component {v}, born at v
    for each neighbour u of v already processed:
        if u and v are in two different components A and B:
            young <- the one born at the higher value            # elder rule
            if f(birth of young) < f(v): record the pair (f(birth of young), f(v))
            merge A and B, keep the older birth
for each remaining component: record (f(birth), +inf)
```

### 4.3 H0 of the terminals (`terminal_h0`)

- **Objective:** the single-linkage hierarchy of the loads under the shortest-path metric, without computing all pairwise distances.
- **Input:** `G`, terminals, optional groups born already merged (for example all the sources as one point).
- **Output:** the death times (the edges of the MST of the metric closure), their sum, the MST edges, the nearest-terminal distance of each terminal, and the Voronoi region of each node.
- **Cost:** `O(m log n)`, by Mehlhorn's theorem: an MST of the graph `G1'` below is an MST of the metric closure.

```text
multi-source Dijkstra from all terminals:  dist[v], region[v]
for each edge (u, v) with region[u] != region[v]:
    candidate edge region[u] - region[v], of length dist[u] + w(u, v) + dist[v]
Kruskal on the candidates, by increasing length
the accepted lengths are the death times of H0
```

### 4.4 H1 of the planarised graph (`h1_planar`)

- **Objective:** tell the holes of the network, which are real alternative routes, from the small triangles of the mesh.
- **Input:** `G` with node positions.
- **Output:** the diagram of the cycles, and for each edge the persistence of the largest cycle it closes.
- **Why planarise:** the drawings are not planar (source links cross other edges), so a virtual vertex is added at each crossing.
- **Filtration:** an edge enters at its length, a face at max(half its perimeter, its longest edge). Persistence ≈ half perimeter of the hole minus the length of the edge closing it.

```text
planarise G;  compute the faces of the drawing
primal:  Kruskal on the sub-edges by increasing value
         edges that merge components kill H0;  the others (positive) create cycles
dual:    union-find on the faces
for each positive edge e, by decreasing value:
    let A, B be the components of the two faces around e
    if A == B: essential cycle, continue
    young <- the component whose oldest face has the smaller value
    record (value(e), value(oldest face of young));  merge A and B
```

### 4.5 Diagram summaries and distances (`vectorize`)

- **Objective:** make diagrams usable by statistics, and compare them.
- **Input:** one or two diagrams.
- **Output:** a number or a vector.

| Function | Output |
|---|---|
| `summary` | number of bars, total, mean, max, quantiles of the persistence, entropy, largest gap ratio |
| `largest_gap` | threshold at the largest ratio between consecutive persistences |
| `persistence_image` | Gaussian kernels at (birth, persistence), weighted by persistence, on a grid |
| `betti_curve` | number of bars alive at each scale |
| `wasserstein` | exact p-Wasserstein distance (Hungarian algorithm on the diagonal-augmented cost matrix) |
| `bottleneck` | exact bottleneck distance (binary search on the costs, plus a perfect matching test) |

### 4.6 ToMATo: clustering by density peaks (`clustering`)

**Reference:** Chazal, Guibas, Oudot, Skraba, *Persistence-Based Clustering in Riemannian Manifolds*, JACM 2013, Algorithm 1.

- **Objective:** cluster the loads into towns, without choosing the number of clusters in advance.
- **Input:** `G`, bandwidth `h`, mass of a load (`power` = `|D|` or `count` = 1), merging parameter `τ`.
- **Output:** a cluster label per node, the persistence diagram, the peaks.
- **Density:** `ρ(v) = Σ_l mass(l) · exp(−(d(v, l) / h)²)`, truncated at `3h`.
- **Cost:** `O(n log n + m α(n))` after the distances.
- **Check:** the code gives exactly the same diagram as the independent H0 of §4.2, on the three bandwidths tested on the 2000-node instance.

```text
sort the nodes by decreasing ρ
for v in that order:
    N <- neighbours of v already processed
    if N is empty:   v is a peak, it starts a cluster
    else:
        attach v to the cluster of its densest neighbour      # hill climbing
        for u in N:
            if the clusters of u and v differ:
                low <- the one with the lower peak
                if ρ(peak of low) − ρ(v) < τ:
                    merge low into the other                   # prominence below τ
```

**Choice of τ.** Run once with `τ = +∞`, which gives the full diagram, then choose `τ`:

- `auto_tau` (largest gap in the prominences) is the heuristic suggested in the paper.
- The **bottleneck bootstrap** of AuToMATo (Huber et al., arXiv:2408.06958; `decomposition.bootstrap_tau`) tests whether peaks are significant:

```text
D <- diagram of ρ
repeat B times:
    resample the loads with replacement, recompute ρ and its diagram D*
    d_i <- bottleneck distance between the finite points of D* and D
τ <- 2 × the (1 − α) quantile of the d_i           # α = 0.35 by default
significant peaks = peaks of prominence >= τ
```

### 4.7 Split tree: the corridors of a cell (`merge_tree`)

- **Objective:** describe how the branches of the network join on their way to the source, to group the loads of the same corridor in the same feeder.
- **Input:** a cell `H`, its source, `f` = distance to the source inside the cell, the loads.
- **Output:** a tree. A leaf is the end of a corridor, an internal node a junction, and each branch carries its loads.
- **Persistence of a branch** = length of the corridor before it joins another one. `simplify(τ)` removes branches shorter than `τ`.

```text
sort the nodes of H by decreasing f
for v in that order:
    merge the components of the already processed neighbours of v
    T <- tree nodes of those components that contain loads
    if |T| >= 2:   create a junction node at f(v), children T
    if v is a load:   attach it to the current branch (new leaf if there is none)
root the tree at the source
```

### 4.8 Lower bounds and necessary conditions (`bounds`)

**Lower bound on `J`.** With the sources contracted to one point, any feasible network is a Steiner tree on the `t` points (loads plus root).

- **Input:** the instance. **Output:** `lb = max(lb_mst, lb_nn)`, valid for the problem without feeder limits. Cost: a few milliseconds even at 4000 nodes.

```text
MST <- sum of the H0 death times of the loads, sources merged     # §4.3
lb_mst <- MST · t / (2 (t − 1))             # MST <= 2 (1 − 1/t) Steiner*
lb_nn  <- 1/2 · sum of the nearest-neighbour distances            # Euler tour
return max(lb_mst, lb_nn)
```

**Necessary conditions** (reach filtration at radius `L_max`).

- **Input:** the instance. **Output:** the list of reasons why it is infeasible (empty: no proof), the number of connected components of the reach graph.

```text
allowed[l] <- sources s with d(s, l) <= L_max
1. a load with no allowed source: infeasible
2. a load with |D| > P_max for every allowed source: infeasible
3. transportation LP: x[l, s] >= 0 on allowed pairs, sum_s x[l, s] = 1,
       |sum_l D[l] x[l, s]| <= min(capacity, n_departs_max · P_max):  infeasible if the LP is
4. for each source s: forced(s) <- loads whose only allowed source is s
       if the lower bound on forced(s) + s exceeds n_departs_max · L_max: infeasible
```

### 4.9 Spectral tools (`spectral`)

The graph is viewed as an electrical network: a line of length `l` has conductance `1 / l`.

| Function | Objective | Input → Output |
|---|---|---|
| `spectral_embedding` | Coordinates that reveal bottlenecks | `G, k` → first `k` non-trivial eigenvectors of the normalised Laplacian, eigenvalues |
| `spectral_edge_costs` | Make cuts fall at bottlenecks | `G, β` → edge attribute `length · (1 + β · min(jump / median jump, 10))` |
| `anchored_spectral_cells` | Spectral clustering anchored on the sources | `G` → `{source: nodes}`, connected |
| `effective_resistance` | Redundancy of the routes to a load | `G`, sources → `R(s, v) ≤ d(s, v)` |
| `kron_reduction` | Reduce the network onto some nodes | `G`, nodes → Schur complement of the Laplacian |

```text
ANCHORED_SPECTRAL_CELLS(G):
    X <- spectral embedding of dimension |S|, rows normalised to length 1
    each node goes to the source nearest in X
    each cell keeps only the connected component of its source
    the other nodes go to the nearest cell, by Dijkstra on the true lengths
```

### 4.10 Exact graph reductions (`reductions`)

- **Objective:** make the graph smaller without losing an optimal solution.
- **Input:** `G`. **Output:** reduced graph, each edge remembering the original path (`expand_edges` maps a solution back).

```text
repeat until nothing changes, for each transit node v:
    degree 1:  remove v                      # a tree never needs a dead end
    degree 2, neighbours a, b:
        replace a - v - b by one edge of length l(a, v) + l(v, b)
        if edge a - b already exists: keep the shorter one
```

### 4.11 Feeder design (`departs`)

- **Objective:** choose the feeders of a cell explicitly, instead of inheriting them from the shape of a single Steiner tree.
- **Input:** a cell, its source, `P_max`, `L_max`, `n_departs_max`, a method.
- **Output:** the edges of the tree of the source, one subtree per group of loads.
- **Estimate of a group `A`:** `est(A) = MST(A) + gate(A)`, with `gate(A)` = distance from the source to `A`. It bounds from above the length of a Steiner tree on `A` plus the source. When two groups merge, `est` changes by `d(A, B) − max(gate(A), gate(B))`, the Esau-Williams **saving**.

| Method | Order of the merges |
|---|---|
| `kruskal` | increasing `d(A, B)`: the H0 filtration of the loads, truncated by the capacities |
| `savings` | decreasing saving |
| `mergetree` | savings, but first inside each branch of the split tree (§4.7), bottom-up |

```text
GROUP_LOADS(cell):
    start with one group per load
    while a pair (A, B) is mergeable with |P(A ∪ B)| <= P_max and est(A ∪ B) <= L_max:
        merge the best pair for the chosen order
    while there are more groups than n_departs_max:
        merge the pair that increases est the least           # even beyond the limits

BUILD_FEEDERS(cell, groups):
    for each group, largest est first:
        Steiner tree (Mehlhorn) on group + source, in the cell minus the nodes already used
        if that is impossible: build in the whole cell, then repair the union (MST, prune)
```

### 4.12 Variable-distance Voronoi (`partition`)

- **Objective:** assign the nodes to the sources so that each source respects its limits.
- **Input:** `G`, initial weights `{source: w}`, a tree builder, a distance attribute.
- **Output:** the best solution found, with cells, trees, certificate, history.

```text
VARIABLE_DISTANCE(w):
    for k = 0 .. n_iter − 1:
        cells <- Voronoi cells with edge length / w_s from source s
        trees <- tree builder on the cells;  cert <- verify.py
        keep the best solution so far;  stop if feasible
        r_s <- worst value / limit among the constraints of source s
        for each overloaded source (r_s > 1):
            w_s <- w_s · max(r_s^(−step · decay^k), 1/2)            # the cell shrinks
        if grow:  for each source with r_s < 1 next to an overloaded cell:
            w_s <- w_s · min((1 / r_s)^(step · decay^k / 2), 2)     # the neighbour takes over
```

**ToMATo initial weights (`tomato_init_weights`).** The clusters of §4.6 are assigned to sources by a transportation LP that minimises power × distance under the capacities. Weights are then adjusted so that each cell holds the power of its clusters. The intent is that a town is not cut by a cell boundary.

### 4.13 Local search (`local_search`)

- **Objective:** improve a solution, in lexicographic order: fewer violations, then less excess, then shorter.
- **Input:** a solution, candidate edges, candidate loads, a time limit.
- **Output:** the improved solution and statistics.

```text
repeat until no move improves or the time limit:
    for each candidate edge (u, v) of a same tree:             # INSERT
        add it, remove the longest key path of the cycle it closes
        keep the move if the score of the source improves
    for each candidate leaf load l:                            # RELOCATE
        detach l and its branch, reconnect it by a shortest path through free nodes
        to the nearest node of each tree, or to the source (new feeder)
        keep the first move that improves the score of the two sources
```

**Where TDA comes in: filters on the candidates.** Insert moves can be limited to edges that close a persistent cycle (§4.4), and relocate moves to loads whose two nearest sources are almost equidistant (`boundary_margin`, the cut locus of the Voronoi cells).

### 4.14 Instance features and configuration selection (`features`)

- **Objective:** predict which configuration will work best on an instance.
- **Input:** an instance, or a feature table with labels and groups.
- **Output:** a vector of features; for `knn_select`, the accuracy against always predicting the majority.
- **Features:** size and slack (`simple_features`); summaries of the four diagrams, the lower bound, the spectral gap and the effective resistance (`tda_features`), normalised by `L_max` and the top density.

```text
KNN_SELECT(X, y, groups, k = 3):
    for each instance i:
        train on the instances of the other groups, standardised on that part
        predict the majority label among the k nearest neighbours of i
    return the accuracy and the accuracy of the majority baseline
```

### 4.15 ToMATo decomposition (`decomposition`)

- **Objective:** split a large instance into smaller, independently solvable subproblems.
- **Input:** instance, bandwidth `h`, mass, `max_sources`, margin (default 0.85).
- **Output:** a partition of the nodes into connected regions, with statistics per region, and the sub-instances as files.
- **Why it departs from the paper:** without a clear diagram, §3 of the paper says `τ` must come from application-specific information. Here that information is electrical. This merging rule is our own adaptation.

```text
atoms <- ToMATo basins with τ = 0, every node labelled
# 1. sources are markers
for margin in (0.70, 0.85):
    by decreasing saddle: an atom without source joins a neighbouring region with a source,
    if the net power of the union stays <= margin · supply
# 2. group regions, by decreasing saddle, up to max_sources per region, within the margin
# 3. repair: while a region has no source, no load, net power > supply, or an unreachable load:
    merge the region of lowest peak into a neighbour
    (prefer a neighbour that makes the union self-sufficient, then the highest saddle)
solve each region separately;  put the trees together;  verify.py on the whole graph
```

Regions partition the nodes, so the union of the trees is node-disjoint by construction.

### 4.16 Persistence-based cuts for the PLNE (`plne_cuts`)

This is the new link between topology and exact optimisation, tested in §5.

- **Objective:** strengthen the LP relaxation exactly where it is weak.
- **Input:** the instance. **Output:** a set of valid inequalities.

**Why.** Experiment 5 showed that the root LP of the tightened model equals the optimum *without* feeder limits: the LP spreads each load across all the feeders, so 8 MW and 50 km never bind.

**The inequality.** Let `S` be a node set without source, and `C` the loads inside it. Each feeder that serves a load of `C` is a tree rooted outside `S`, so it has its own arc entering `S`:

```text
sum over feeders k, sum over arcs a entering S of  x[a, k]  >=  r(S)
```

for any lower bound `r(S)` on the number of feeders that serve `C`. We use the maximum of:

- power: `⌈ max(P(C) − |productions outside C|, −P(C) − consumptions outside C) / P_max ⌉`. A feeder can mix loads of `C` with loads outside it, so what could compensate is subtracted.
- length: `⌈ moat(C) / L_max ⌉`, where `moat(C) = ½ (MST + longest MST edge)` is a lower bound of the Steiner tree on `C` plus the contracted sources.

**Validity does not depend on the choice of `S`.** Persistence only answers the question "which of the `2^n` sets are worth trying?". Each hierarchy provides about one set per load.

```text
sets <- {}
# H0 family: single-linkage of the loads
for each cluster in the merge order of the loads (§4.3):
    add S = union of the Voronoi regions of its loads
# ToMATo family
for each bandwidth h and mass:
    for each union of basins, by decreasing saddle:
        add S = that union, sources removed
add S = V without the sources
keep the sets with r(S) >= 2 and add one inequality each
```

## 5. Results

Full tables and limits: [Experiment 7](../experiments/experiment_7_tda_heuristics_and_plne.md).

**Heuristics on the 54 Petit instances.** The gap is measured on the Petit_N60 instances where the heuristic is valid.

| Heuristic | Valid (of 54) | Mean time (s) | Mean gap on N60 |
|---|---|---|---|
| Portfolio of all heuristics | 51 | 3.13 | 4.2% |
| Merge tree + local search | 50 | 0.20 | 3.9% |
| Merge tree feeders | 50 | 0.18 | 9.6% |
| Voronoi variable distance | 23 | 0.44 | 10.3% |

**PLNE variants on the 18 Petit_N60 instances** (1 thread, limit of 120 Gurobi work units).

| Variant | Proven optimal | Mean work |
|---|---|---|
| Tightened model alone | 10 | 57.5 |
| + warm start | 11 | 53.6 |
| + persistence cuts | 13 | 44.7 |
| + warm start and cuts | 13 | 41.2 |
| Restricted edges, then warm start and cuts | 14 | 38.8 |

**Decomposition.** On one TrèsGrand instance, the bootstrap finds a single significant density peak, and the regions stay too large (up to 392 loads). See `output/tomato_decomposition/`.
