# Exact methods explained: Gurobi PLNE, dynamic program, duality

This note explains the three methods used to find or bound the optimal network. Each part gives a plain-English explanation, pseudo-code, and real code taken from the project.

Results on the instances are in [Experiment 5](../experiments/experiment_5_gurobi_lower_bounds_petit_n60.md).

## 1. The problem they all solve

We have a graph of candidate lines with lengths `L_e` (km). Nodes are **sources**, **loads** (consumers `D > 0`, producers `D < 0`, and consumers of 0 MW) or **transit** nodes. We must build a set of lines of minimal total length `J` such that:

- every load is connected to exactly one source, in a radial network (a forest, no cycle, no node fed by two sources);
- no tree passes through a second source.

A **feeder** (départ) is the sub-tree hanging from one neighbour of a source. Each one must satisfy:

| Constraint | Scope | Limit on the test instances |
|---|---|---|
| Total length of the feeder | each feeder | ≤ 50 km |
| Absolute net power of the feeder | each feeder | ≤ 8 MW |
| Number of feeders | each source | ≤ 8 |
| Absolute net power of the source | each source | ≤ 70 MW |

These are the constraints checked by `verify.py`. The per-line power limit and the P×L limit in the instance files are not enforced by any of the three methods.

### The three methods at a glance

| | Gurobi PLNE | Dynamic program | Dual bound |
|---|---|---|---|
| What it gives | A network and a proof of optimality (if it finishes) | A lower bound, exact when its feeders do not overlap | A lower bound with a certificate |
| Main idea | Yes/no variables, solved by branch and bound | Solve small subsets of loads first, reuse them | Give each load a price, under the cost of every allowed feeder |
| Cost | Unpredictable, solver-dependent | About 3^(loads) operations, 2^(loads) memory | One LP with up to one column per feeder |
| Works on large instances | Partially, time grows quickly | No (see part 6) | Not as implemented, needs column generation |

## 2. Gurobi PLNE

### Plain English

We describe the network as a list of yes/no questions: "is this line built for this feeder?", "is this load served by this feeder?". Linear constraints say which combinations are valid, and the objective adds up the length of the built lines.

Gurobi solves this with **branch and bound**. It first ignores the yes/no requirement and lets variables take fractional values, such as "half of this line is built". That easier problem (the **LP relaxation**) gives a lower bound. If the answer has a fractional variable, Gurobi splits into two cases, variable = 0 and variable = 1, and repeats. A case can be dropped when its bound is no better than the best network already found.

The solver stops when the lower bound meets the best network found (the gap is zero), which proves optimality.

### Model

Each source has several feeders `k`. For a feeder `k` of source `s(k)`:

| Variable | Type | Meaning |
|---|---|---|
| `y[v,k]` | binary | load `v` is served by feeder `k` |
| `z[e,k]` | binary | line `e` is built for feeder `k` |
| `x[(i,j),k]` | binary | arc `i → j` is used in feeder `k`, oriented away from the source |
| `w[i,k]` | binary | node `i` belongs to the tree of feeder `k` |
| `f[(i,j),k]` | continuous ≥ 0 | flow sending one unit to each node of the tree, to force connectivity |

```text
minimize    J = sum over e, k of  L[e] * z[e,k]

subject to  for each load v:          sum over k of y[v,k] = 1
            for each feeder k:
                -Pmax <= sum over v of D[v] * y[v,k] <= Pmax      (feeder power)
                sum over e of L[e] * z[e,k] <= Lmax               (feeder length)
                at most one arc leaves the source                  (one feeder = one head)
                a tree rooted at s(k):
                    x[(i,j),k] + x[(j,i),k] = z[{i,j},k]          (a line has one direction)
                    arcs entering s(k) = 0, and entering i = w[i,k]   (one parent per node)
                    flow in - flow out at i = w[i,k]              (each node consumes one unit)
                    f[a,k] <= (|V| - 1) * x[a,k]                  (flow only on open arcs)
                    y[v,k] <= w[v,k]                              (a served load is in the tree)
            for each node i:              sum over k of w[i,k] <= 1   (feeders are node-disjoint)
            for each source r:            |sum of D[v]*y[v,k] for k of r| <= capacity[r]
```

### Pseudo-code

```text
BUILD_PLNE(instance):
    create the variables y, z, x, w, f for each feeder
    add the objective and the constraints above
    if tighten:
        keep at most min(n_departs_max, degree of s) feeders per source s
        add the valid inequalities of TIGHTEN (below)
    return the model

BRANCH_AND_BOUND(model, stop_value):
    best_network <- none,  best_value <- +infinity
    open_nodes <- { the original problem }
    while open_nodes is not empty and time remains and best_value > stop_value:
        take a node N from open_nodes
        solve the LP relaxation of N, giving value v and solution x*
        if v >= best_value: discard N                      # cannot beat the best network
        else if x* is integer:
            best_network <- x*, best_value <- v            # a real network
        else:
            add cuts, run heuristics (may update best_network)
            pick a fractional binary variable b
            add N + {b = 0} and N + {b = 1} to open_nodes
    lower_bound <- smallest LP value among open_nodes
```

### Real code: the variables and objective ([plne_gurobi.py:73](../plne_gurobi.py#L73))

```python
    m = gp.Model("network_planning")
    y = add_vars(m, "y", [(v, k) for v in Vc for k in K], GRB.BINARY)
    z = add_vars(m, "z", [(e, k) for e in E for k in K], GRB.BINARY)
    x = add_vars(m, "x", [(a, k) for a in A for k in K], GRB.BINARY)
    w = add_vars(m, "w", [(i, k) for i in V for k in K], GRB.BINARY)
    f = add_vars(m, "f", [(a, k) for a in A for k in K])
    length = {k: gp.quicksum(L[e] * z[e, k] for e in E) for k in K}

    # (2.8) objectif
    m.setObjective(gp.quicksum(length[k] for k in K), GRB.MINIMIZE)
```

### Real code: the tree of each feeder

```python
        # (2.5) structure d'arbre enraciné en s(k)
        m.addConstr(w[r, k] == 1)
        for (i, j) in E:
            m.addConstr(x[(i, j), k] + x[(j, i), k] == z[(i, j), k])
        for i in V:
            inflow = gp.quicksum(x[(j, i), k] for j in G.neighbors(i))
            if i == r:
                m.addConstr(inflow == 0)
            else:
                m.addConstr(inflow == w[i, k])
                m.addConstr(gp.quicksum(f[(j, i), k] for j in G.neighbors(i))
                            - gp.quicksum(f[(i, j), k] for j in G.neighbors(i)) == w[i, k])
        for a in A:
            m.addConstr(f[a, k] <= (nV - 1) * x[a, k])
```

The flow is what prevents disconnected loops: each node of the tree consumes one unit, and the units can only reach it from the source along open arcs.

### Tightening the model

The model above is correct but its LP relaxation is weak. The LP can open an arc at `x = 1/59` and still carry one unit of flow, so it pays very little length for connectivity. The feeder limits also barely bind, because the LP can spread each load across all feeders. On the 18 Petit_N60 instances, the root bound of the tightened model equals the optimum **without any feeder limit**.

The tightening adds constraints that keep every optimal network but remove many fractional solutions:

```text
TIGHTEN:
    for each feeder k and arc a:           f[a,k] >= x[a,k]
    for each feeder k and transit node i:  (children of i in k) >= w[i,k]
    for each load t:
        add g[a] in [0,1], the flow of one unit from the sources to t
        g[a] <= sum over k of x[a,k]            # flow only on built arcs
        conservation of g from the sources to t
```

Why each one is safe:

- `f ≥ x`: an open arc feeds a node, and that node consumes at least one unit.
- A transit node with no child is a useless branch and can be removed to shorten the network, so an optimal network has none.
- A built network contains a path from a source to every load, so the unit of flow `g` always exists.

The last family is the strong one. A fractional opening of arcs must now carry a full unit to each load, which removes the `1/59` effect.

### Real code: the tightening ([plne_gurobi.py:180](../plne_gurobi.py#L180))

```python
    G, x, w, f, K = inst["G"], var["x"], var["w"], var["f"], inst["feeders"]
    sources = set(inst["sources"])
    arcs = [(i, j) for (i, j) in var["E"]] + [(j, i) for (i, j) in var["E"]]
    for k in K:
        for a in arcs:
            m.addConstr(f[a, k] >= x[a, k])
        for i in inst["transship"]:
            m.addConstr(gp.quicksum(x[(i, j), k] for j in G.neighbors(i)) >= w[i, k])
    for t in inst["prosumers"]:
        usable = [a for a in arcs if a[1] not in sources and a[0] != t]
        g = add_vars(m, f"g{t}", usable, ub=1.0)
        for a in usable:
            m.addConstr(g[a] <= gp.quicksum(x[a, k] for k in K))
        m.addConstr(gp.quicksum(g[a] for a in usable if a[0] in sources) == 1)
        for v in G.nodes:
            if v not in sources:
                m.addConstr(gp.quicksum(g[(v, j)] for j in G.neighbors(v) if (v, j) in g)
                            - gp.quicksum(g[(j, v)] for j in G.neighbors(v) if (j, v) in g) == (-1 if v == t else 0))
```

### Early stop at the dynamic-program bound

A valid network whose length equals a lower bound is optimal, because no network can be shorter. So Gurobi can stop as soon as it finds such a network, without closing its own gap.

`BestObjStop` tells Gurobi to stop as soon as the best network is at or below a value ([plne_gurobi.py:293](../plne_gurobi.py#L293)):

```python
    if dp_bound is not None:
        params = dict(params or {}, BestObjStop=dp_bound + DP_TOL / 2)
        res['dp_bound'] = dp_bound
    res.update(optimize(m, time_limit, threads, log_file, console, params))
```

The claim is only made after `verify.py` has accepted the network ([plne_gurobi.py:222](../plne_gurobi.py#L222)):

```python
    lb = res.get('dp_bound')
    if lb is not None and cert['feasible'] and cert['J'] <= lb + DP_TOL:
        res.update(status='OPTIMAL_BY_DP', bound=lb, gap=0.0)
```

## 3. Dynamic program

### Plain English

Two ideas are combined.

**Idea 1, the shortest tree for each group of loads.** Imagine a tree connecting a group of loads `S` and a vertex `v`. From `v`, the tree goes along a path to a vertex `u` where it splits into two parts, each serving part of the group. So the best tree for a large group is built from the best trees of two smaller groups. If those answers are stored, nothing is recomputed. This is the **Dreyfus-Wagner** algorithm.

**Idea 2, the best split of loads into feeders.** From the shortest trees, we know the shortest possible feeder for each group, and whether it respects the power and length limits. We then look for the cheapest way to split all loads into allowed feeders.

Each group is stored as a binary number. With loads `A, B, C`, the number `0b101` means the group `{A, C}`. This turns "the best tree for this group" into an array lookup.

**Example with 3 loads.** The groups are built from smallest to largest: `{A}`, `{B}`, `{C}`, then `{A,B}`, `{A,C}`, `{B,C}`, then `{A,B,C}`. To get `{A,B,C}` at vertex `v`, the program compares `{A} + {B,C}`, `{A,B} + {C}` and `{A,C} + {B}`, all of which are already known.

### Pseudo-code

```text
Input: graph G, loads L (k = |L|), sources R
H <- G without the sources
D <- all-pairs shortest distances in H

STEP 1, shortest trees (Dreyfus-Wagner):
    for each load i and vertex v:   dp[{i}][v] <- D[i][v]
    for each group S, smallest first, with |S| >= 2:
        for each vertex v:
            merged[v] <- min over S1 inside S (S1 contains the first load of S, S1 != S)
                          of  dp[S1][v] + dp[S \ S1][v]
        for each vertex v:
            dp[S][v] <- min over u of  merged[u] + D[u][v]

STEP 2, shortest feeder for each group and source r:
    cost[r][S] <- min over neighbours h of r of  L(r,h) + dp[S][h]
    allowed[r][S] <- |P(S)| <= Pmax  and  cost[r][S] <= Lmax

STEP 3, best way for source r to serve a group M:
    F[r][M] <- cost[r][M] if allowed, else +infinity
    repeat min(n_departs_max, degree of r) - 1 times:
        for each group M:
            F[r][M] <- min(F[r][M],
                           min over allowed S inside M, containing the first load of M,
                               of  cost[r][S] + previous F[r][M \ S])
    F[r][M] <- +infinity if |P(M)| > capacity[r]

STEP 4, assignment of the loads to the sources:
    B[empty] <- 0
    for each source r:
        B[M] <- min over S inside M of  F[r][S] + previous B[M \ S]

Return B[all loads]
```

### Why the result is a lower bound

Take any valid network, and let `F₁ … F_m` be its feeders, serving groups of loads `S₁ … S_m`. Every load is in one group.

1. `cost[S_i]` is the shortest possible feeder for `S_i`, so `length(F_i) ≥ cost[S_i]`.
2. The network is valid, so `S_i` satisfies the power limit, and `cost[S_i] ≤ length(F_i) ≤ L_max`. So `S_i` is allowed.
3. The groups form a split considered by steps 3 and 4, so `Σ cost[S_i] ≥ B[all loads]`.

Therefore `J ≥ B[all loads]` for every valid network.

The program builds each feeder independently, so **two feeders can use the same node**. If its best split has no overlap, its feeders form a valid network and the bound is the true optimum. If they overlap, it is only a lower bound.

### Real code: the shortest trees ([subset_dp.py:90](../subset_dp.py#L90))

```python
    dp = np.full((full + 1, n), np.inf)
    for i, t in enumerate(loads):
        dp[1 << i] = D[index[t]]
    for mask in range(1, full + 1):
        if mask & (mask - 1):
            dp[mask] = (merge(dp, mask)[None, :] + D).min(axis=1)
```

`mask & (mask - 1)` is non-zero when `mask` has at least two loads. `merge` computes all the combinations of two smaller groups ([subset_dp.py:60](../subset_dp.py#L60)):

```python
    subs = low_submasks(mask)
    best = np.full(dp.shape[1], np.inf)
    for start in range(0, len(subs), BATCH):
        part = subs[start:start + BATCH]
        best = np.minimum(best, (dp[part] + dp[mask ^ part]).min(axis=0))
    return best
```

`low_submasks` keeps the lowest load in `S1`, so `{A} + {B}` and `{B} + {A}` are not computed twice. `mask ^ part` is the complementary group.

### Real code: feeders, sources and the final bound ([subset_dp.py:90](../subset_dp.py#L90))

```python
    for r in sources:
        neighbours = [h for h in G[r] if h in index]
        lim = source_limits(graph, r)
        c = (dp[:, [index[h] for h in neighbours]] + np.array([G[r][h]['weight'] for h in neighbours])).min(axis=1)
        c[0] = np.inf
        allowed = (np.abs(power) <= lim['p_max_depart'] + 1e-9) & (c <= lim['l_max_depart'] + 1e-9)
        allowed[0] = False
        F = source_table(c, allowed, int(min(lim['n_departs_max'], len(neighbours))), power, lim['capacity'])
        new = np.full(full + 1, np.inf)
        for mask in range(full + 1):
            subs = all_submasks(mask)
            new[mask] = float((F[subs] + best[mask ^ subs]).min())
        best = new
```

`source_table` adds one feeder at a time, up to the maximum number of feeders ([subset_dp.py:70](../subset_dp.py#L70)):

```python
    F = np.where(allowed, c, np.inf)
    F[0] = 0.0
    for _ in range(n_max - 1):
        previous = F.copy()
        for mask in range(1, len(F)):
            if mask & (mask - 1) == 0:
                continue
            subs = low_submasks(mask)
            subs = subs[allowed[subs]]
            if len(subs):
                F[mask] = min(F[mask], float((c[subs] + previous[mask ^ subs]).min()))
        if np.array_equal(F, previous):
            break
```

## 4. Duality

### Plain English

Duality gives a lower bound **by setting prices**. Give each load `i` a price `π_i`. Require that, for every allowed feeder, the sum of the prices of its loads is at most the cost of that feeder. Then no network can be cheaper than the sum of all prices. This is because a network is a set of feeders, each costing at least the prices of its loads, and each load is in exactly one feeder.

The best prices give the best bound: we maximize the sum of the prices. This is a linear program, and a feasible set of prices is a **certificate**. Anyone can recheck it by testing the inequalities, without trusting the solver.

Duality is also what Gurobi uses: the LP value at a node of its search tree is the optimal value of the dual of that LP. So Gurobi's lower bound and our dual bound are the same idea, applied to different formulations.

### The feeder LP and its dual

Primal: one variable `λ[r,S]` for each source `r` and allowed feeder group `S`, equal to 1 if the feeder is used.

```text
minimize    sum over (r,S) of  cost[r][S] * λ[r,S]

subject to  for each load i:   sum over (r,S) with i in S of  λ[r,S] = 1        (price π_i, free)
            for each source r: sum over S of λ[r,S] <= nmax_r                   (price μn_r <= 0)
                               sum over S of P(S) * λ[r,S] <= cap_r             (price μp_r <= 0)
                              -sum over S of P(S) * λ[r,S] <= cap_r             (price μm_r <= 0)
            λ >= 0
```

Its dual, the **lower bound** we compute:

```text
maximize    sum over i of π_i  +  sum over r of ( nmax_r * μn_r + cap_r * μp_r + cap_r * μm_r )

subject to  for each source r and allowed S:
                sum over i in S of π_i  +  μn_r  +  P(S) * (μp_r - μm_r)  <=  cost[r][S]
            μn, μp, μm <= 0
```

### Why it is a lower bound

Any valid network gives values of `λ` satisfying the primal constraints. By LP weak duality, its cost is at least the value of any dual-feasible solution. So the dual value is a valid lower bound.

The LP does not require `λ` to be integer, so a fractional split of loads can cost less than any real one. This is the **duality gap**. On the 18 Petit_N60 instances, the gap was zero on 14.

### Pseudo-code

```text
build one column for each (source r, allowed group S), of cost cost[r][S]
solve the LP, read the prices π, μ from the dual solution
lower_bound <- sum of π_i + sum over r of ( nmax_r*μn_r + cap_r*μp_r + cap_r*μm_r )

CERTIFICATE (no trust in the solver):
    for each column (r,S):
        reduced_cost <- cost[r][S] - (prices of the loads of S) - μn_r - P(S)*(μp_r - μm_r)
        check reduced_cost >= -epsilon
```

### Real code

This part was run with a temporary script (`n60_bounds.py`, in the session scratchpad), and is **not in the repository**. The prices and the certificate:

```python
        lp = linprog(c_vec, A_eq=A_eq, b_eq=np.ones(k), A_ub=A_ub, b_ub=np.array(b_ub, float), bounds=(0, None), method='highs')
        res['lp_status'] = lp.message[:80]
        if lp.status == 0:
            pi, mu = lp.eqlin.marginals, lp.ineqlin.marginals
            rc = c_vec - A_eq.T.dot(pi) - A_ub.T.dot(mu)
            res.update(lp_bound=float(lp.fun), lp_fractional=int(((lp.x > 1e-9) & (lp.x < 1 - 1e-9)).sum()),
                       n_columns=j, dual_bound=float(pi.sum() + np.dot(mu, b_ub)),
                       dual_max_violation=float(max(0.0, -rc.min())), dual_sign_ok=bool((mu <= 1e-9).all()))
```

`rc` holds the reduced costs, and `dual_max_violation` is the largest amount by which a feeder violates its inequality. The largest value observed on the 18 instances was about 10⁻¹⁴, which is rounding error.

### Scaling it up: column generation (not implemented)

Listing every allowed feeder is only possible with few loads, because there are up to 2^(loads) of them per source. **Column generation** solves the same LP while creating feeders only when needed:

```text
start with a small set of feeders
repeat:
    solve the LP on the current feeders, read the prices π, μ
    look for a feeder whose prices exceed its cost (reduced cost < 0)
    if none exists: the LP is solved, stop
    add that feeder to the set
```

The search for a new feeder is itself a hard problem, a Steiner tree with power and length limits. During the iterations the LP value is not yet a valid lower bound, and becomes one only at the end, or after a Lagrangian correction.

## 5. How the three relate

| | Gurobi PLNE | Dynamic program | Dual |
|---|---|---|---|
| Considers | Everything in `verify.py` | Everything except node-disjointness between feeders | Same relaxation as the dynamic program, and integrality is dropped too |
| Works on | Lines, edge by edge | Whole groups of loads | Whole feeders |
| Gives a network | Yes | Sometimes (if no overlap) | No |
| Its value | Equal to the optimum if it finishes | ≤ optimum | ≤ dynamic program ≤ optimum |

The dual is always at most the dynamic program: the dual relaxes one more thing (integrality), and the dynamic program solves that integer problem exactly.

The reason these bounds beat Gurobi's on these instances is the same as in Experiment 5: describing feeders line by line lets the LP blend feeders together, which hides the power and length limits. Working with whole feeders removes that freedom.

## 6. Why the dynamic program stops at 20 loads

`DP_MAX_LOADS = 20` in [plne_gurobi.py](../plne_gurobi.py#L63) is a **practical limit chosen for this project**, not a mathematical one. The program would run for 21 or 22 loads, but its cost becomes unreasonable.

**Time grows by 3 for each extra load.** Each merge chooses, for every load, one of three states: outside the group `S`, in the first part `S1`, or in the second part `S ∖ S1`. That is about `3^k` choices. Each one is an addition over all vertices, and the feeder and source steps add similar work.

**Memory doubles with each extra load.** The table `dp` has one row per subset, so `2^k` rows of `|V|` numbers at 8 bytes each.

Time estimates start from a measured 584 s for 19 loads (`Petit_N75_S3_0.25_Hiver`) and multiply by 3 per extra load, so they vary with the number of sources and feeders.

| Loads | Pairs of groups (≈ 3^k / 2) | Memory of `dp`, 90 vertices | Time |
|---|---|---|---|
| 16 | 2.2 × 10⁷ | 47 MB | 12 to 28 s (measured) |
| 18 | 1.9 × 10⁸ | 189 MB | **162 to 246 s (measured)** |
| 19 | 5.8 × 10⁸ | 377 MB | **584 s (measured)** |
| 20 | 1.7 × 10⁹ | 755 MB | ≈ 29 min |
| 22 | 1.6 × 10¹⁰ | 3.0 GB | ≈ 4.4 h |
| 25 | 4.2 × 10¹¹ | 24 GB | ≈ 5 days |
| 30 | 1.0 × 10¹⁴ | 773 GB | ≈ 3 years |

**Why 20 specifically.** At 20 loads, the program already takes longer than the 20-minute limit given to Gurobi, so a larger value would spend more time on the bound than on the solve. At 22 loads it needs about 4.4 hours and 3 GB per run, which is too much when 3 runs share the machine.

This is why Moyen instances, with 22 to 62 loads, have no bound from this program. Going further needs a method that does not enumerate all subsets, such as the column generation in part 4.
