"""
Lower bound on J and necessary conditions of feasibility, from the H0
persistence of the loads.

Lower bound (lower_bound). Contract all the sources into one root: a feasible
forest becomes a Steiner tree of the loads + root, so J* >= Steiner*. Let t be
the number of points (loads + root) and MST the weight of the minimum spanning
tree of their metric closure (= sum of the H0 death times, terminal_h0):
  - MST <= 2 (1 - 1/t) Steiner*      =>  J* >= MST t / (2 (t - 1))
  - an Euler tour of the optimal tree visits every point and has length
    2 Steiner*, each step being >= the distance to the nearest other point
                                     =>  J* >= 1/2 sum of the nearest-neighbour distances
Both are valid; the bound is their max. They ignore the capacities and the
limits on the départs: they bound the uncapacitated problem.

Necessary conditions (necessary_conditions). Reach filtration: the union of
the balls of radius l_max_depart around the sources (H0 at r = l_max). A load
can only be served by the sources whose ball contains it, so
  1. a load reachable by no source is infeasible (already in verify.py);
  2. a load with |demande| > p_max_depart of every source that reaches it is infeasible;
  3. transportation relaxation: fractional assignment x[l, s] >= 0 of the
     loads to the sources that reach them, sum_s x[l, s] = 1, and
     |sum_l demande(l) x[l, s]| <= min(capacity, n_departs_max p_max_depart).
     If this LP is infeasible, so is the instance. It subsumes the global
     power balance of verify.necessary_conditions;
  4. length budget: the loads reachable only by s (forced loads) are covered
     by at most n_departs_max départs of length <= l_max, so the lower bound
     above on loads + s must be <= n_departs_max l_max.
"""
import math

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix

from .persistence import shortest_distances, terminal_h0

LOAD_LABELS = ('Conso', 'Prod')
TOL = 1e-6


def _bound_from_h0(res: dict) -> dict:
    t = res['n_labels']
    lb_mst = res['total'] * t / (2 * (t - 1)) if t > 1 else 0.0
    nn = [d for d in res['nn'].values() if math.isfinite(d)]
    lb_nn = 0.5 * sum(nn) if t > 1 else 0.0
    return {'lb': max(lb_mst, lb_nn), 'lb_mst': lb_mst, 'lb_nn': lb_nn, 'mst': res['total'],
            'n_points': t}


def lower_bound(graph, weight: str = "weight") -> dict:
    """Lower bound on J (see module docstring), with the H0 result under key 'h0'."""
    G = graph.G
    sources = [n for n, lab in G.nodes('label') if lab == 'Source']
    loads = [n for n, lab in G.nodes('label') if lab in LOAD_LABELS]
    res = terminal_h0(G, loads, [sources], weight)
    return {**_bound_from_h0(res), 'h0': res}


def reach(graph, weight: str = "weight") -> tuple:
    """(sources, loads, D, allowed): D[s][l] distance, allowed[l] = sources whose l_max ball contains l."""
    from verify import source_limits
    G = graph.G
    sources = [n for n, lab in G.nodes('label') if lab == 'Source']
    loads = [n for n, lab in G.nodes('label') if lab in LOAD_LABELS]
    _, index, D = shortest_distances(G, sources, weight)
    lims = {s: source_limits(graph, s) for s in sources}
    dist = {s: {l: float(D[i, index[l]]) for l in loads} for i, s in enumerate(sources)}
    allowed = {l: [s for s in sources if dist[s][l] <= lims[s]['l_max_depart'] + TOL] for l in loads}
    return sources, loads, dist, allowed, lims


def transportation_lp(graph, sources, loads, allowed, lims) -> int:
    """Status of the transportation relaxation (0 feasible, 2 infeasible, see scipy linprog)."""
    G = graph.G
    pairs = [(l, s) for l in loads for s in allowed[l]]
    if not pairs:
        return 0 if not loads else 2
    col = {p: k for k, p in enumerate(pairs)}
    A_eq = lil_matrix((len(loads), len(pairs)))
    for i, l in enumerate(loads):
        for s in allowed[l]:
            A_eq[i, col[l, s]] = 1.0
    A_ub = lil_matrix((2 * len(sources), len(pairs)))
    b_ub = np.zeros(2 * len(sources))
    for j, s in enumerate(sources):
        cap = min(lims[s]['capacity'], lims[s]['n_departs_max'] * lims[s]['p_max_depart'])
        b_ub[2 * j] = b_ub[2 * j + 1] = cap if math.isfinite(cap) else 1e12
    for (l, s), k in col.items():
        j = sources.index(s)
        d = G.nodes[l]['conso_ete']
        A_ub[2 * j, k] = d
        A_ub[2 * j + 1, k] = -d
    res = linprog(np.zeros(len(pairs)), A_ub=A_ub.tocsr(), b_ub=b_ub, A_eq=A_eq.tocsr(),
                  b_eq=np.ones(len(loads)), bounds=(0, 1), method='highs')
    return res.status


def necessary_conditions(graph, weight: str = "weight") -> dict:
    """
    Proofs of infeasibility (see module docstring).
    Returns {reasons: list of messages (empty: no proof), lp_status, n_forced: {s: n},
    budget: {s: (lower bound of the forced loads, n_departs_max l_max)},
    n_reach_components: connected components of the reach graph (loads - sources)}.
    """
    from .union_find import UnionFind
    G = graph.G
    sources, loads, dist, allowed, lims = reach(graph, weight)
    reasons = []
    far = [l for l in loads if not allowed[l]]
    if far:
        reasons.append(f"charges hors de portée l_max de toutes les sources : {far}")
    heavy = [l for l in loads if allowed[l] and
             all(abs(G.nodes[l]['conso_ete']) > lims[s]['p_max_depart'] + TOL for s in allowed[l])]
    if heavy:
        reasons.append(f"charges plus fortes que p_max_depart de toutes les sources qui les atteignent : {heavy}")
    status = 2 if far else transportation_lp(graph, sources, loads, allowed, lims)
    if not far and status == 2:
        reasons.append("relaxation de transport infaisable : aucune affectation fractionnaire des charges "
                       "aux sources qui les atteignent ne respecte min(capacité, n_departs_max x p_max_depart)")
    budget, n_forced = {}, {}
    for s in sources:
        forced = [l for l in loads if allowed[l] == [s]]
        n_forced[s] = len(forced)
        if not forced:
            continue
        lb = _bound_from_h0(terminal_h0(G, forced, [[s]], weight))['lb']
        cap = lims[s]['n_departs_max'] * lims[s]['l_max_depart']
        budget[s] = (lb, cap)
        if lb > cap + TOL:
            reasons.append(f"source {s} : ses {len(forced)} charges forcées demandent au moins {lb:.1f} km "
                           f"de lignes > n_departs_max x l_max = {cap:.1f} km")
    uf = UnionFind(list(sources) + list(loads))
    for l in loads:
        for s in allowed[l]:
            uf.union(l, s)
    return {'reasons': reasons, 'lp_status': int(status), 'n_forced': n_forced, 'budget': budget,
            'n_reach_components': len(uf.groups())}
