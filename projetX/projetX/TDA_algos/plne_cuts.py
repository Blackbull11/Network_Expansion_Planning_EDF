"""
Valid inequalities for the PLNE, on load sets proposed by the persistence hierarchies.

Exact solvers are slow here because the LP blends the feeders: each load is
spread over all the feeders, so the 8 MW and 50 km limits do not bind at the
root (experiment 5). Rounded capacity cuts bring them back. For a node set S
without source and C the loads inside S, every feeder that serves a load of C
is a tree rooted outside S, so it has its own arc entering S:

    sum over the feeders k and the arcs a entering S of x[a, k]  >=  r(S)

for any lower bound r(S) on the number of feeders that serve C:
  power   : r >= ceil(max(P(C) - |productions outside C|, -P(C) - consumptions outside C) / p_max)
            (a feeder can mix loads of C with loads outside C: the production
            outside C could compensate, so it is subtracted)
  length  : r >= ceil(Steiner(C + sources contracted) / l_max), with the moat
            lower bound 1/2 (MST + longest MST edge) of Steiner (experiment 4)
A cut is only useful if r >= 2. They are valid for any S: persistence only
chooses which S to try among the 2^n possible ones.

Families of sets (each one a hierarchy of nested clusters, so about one set per load):
  h0     : clusters of the single-linkage hierarchy of the loads (H0 of the
           shortest-path metric, terminal_h0), S = union of their Voronoi regions
  tomato : unions of ToMATo basins, in decreasing saddle order, for several bandwidths
"""
import math

from .clustering import power_density, tomato
from .decomposition import Regions, fill_unlabelled
from .persistence import terminal_h0
from .union_find import UnionFind

LOAD_LABELS = ('Conso', 'Prod')
EPS = 1e-9


def required_feeders(graph, loads: list, limits: dict) -> int:
    """r(S) for the loads `loads` of S: lower bound on the number of feeders that serve them (see module docstring)."""
    G = graph.G
    p = sum(G.nodes[l]['conso_ete'] for l in loads)
    inside = set(loads)
    outside = [G.nodes[v]['conso_ete'] for v, lab in G.nodes('label') if lab in LOAD_LABELS and v not in inside]
    prod_out = -sum(d for d in outside if d < 0)
    cons_out = sum(d for d in outside if d > 0)
    r_power = math.ceil(max(p - prod_out, -p - cons_out) / limits['p_max'] - EPS)
    sources = [n for n, lab in G.nodes('label') if lab == 'Source']
    h0 = terminal_h0(G, loads, [sources])
    moat = 0.5 * (h0['total'] + float(h0['deaths'].max())) if len(h0['deaths']) else 0.0
    return max(r_power, math.ceil(moat / limits['l_max'] - EPS), 1)


def _limits(graph) -> dict:
    from verify import source_limits
    lims = [source_limits(graph, s) for s, lab in graph.G.nodes('label') if lab == 'Source']
    return {'p_max': max(l['p_max_depart'] for l in lims), 'l_max': max(l['l_max_depart'] for l in lims)}


def _record(graph, out: dict, nodes: set, family: str, limits: dict, sources: set) -> None:
    S = frozenset(nodes - sources)
    loads = [v for v in S if graph.G.nodes[v]['label'] in LOAD_LABELS]
    if S and S not in out and loads:
        r = required_feeders(graph, loads, limits)
        if r >= 2:
            out[S] = {'nodes': S, 'r': r, 'family': family, 'n_loads': len(loads)}


def h0_sets(graph, out: dict = None) -> dict:
    """{frozenset S: cut} for the clusters of the single-linkage hierarchy of the loads, and the whole graph."""
    G = graph.G
    out = {} if out is None else out
    limits = _limits(graph)
    sources = {v for v, lab in G.nodes('label') if lab == 'Source'}
    loads = [v for v, lab in G.nodes('label') if lab in LOAD_LABELS]
    if len(loads) < 2:
        return out
    h0 = terminal_h0(G, loads)
    members = {}                      # {label: nodes of its Voronoi region}
    for v, lab in h0['region'].items():
        members.setdefault(lab, set()).add(v)
    uf = UnionFind(members)
    for a, b, _ in sorted(h0['mst_edges'], key=lambda e: e[2]):
        ra, rb = uf.find(a), uf.find(b)
        root = uf.union(ra, rb)
        other = rb if root == ra else ra
        members[root] = members[root] | members.pop(other)
        _record(graph, out, set(members[root]), 'h0', limits, sources)
    _record(graph, out, set(G), 'h0', limits, sources)
    return out


def tomato_sets(graph, bandwidths=(3.0, 5.0, 7.0, 10.0), masses=('power', 'count'), out: dict = None) -> dict:
    """{frozenset S: cut} for the unions of ToMATo basins, in decreasing saddle order."""
    G = graph.G
    out = {} if out is None else out
    limits = _limits(graph)
    sources = {v for v, lab in G.nodes('label') if lab == 'Source'}
    for h in bandwidths:
        for mass in masses:
            rho, _ = power_density(G, h, mass=mass)
            regs = Regions(graph, fill_unlabelled(G, tomato(G, rho, 0.0, min_density=0.0)['label']), rho)
            pairs = sorted(((s, a, b) for a, nbrs in regs.saddle.items() for b, s in nbrs.items() if str(a) < str(b)),
                           key=lambda p: (-p[0], str(p[1]), str(p[2])))
            members = {r: set(ns) for r, ns in regs.members.items()}
            for ns in list(members.values()):
                _record(graph, out, set(ns), 'tomato', limits, sources)
            uf = UnionFind(members)
            for _, a, b in pairs:
                ra, rb = uf.find(a), uf.find(b)
                if ra == rb:
                    continue
                root = uf.union(ra, rb)
                other = rb if root == ra else ra
                members[root] = members[root] | members.pop(other)
                _record(graph, out, set(members[root]), 'tomato', limits, sources)
    return out


def add_cuts(m, var: dict, inst: dict, cuts: dict) -> int:
    """Add the cuts {S: {'r': ...}} to the Gurobi model; returns how many were added."""
    import gurobipy as gp
    x, K = var["x"], inst["feeders"]
    arcs = [(i, j) for (i, j) in var["E"]] + [(j, i) for (i, j) in var["E"]]
    for n, (S, cut) in enumerate(cuts.items()):
        entering = [a for a in arcs if a[1] in S and a[0] not in S]
        m.addConstr(gp.quicksum(x[a, k] for a in entering for k in K) >= cut['r'], name=f"tda_cut_{n}")
    return len(cuts)
