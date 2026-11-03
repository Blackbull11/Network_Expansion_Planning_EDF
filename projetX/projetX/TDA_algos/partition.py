"""
Variable-distance Voronoi cells (voronoi_variable_distance.py) generalised:
  - any initial weights {source: w} (e.g. tomato_init_weights);
  - any edge attribute for the distance (e.g. spectral.spectral_edge_costs);
  - any tree builder (graph, cells) -> {source: edges}: Steiner tree per cell
    (steiner_builder) or designed départs (departs_builder);
  - grow=True: a source with r_s < 1 next to an overloaded cell also has its
    weight multiplied by (1 / r_s) ** (step decay^k / 2) (at most x2), so that
    it takes over the boundary nodes the overloaded cell gives up.

tomato_init_weights: weights such that the Voronoi cells hold the ToMATo
clusters of the loads (clustering.py) assigned to each source by a
transportation LP (clusters are not cut, capacities are respected).
"""
import math

import numpy as np
from scipy.optimize import linprog

from voronoi import WEIGHTINGS, source_weight, voronoi_cells_SB
from steiner import steiner_trees_by_cells
from verify import check_solution, solution_key, source_limits
from voronoi_variable_distance import load_ratio, MAX_SHRINK

from .clustering import cluster_loads
from .departs import departs_by_cells
from .persistence import shortest_distances

MAX_GROW = 2.0


def sources_of(G) -> list:
    return [n for n, lab in G.nodes('label') if lab == 'Source']


def steiner_builder(graph, cells) -> dict:
    return steiner_trees_by_cells(graph.G, cells, log=False)


def departs_builder(method: str = 'mergetree', tau: float = 0.0):
    def build(graph, cells):
        return departs_by_cells(graph, cells, method, tau)
    build.__name__ = f"departs_{method}"
    return build


def standard_inits(graph) -> dict:
    """{name: {source: w}} for the 3 fixed weightings of voronoi.py."""
    return {name: {s: source_weight(graph.G, s, name) for s in sources_of(graph.G)} for name in WEIGHTINGS}


def cell_adjacency(G, cells: dict) -> dict:
    """{source: set of the sources whose cell touches its cell}."""
    owner = {v: s for s, nodes in cells.items() if s != 'unreachable' for v in nodes}
    adj = {s: set() for s in cells if s != 'unreachable'}
    for u, v in G.edges:
        a, b = owner.get(u), owner.get(v)
        if a is not None and b is not None and a != b:
            adj[a].add(b)
            adj[b].add(a)
    return adj


def variable_distance(graph, init_weights: dict, tree_builder=steiner_builder, weight: str = "weight",
                      n_iter: int = 60, step: float = 0.5, decay: float = 0.95, grow: bool = False) -> dict:
    """
    Variable-distance loop (see module docstring). Returns the best solution
    met: {weights, cells, trees, cert, iteration, n_iter, history: list of
    (feasible, number of violations, J) per iteration}.
    """
    G = graph.G
    sources = sources_of(G)
    w = dict(init_weights)
    best, history = None, []
    for k in range(n_iter):
        _, cells, _ = voronoi_cells_SB(G, sources, weight=weight, weighting=w)
        trees = tree_builder(graph, cells)
        cert = check_solution(graph, trees)
        history.append((cert['feasible'], sum(v[0] for v in cert['violations'].values()), cert['J']))
        if best is None or solution_key(cert) < solution_key(best['cert']):
            best = {'weights': dict(w), 'cells': cells, 'trees': trees, 'cert': cert, 'iteration': k}
        if cert['feasible']:
            break
        ratios = {s: load_ratio(src) for s, src in cert['sources'].items()}
        over = {s: r for s, r in ratios.items() if r > 1}
        if not over:
            break
        factor = step * decay ** k
        for s, r in over.items():
            w[s] *= max(r ** -factor, MAX_SHRINK)
        if grow:
            adj = cell_adjacency(G, cells)
            for s, r in ratios.items():
                if r < 1 and adj.get(s, set()) & over.keys():
                    w[s] *= min((1 / max(r, 1e-3)) ** (factor / 2), MAX_GROW)
    best['n_iter'] = k + 1
    best['history'] = history
    return best


def best_of(graph, inits: dict, **kwargs) -> dict:
    """variable_distance from each init of {name: weights}; the best run, with 'init' and 'runs'."""
    runs = {}
    for name, w in inits.items():
        runs[name] = variable_distance(graph, w, **kwargs)
        runs[name]['init'] = name
    best = min(runs.values(), key=lambda r: solution_key(r['cert']))
    return {**best, 'runs': runs}


def assign_clusters(graph, clusters: dict, weight: str = "weight") -> dict:
    """
    {cluster peak: source}: transportation LP minimising sum of gross power x
    distance (peak, source), with |net power of the source| <= min(capacity,
    n_departs_max p_max_depart), rounded to the largest fraction. Falls back
    to the nearest source if the LP is infeasible.
    """
    G = graph.G
    sources = sources_of(G)
    peaks = [p for p in clusters if p != -1]
    if not peaks:
        return {}
    _, index, D = shortest_distances(G, sources, weight)
    dist = np.array([[D[j, index[p]] for j in range(len(sources))] for p in peaks])
    dist[~np.isfinite(dist)] = 1e9
    gross = np.array([max(clusters[p]['gross'], 1e-6) for p in peaks])
    net = np.array([clusters[p]['power'] for p in peaks])
    nc, ns = len(peaks), len(sources)
    c = (gross[:, None] * dist).ravel()
    A_eq = np.zeros((nc, nc * ns))
    for i in range(nc):
        A_eq[i, i * ns:(i + 1) * ns] = 1
    A_ub = np.zeros((2 * ns, nc * ns))
    b_ub = np.zeros(2 * ns)
    for j, s in enumerate(sources):
        lim = source_limits(graph, s)
        cap = min(lim['capacity'], lim['n_departs_max'] * lim['p_max_depart'])
        A_ub[2 * j, j::ns] = net
        A_ub[2 * j + 1, j::ns] = -net
        b_ub[2 * j] = b_ub[2 * j + 1] = cap if math.isfinite(cap) else 1e12
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=np.ones(nc), bounds=(0, 1), method='highs')
    x = res.x.reshape(nc, ns) if res.status == 0 else -dist
    return {p: sources[int(np.argmax(x[i]))] for i, p in enumerate(peaks)}


def tomato_init_weights(graph, n_fit: int = 15, weight: str = "weight", clustering: dict = None) -> dict:
    """
    Weights {source: w} whose Voronoi cells hold the gross power of the
    clusters assigned to each source: w_s *= ((target_s + d) / (mass_s + d)) ** 0.5
    (clipped to [1/2, 2]) for n_fit iterations, mass_s = gross power of the
    loads of the cell. Returns {'weights', 'assignment', 'clustering', 'target'}.
    """
    G = graph.G
    sources = sources_of(G)
    if clustering is None:
        clustering = cluster_loads(G, weight=weight)
    clusters = clustering['clusters']
    assignment = assign_clusters(graph, clusters, weight)
    target = {s: 0.0 for s in sources}
    for p, s in assignment.items():
        target[s] += clusters[p]['gross']
    if -1 in clusters:   # loads out of every cluster: nearest source
        _, index, D = shortest_distances(G, sources, weight)
        for l in clusters[-1]['loads']:
            target[sources[int(np.argmin(D[:, index[l]]))]] += max(G.nodes[l]['conso_ete'], 0.0)
    loads = [n for n, lab in G.nodes('label') if lab in ('Conso', 'Prod')]
    delta = 0.05 * (sum(target.values()) / max(len(sources), 1)) + 1e-9
    w = {s: 1.0 for s in sources}
    for _ in range(n_fit):
        _, cells, _ = voronoi_cells_SB(G, sources, weight=weight, weighting=w)
        owner = {v: s for s, nodes in cells.items() for v in nodes}
        mass = {s: 0.0 for s in sources}
        for l in loads:
            if owner.get(l) in mass:
                mass[owner[l]] += max(G.nodes[l]['conso_ete'], 0.0)
        for s in sources:
            w[s] *= min(max(((target[s] + delta) / (mass[s] + delta)) ** 0.5, 0.5), 2.0)
    return {'weights': w, 'assignment': assignment, 'clustering': clustering, 'target': target}
