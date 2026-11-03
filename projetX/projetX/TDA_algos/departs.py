"""
Design of the départs of each cell, instead of letting them follow the shape
of a single Steiner tree (steiner.py).

1. Group the loads of the cell: each group will be one départ. A group A has
     P(A)   = sum of the demands (|P(A)| <= p_max_depart required),
     est(A) = MST(A) + gate(A), MST of A in the metric closure of the cell
              and gate(A) = distance from the source to A.
   est(A) is an upper bound of the length of the Mehlhorn Steiner tree of
   A + source in the cell (that tree is no longer than the MST of the metric
   closure, itself <= est(A)): a group with est(A) <= l_max_depart gives a
   départ that respects the length limit, as long as its tree does not have
   to avoid the nodes of the départs built before it.
   Merging A and B: MST(A u B) <= MST(A) + MST(B) + d(A, B), gate = min of
   the gates, so the merged estimate saves max(gate) - d(A, B) (Esau-Williams
   saving). d(A, B) = single-linkage distance, updated by min on merges.

   Grouping orders:
     'kruskal'   : merges by increasing d(A, B), i.e. the H0 filtration of the
                   loads truncated by the capacities (capacitated single linkage);
     'savings'   : merges by decreasing saving (Esau-Williams);
     'mergetree' : savings, but merges are first done inside each branch of
                   the split tree of the distance to the source (merge_tree.py),
                   bottom-up: the loads of one corridor are grouped together
                   before corridors are mixed.
   If there are more groups than n_departs_max, the pairs that increase est
   the least are merged, even beyond the limits (verify.py judges).
2. Build one Steiner tree per group (largest first) in the cell minus the
   nodes of the départs already built, so that départs are disjoint. If a
   group cannot reach the source that way, its tree is built in the whole cell
   and the union is turned back into a tree (MST of the union, pruned).
"""
import math

import networkx as nx
import numpy as np
from networkx.algorithms.approximation import steiner_tree

from .merge_tree import split_tree
from .persistence import shortest_distances

LOAD_LABELS = ('Conso', 'Prod')
METHODS = ('kruskal', 'savings', 'mergetree')


class Cell:
    """Distances between the source and the loads of a cell, and group bookkeeping."""

    def __init__(self, G, nodes, source, weight="weight"):
        self.G, self.source, self.weight = G, source, weight
        self.nodes = set(nodes)
        self.H = G.subgraph(self.nodes)
        self.loads = [v for v in self.nodes if G.nodes[v]['label'] in LOAD_LABELS]
        pts = [source] + self.loads
        nodes_, index, D = shortest_distances(self.H, pts, weight, nodes=list(self.nodes))
        self.f = {v: float(D[0, index[v]]) for v in nodes_}         # distance to the source
        self.Dt = D[1:, [index[l] for l in self.loads]]             # load x load distances
        self.gate0 = D[0, [index[l] for l in self.loads]]           # source -> load
        self.pos = {l: i for i, l in enumerate(self.loads)}
        self.demand = np.array([G.nodes[l]['conso_ete'] for l in self.loads], float)


class Group:
    __slots__ = ('terms', 'P', 'mst', 'gate', 'row')

    def __init__(self, cell: Cell, i: int):
        self.terms = [i]
        self.P = cell.demand[i]
        self.mst = 0.0
        self.gate = cell.gate0[i]
        self.row = cell.Dt[i].copy()

    @property
    def est(self):
        return self.mst + self.gate


def _dist(a: Group, b: Group) -> float:
    return float(a.row[b.terms].min())


def _merge(a: Group, b: Group) -> Group:
    d = _dist(a, b)
    a.terms = a.terms + b.terms
    a.P += b.P
    a.mst += b.mst + d
    a.gate = min(a.gate, b.gate)
    a.row = np.minimum(a.row, b.row)
    return a


def _fits(a: Group, b: Group, pmax: float, lmax: float) -> bool:
    d = _dist(a, b)
    return abs(a.P + b.P) <= pmax + 1e-9 and a.mst + b.mst + d + min(a.gate, b.gate) <= lmax + 1e-9


def _key(a: Group, b: Group, order: str):
    """Priority of merging a and b (smaller first), None if the merge is not wanted."""
    d = _dist(a, b)
    if order == 'kruskal':
        return d
    saving = max(a.gate, b.gate) - d
    return -saving if saving > 0 else None


def merge_greedy(groups: list, pmax: float, lmax: float, order: str = 'savings') -> list:
    """
    Merge feasible pairs, best first ('savings': largest positive saving,
    'kruskal': smallest distance), with a lazy heap: entries of a group that
    has changed since are skipped, and the pairs of a merged group are pushed again.
    """
    import heapq
    alive = dict(enumerate(groups))
    version = {i: 0 for i in alive}
    heap = []

    def push(i, j):
        a, b = alive[i], alive[j]
        if _fits(a, b, pmax, lmax):
            k = _key(a, b, order)
            if k is not None:
                heapq.heappush(heap, (k, i, j, version[i], version[j]))

    ids = list(alive)
    for x, i in enumerate(ids):
        for j in ids[x + 1:]:
            push(i, j)
    while heap:
        _, i, j, vi, vj = heapq.heappop(heap)
        if i not in alive or j not in alive or version[i] != vi or version[j] != vj:
            continue
        alive[i] = _merge(alive[i], alive.pop(j))
        version[i] += 1
        for x in alive:
            if x != i:
                push(min(i, x), max(i, x))
    return list(alive.values())


def enforce_nmax(groups: list, nmax: float) -> list:
    """Merge the pairs that increase est the least until there are at most nmax groups."""
    groups = list(groups)
    while len(groups) > max(nmax, 1):
        best, best_inc = None, math.inf
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                a, b = groups[i], groups[j]
                inc = _dist(a, b) + min(a.gate, b.gate) - a.gate - b.gate
                if inc < best_inc:
                    best, best_inc = (i, j), inc
        i, j = best
        groups[i] = _merge(groups[i], groups[j])
        del groups[j]
    return groups


def group_loads(cell: Cell, pmax: float, lmax: float, nmax: float, method: str = 'mergetree',
                tau: float = 0.0) -> list:
    """Groups of loads (lists of nodes) of the cell, one per départ."""
    if not cell.loads:
        return []
    singles = [Group(cell, i) for i in range(len(cell.loads))]
    if method == 'kruskal':
        groups = merge_greedy(singles, pmax, lmax, order='kruskal')
    elif method == 'savings':
        groups = merge_greedy(singles, pmax, lmax, order='savings')
    elif method == 'mergetree':
        mt = split_tree(cell.H, cell.source, cell.f, cell.loads).simplify(tau)
        by_node = {l: singles[cell.pos[l]] for l in cell.loads}
        open_groups = {}
        for k in mt.postorder():
            node = mt.nodes[k]
            gs = [g for c in node['children'] for g in open_groups.pop(c, [])]
            gs += [by_node[l] for l in node['own'] if l in by_node]
            open_groups[k] = merge_greedy(gs, pmax, lmax, order='savings')
        groups = open_groups.get(mt.root, [])
        seen = {i for g in groups for i in g.terms}   # loads out of the tree (not reachable): alone
        groups += [Group(cell, i) for i in range(len(cell.loads)) if i not in seen]
    else:
        raise ValueError(f"Unknown method {method!r}, expected one of {METHODS}")
    groups = enforce_nmax(groups, nmax)
    return [[cell.loads[i] for i in g.terms] for g in sorted(groups, key=lambda g: -g.est)]


def _steiner_edges(H, terminals, weight):
    if len(terminals) == 1:
        return []
    T = steiner_tree(H, terminals, weight=weight, method='mehlhorn')
    return [tuple(e) for e in T.edges]


def prune_tree(T: nx.Graph, keep) -> nx.Graph:
    """Remove the leaves that are not in `keep`, repeatedly (in place)."""
    leaves = [v for v in T if T.degree(v) <= 1 and v not in keep]
    while leaves:
        v = leaves.pop()
        if v not in T:
            continue
        nbrs = list(T[v])
        T.remove_node(v)
        leaves += [u for u in nbrs if T.degree(u) <= 1 and u not in keep]
    return T


def build_departs(cell: Cell, groups: list) -> list:
    """Edges of the tree of the cell made of one Steiner tree per group (see module docstring)."""
    s, H, w = cell.source, cell.H, cell.weight
    used = set()
    T = nx.Graph()
    T.add_node(s)
    fallback = False
    for g in groups:
        avail = H.subgraph((cell.nodes - used) | {s})
        comp = nx.node_connected_component(avail, s)
        if all(t in comp for t in g):
            edges = _steiner_edges(avail.subgraph(comp), [s] + g, w)
        else:
            edges = _steiner_edges(H, [s] + g, w)
            fallback = True
        T.add_edges_from(edges)
        used |= {v for e in edges for v in e} - {s}
    if fallback or not nx.is_tree(T):
        for u, v in T.edges:
            T[u][v][w] = cell.G[u][v][w]
        T = nx.minimum_spanning_tree(T, weight=w)
        prune_tree(T, set(cell.loads) | {s})
    return sorted((min(u, v), max(u, v)) for u, v in T.edges)


def departs_by_cells(graph, cells: dict, method: str = 'mergetree', tau: float = 0.0,
                     weight: str = "weight") -> dict:
    """
    {source: tree edges} with départs designed by `method` (see module
    docstring) in each cell. Same format as steiner.steiner_trees_by_cells.
    """
    from verify import source_limits
    trees = {}
    for s, nodes in cells.items():
        if s == 'unreachable':
            continue
        lim = source_limits(graph, s)
        cell = Cell(graph.G, nodes, s, weight)
        groups = group_loads(cell, lim['p_max_depart'], lim['l_max_depart'], lim['n_departs_max'],
                             method, tau)
        trees[s] = build_departs(cell, groups) if groups else []
    return trees
