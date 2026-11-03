"""
Local search on a solution {source: tree edges}, every move being accepted
only if it improves the score of the sources it touches, lexicographically
(number of violated constraints, total relative excess, length) - the order
of verify.solution_key, source by source.

Moves
  insert   : add an edge e = (u, v) of the graph between two nodes of the same
             tree; it closes a cycle in the tree. Remove the longest key path
             of that cycle (a path whose inner nodes are degree-2 transit
             nodes) if it is longer than e.
             Candidates: all such edges, or, with an H1 filter, only the
             edges that close a persistent cycle of the planarised graph
             (persistence.h1_planar): loops around holes of the network, the
             real alternative routes, instead of the mesh triangles.
  relocate : detach a leaf load and its branch (up to the first key node) and
             reconnect it, by a shortest path through free nodes, to the
             closest node of each tree (its own tree included: attaching to
             the source directly opens a new départ). Candidates: all leaf
             loads, or, with a margin filter, the loads whose two nearest
             sources are almost equidistant (the cut locus of the Voronoi
             cells, where the partition is the least certain).
"""
import heapq
import math
import time

import networkx as nx
import numpy as np

from verify import departs_of_tree, excess, source_limits

from .persistence import shortest_distances

LOAD_LABELS = ('Conso', 'Prod')
EPS = 1e-9


def source_score(graph, s, T: nx.Graph) -> tuple:
    """(violations, excess, length) of the tree T of source s, as in verify.check_solution."""
    G = graph.G
    lim = source_limits(graph, s)
    edges = list(T.edges)
    departs = departs_of_tree(G, s, edges) if edges else []
    n_viol, exc, length = 0, 0.0, 0.0
    for d in departs:
        e_L, e_P = excess(d['length'], lim['l_max_depart']), excess(abs(d['power']), lim['p_max_depart'])
        n_viol += (e_L > 0) + (e_P > 0)
        exc += e_L + e_P
        length += d['length']
    power = sum(G.nodes[n]['conso_ete'] for n in T if n != s)
    e_src, e_n = excess(abs(power), lim['capacity']), excess(len(departs), lim['n_departs_max'])
    n_viol += (e_src > 0) + (e_n > 0)
    return n_viol, exc + e_src + e_n, length


def _add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _better(new, old) -> bool:
    """new < old lexicographically on (violations, excess, length), with tolerances."""
    if new[0] != old[0]:
        return new[0] < old[0]
    if abs(new[1] - old[1]) > 1e-9:
        return new[1] < old[1]
    return new[2] < old[2] - 1e-9


class Forest:

    def __init__(self, graph, trees: dict):
        self.graph, self.G = graph, graph.G
        self.loads = {v for v, lab in self.G.nodes('label') if lab in LOAD_LABELS}
        self.sources = {v for v, lab in self.G.nodes('label') if lab == 'Source'}
        self.T = {}
        for s, edges in trees.items():
            T = nx.Graph()
            T.add_node(s)
            T.add_edges_from(edges)
            self.T[s] = T
        self.owner = {v: s for s, T in self.T.items() for v in T}
        self.score = {s: source_score(graph, s, T) for s, T in self.T.items()}

    def trees(self) -> dict:
        return {s: sorted((min(u, v), max(u, v)) for u, v in T.edges) for s, T in self.T.items()}

    def is_key(self, T, v) -> bool:
        return v in self.loads or v in self.sources or T.degree(v) != 2

    def _commit(self, s, T, score):
        for v in self.T[s]:
            if self.owner.get(v) == s:
                del self.owner[v]
        self.T[s] = T
        for v in T:
            self.owner[v] = s
        self.score[s] = score

    # ------------------------------------------------------------------ insert
    def try_insert(self, u, v) -> bool:
        s = self.owner.get(u)
        if s is None or self.owner.get(v) != s or self.T[s].has_edge(u, v):
            return False
        T = self.T[s]
        path = nx.shortest_path(T, u, v)
        keys = [i for i, x in enumerate(path) if i in (0, len(path) - 1) or self.is_key(T, x)]
        best, best_len = None, 0.0
        for a, b in zip(keys, keys[1:]):
            seg = path[a:b + 1]
            length = sum(self.G[x][y]['weight'] for x, y in zip(seg, seg[1:]))
            if length > best_len:
                best, best_len = seg, length
        if best is None or best_len <= self.G[u][v]['weight'] + EPS:
            return False
        T2 = T.copy()
        T2.remove_edges_from(zip(best, best[1:]))
        T2.remove_nodes_from([x for x in best[1:-1]])
        T2.add_edge(u, v)
        if not nx.is_tree(T2):
            return False
        score = source_score(self.graph, s, T2)
        if _better(score, self.score[s]):
            self._commit(s, T2, score)
            return True
        return False

    # ---------------------------------------------------------------- relocate
    def leaf_branch(self, T, v) -> list:
        """[v, ..., a]: the leaf v up to its first key node a (excluded from removal)."""
        branch = [v]
        prev, cur = None, v
        while True:
            nxt = [x for x in T[cur] if x != prev]
            if len(nxt) != 1:
                break
            prev, cur = cur, nxt[0]
            branch.append(cur)
            if self.is_key(T, cur):
                break
        return branch

    def try_relocate(self, v, cutoff=math.inf) -> bool:
        s = self.owner.get(v)
        if s is None or v not in self.loads or self.T[s].degree(v) != 1:
            return False
        T = self.T[s]
        branch = self.leaf_branch(T, v)
        removed = set(branch[:-1])
        T_s = T.copy()
        T_s.remove_nodes_from(removed)
        # shortest paths from v through free nodes to the closest node of each tree
        free = lambda x: (x not in self.owner or x in removed) and x not in self.sources
        dist, prev, hit = {v: 0.0}, {}, {}
        heap = [(0.0, str(v), v)]
        while heap and len(hit) < len(self.T):
            d, _, x = heapq.heappop(heap)
            if d > dist[x] or d > cutoff:
                continue
            o = self.owner.get(x)
            if x != v and x not in removed and o is not None:
                hit.setdefault(o, x)
                continue          # never go through a tree
            for y, data in self.G[x].items():
                if not (free(y) or (self.owner.get(y) is not None and y not in removed)):
                    continue
                nd = d + data['weight']
                if nd < dist.get(y, math.inf):
                    dist[y], prev[y] = nd, x
                    heapq.heappush(heap, (nd, str(y), y))
        old = self.score[s]
        base_s = source_score(self.graph, s, T_s)
        for t, a in hit.items():
            path = [a]
            while path[-1] != v:
                path.append(prev[path[-1]])
            if t == s:
                T_t = T_s.copy()
                T_t.add_edges_from(zip(path, path[1:]))
                if not nx.is_tree(T_t):
                    continue
                score = source_score(self.graph, s, T_t)
                if _better(score, old):
                    self._commit(s, T_t, score)
                    return True
            else:
                T_t = self.T[t].copy()
                T_t.add_edges_from(zip(path, path[1:]))
                if not nx.is_tree(T_t):
                    continue
                score_t = source_score(self.graph, t, T_t)
                if _better(_add(base_s, score_t), _add(old, self.score[t])):
                    self._commit(s, T_s, base_s)
                    self._commit(t, T_t, score_t)
                    return True
        return False


def boundary_margin(graph, weight: str = "weight") -> dict:
    """{load: distance to its 2nd nearest source - distance to its nearest source}."""
    G = graph.G
    sources = [n for n, lab in G.nodes('label') if lab == 'Source']
    loads = [n for n, lab in G.nodes('label') if lab in LOAD_LABELS]
    _, index, D = shortest_distances(G, sources, weight)
    out = {}
    for l in loads:
        d = np.sort(D[:, index[l]])
        out[l] = float(d[1] - d[0]) if len(d) > 1 else math.inf
    return out


def candidate_insert_edges(graph, h1: dict = None, quantile: float = None) -> list:
    """
    Edges for the insert move: all the edges if h1 is None; else the edges
    closing a cycle of the planarised graph whose persistence is in the top
    `quantile` fraction (default: above the largest gap of the persistences).
    """
    if h1 is None:
        return list(graph.G.edges)
    from .vectorize import largest_gap
    pers = h1['edge_persistence']
    if not pers:
        return []
    values = np.array(list(pers.values()))
    thr = np.quantile(values, 1 - quantile) if quantile is not None else largest_gap(values)[0]
    return [e for e, p in pers.items() if p >= thr]


def candidate_relocations(graph, quantile: float = None) -> list:
    """Loads to relocate: all if quantile is None, else those of smallest margin (fraction `quantile`)."""
    margin = boundary_margin(graph)
    if quantile is None:
        return list(margin)
    thr = np.quantile(list(margin.values()), quantile)
    return [l for l, m in margin.items() if m <= thr]


def local_search(graph, trees: dict, insert_edges=None, relocate_loads=None, max_rounds: int = 20,
                 time_limit: float = 120.0, cutoff: float = math.inf) -> tuple:
    """
    Apply the moves until no move improves or the limits are reached.
    `insert_edges` / `relocate_loads`: candidates (None: move not used).
    Returns (trees, stats {rounds, tried, accepted_insert, accepted_relocate, time}).
    """
    F = Forest(graph, trees)
    t0 = time.time()
    stats = {'rounds': 0, 'tried': 0, 'accepted_insert': 0, 'accepted_relocate': 0}
    for _ in range(max_rounds):
        stats['rounds'] += 1
        improved = False
        for u, v in insert_edges or []:
            stats['tried'] += 1
            if F.try_insert(u, v):
                stats['accepted_insert'] += 1
                improved = True
            if time.time() - t0 > time_limit:
                break
        for l in relocate_loads or []:
            stats['tried'] += 1
            if F.try_relocate(l, cutoff):
                stats['accepted_relocate'] += 1
                improved = True
            if time.time() - t0 > time_limit:
                break
        if not improved or time.time() - t0 > time_limit:
            break
    stats['time'] = time.time() - t0
    return F.trees(), stats
