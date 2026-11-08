"""
Decomposition of a network planning instance into smaller subproblems, from the
ToMATo clusters of the load density (clustering.tomato = Algorithm 1 of Chazal,
Guibas, Oudot, Skraba, JACM 2013).

1. Atoms: basins of the density peaks kept by ToMATo with merging parameter tau
   (tau = 0: plain hill-climbing basins). Every node of the graph gets an atom.
   The final filter of the paper (dropping clusters whose peak is below tau) is
   not applied: a decomposition has to cover every node.
2. tau cannot always come from the diagram. Section 3 of the paper: without a
   well-separated diagram, the choice of tau "requires additional
   application-specific information". bootstrap_tau is the bottleneck bootstrap
   of AuToMATo (arXiv:2408.06958), and tells
   whether the diagram has significant peaks at all.
3. The application-specific information is electrical: ToMATo's merging is
   replayed on the atoms, by decreasing saddle height (the order of the
   superlevel-set filtration), with the test on the prominence replaced by a
   test on the union. The sources act as markers: first, only regions holding
   a source absorb atoms, with net power within margin x supply (sum of
   min(capacity, n_departs_max x p_max_depart)); then regions are grouped, up
   to `max_sources` sources each.
4. A region is self-sufficient if it has a source and a load, its net power
   is at most its supply, and each load can reach a source inside the region within
   l_max_depart. Those that are not are merged with a neighbour, youngest
   (lowest) peak first, preferring a neighbour that makes the union
   self-sufficient, then the highest saddle.
5. Each region is a standalone subproblem (induced subgraph, its own sources).
   Regions partition the nodes, so the union of feasible trees is node-disjoint:
   verify.check_solution on the full graph remains the only judge.
"""
import heapq
import json
import math
import time
from pathlib import Path

import networkx as nx
import numpy as np

from voronoi import voronoi_cells_SB
from verify import check_solution, solution_key, source_limits

from . import local_search as ls
from .clustering import loads_of, power_density, tomato
from .partition import best_of, departs_builder, sources_of, standard_inits, steiner_builder
from .persistence import shortest_distances
from .union_find import UnionFind
from .vectorize import bottleneck

BUILDERS = {'departs_savings': departs_builder('savings'),
            'departs_mergetree': departs_builder('mergetree'),
            'steiner': steiner_builder}


class SubGraph:
    """The attributes of nx_graph.GraphData read by verify.py, partition.py and local_search.py."""

    def __init__(self, graph, nodes):
        self.G = graph.G.subgraph(nodes).copy()
        self.nodes = graph.nodes


# ----------------------------------------------------------------------------
# Choice of tau
# ----------------------------------------------------------------------------

def prominences(dgm: np.ndarray) -> np.ndarray:
    """Prominences of the finite points of a superlevel diagram, decreasing."""
    fin = dgm[np.isfinite(dgm[:, 1])]
    return np.sort(fin[:, 0] - fin[:, 1])[::-1]


def bootstrap_tau(G, bandwidth: float, mass: str = "power", weight: str = "weight", alpha: float = 0.35,
                  n_boot: int = 200, seed: int = 0) -> dict:
    """
    Bottleneck bootstrap (AuToMATo): the loads are resampled with replacement,
    the density and its diagram (finite points only) are recomputed, and
    tau = 2 x the (1 - alpha) quantile of the bottleneck distances to the
    original diagram. A peak of prominence >= tau is significant.
    Returns {tau, n_significant (essential peak included), dgm, distances}.
    """
    loads = loads_of(G)
    nodes, index, D = shortest_distances(G, loads, weight, limit=3 * bandwidth)
    K = np.exp(-(D / bandwidth) ** 2)   # exp(-inf) = 0 beyond the cutoff
    w = np.array([abs(G.nodes[l]['conso_ete']) if mass == 'power' else 1.0 for l in loads])

    def diagram(masses):
        r = masses @ K
        return tomato(G, {v: float(r[index[v]]) for v in nodes})['dgm']

    dgm = diagram(w)
    finite = dgm[np.isfinite(dgm[:, 1])]
    rng = np.random.default_rng(seed)
    n = len(loads)
    dists = []
    for _ in range(n_boot):
        counts = rng.multinomial(n, np.full(n, 1 / n))
        d = diagram(w * counts)
        dists.append(bottleneck(d[np.isfinite(d[:, 1])], finite))
    tau = 2 * float(np.quantile(dists, 1 - alpha))
    return {'tau': tau, 'n_significant': 1 + int((prominences(dgm) >= tau).sum()), 'dgm': dgm,
            'distances': dists}


# ----------------------------------------------------------------------------
# Regions
# ----------------------------------------------------------------------------

def fill_unlabelled(G, label: dict) -> dict:
    """Nodes of label -1 (density 0) take the label of an already labelled neighbour, breadth-first."""
    label = dict(label)
    frontier = [v for v in G if label[v] == -1 and any(label[u] != -1 for u in G[v])]
    while frontier:
        for v in frontier:
            label[v] = next(label[u] for u in G[v] if label[u] != -1)
        frontier = [v for v in G if label[v] == -1 and any(label[u] != -1 for u in G[v])]
    return label


def model_size(n_loads: int, n_edges: int, n_nodes: int, n_feeders: int) -> int:
    """Number of variables of plne_reseau.build_model: K (y: L, z: E, x: 2E, f: 2E, w: V)."""
    return n_feeders * (n_loads + 5 * n_edges + n_nodes)


class Regions:
    """Partition of the nodes into regions, with the checks of self-sufficiency."""

    def __init__(self, graph, label: dict, rho: dict, margin: float = 0.85, weight: str = "weight"):
        self.graph, self.G, self.rho, self.margin, self.weight = graph, graph.G, rho, margin, weight
        self.sources = sources_of(self.G)
        self.lims = {s: source_limits(graph, s) for s in self.sources}
        self.d_source = nx.multi_source_dijkstra_path_length(self.G, self.sources, weight=weight)
        self.members, self.saddle = {}, {}
        for v, r in label.items():
            self.members.setdefault(r, set()).add(v)
            self.saddle.setdefault(r, {})
        for u, v in self.G.edges:
            a, b = label[u], label[v]
            if a != b:
                s = min(rho[u], rho[v])
                self.saddle[a][b] = self.saddle[b][a] = max(self.saddle[a].get(b, -math.inf), s)
        self.n_atoms = len(self.members)
        self.peak = {r: max(rho[v] for v in ns) for r, ns in self.members.items()}
        self.merges = []

    def stats(self, r) -> dict:
        return self.stats_of(self.members[r])

    def stats_of(self, ns) -> dict:
        H = self.G.subgraph(ns)
        srcs = [v for v in ns if self.G.nodes[v]['label'] == 'Source']
        loads = [v for v in ns if self.G.nodes[v]['label'] in ('Conso', 'Prod')]
        power = sum(self.G.nodes[l]['conso_ete'] for l in loads)
        supply = sum(min(self.lims[s]['capacity'], self.lims[s]['n_departs_max'] * self.lims[s]['p_max_depart'])
                     for s in srcs)
        feeders = int(sum(self.lims[s]['n_departs_max'] for s in srcs))
        far = []
        if srcs:
            d = nx.multi_source_dijkstra_path_length(H, srcs, weight=self.weight)
            # a load already out of reach in G cannot do better than its own distance in G
            far = [l for l in loads
                   if d.get(l, math.inf) > max(self.lims[srcs[0]]['l_max_depart'], 1.1 * self.d_source[l]) + 1e-9]
        return {'n_nodes': len(ns), 'n_edges': H.number_of_edges(), 'n_sources': len(srcs), 'n_loads': len(loads),
                'power': power, 'supply': supply, 'usage': abs(power) / supply if supply else math.inf,
                'n_far': len(far), 'n_feeders': feeders,
                'n_vars': model_size(len(loads), H.number_of_edges(), len(ns), feeders)}

    def problems(self, r) -> list:
        """Reasons why region r is not self-sufficient (empty list: it is)."""
        return self.problems_of(self.members[r])

    def problems_of(self, ns) -> list:
        st = self.stats_of(ns)
        if st['n_sources'] == 0:
            return ['no source']
        if st['n_loads'] == 0:   # its sources would be stranded
            return ['no load']
        out = []
        if st['usage'] > 1.0 + 1e-9:   # margin is a planning target of the merges, 1 the physical limit
            out.append(f"usage {st['usage']:.2f}")
        if st['n_far']:
            out.append(f"{st['n_far']} loads out of reach")
        return out

    def merge(self, a, b) -> None:
        """Region b is absorbed by region a; saddle(a, c) = max(saddle(a, c), saddle(b, c))."""
        self.members[a] |= self.members.pop(b)
        self.peak[a] = max(self.peak[a], self.peak.pop(b))
        for c, s in self.saddle.pop(b).items():
            del self.saddle[c][b]
            if c != a:
                self.saddle[a][c] = self.saddle[c][a] = max(self.saddle[a].get(c, -math.inf), s)

    def _aggregate(self, r) -> tuple:
        """(n_sources, n_loads, power, supply) of region r."""
        srcs = [v for v in self.members[r] if self.G.nodes[v]['label'] == 'Source']
        loads = [v for v in self.members[r] if self.G.nodes[v]['label'] in ('Conso', 'Prod')]
        return (len(srcs), len(loads), sum(self.G.nodes[l]['conso_ete'] for l in loads),
                sum(min(self.lims[s]['capacity'], self.lims[s]['n_departs_max'] * self.lims[s]['p_max_depart'])
                    for s in srcs))

    def grow_from_sources(self, margins=(0.7, 0.85)) -> None:
        """
        Marker-controlled merging: only a region with a source can absorb a
        region without one, by decreasing saddle, within margin x its supply.
        Without markers, source-free atoms merge with each other up to the
        capacity limit and strand the sources between them. Each value of
        `margins` is one pass: the first one spreads the atoms between the
        sources, the next ones fill them up.
        """
        agg = {r: self._aggregate(r) for r in self.members}
        for margin in margins:
            heap = [(-s, str(a), str(b), a, b) for a in self.members if agg[a][0]
                    for b, s in self.saddle[a].items() if not agg[b][0]]
            heapq.heapify(heap)
            while heap:
                neg, _, _, a, b = heapq.heappop(heap)
                # the saddle of a pair can only rise when a absorbs a neighbour: stale entries are skipped
                if a not in self.members or b not in self.members or agg[b][0] or self.saddle[a].get(b) != -neg:
                    continue
                n_s, n_l, power, supply = (x + y for x, y in zip(agg[a], agg[b]))
                if abs(power) > margin * supply + 1e-9:
                    continue
                self.merge(a, b)
                agg[a] = (n_s, n_l, power, supply)
                del agg[b]
                self.merges.append({'region': str(b), 'into': str(a), 'saddle': -neg, 'reason': ['marker']})
                for c, s in self.saddle[a].items():
                    if not agg[c][0]:
                        heapq.heappush(heap, (-s, str(a), str(c), a, c))

    def merge_along_persistence(self, max_sources: float = 1, max_loads: float = math.inf) -> None:
        """
        ToMATo's merging replayed on the regions, by decreasing saddle, a merge
        being allowed if the union has at most `max_sources` sources and at
        most `max_loads` loads, and its net power is within margin x its supply
        (x the largest supply of one source if it has no source yet). A
        rejected merge stays rejected (regions only grow), so one pass is enough.
        """
        unit = max((min(l['capacity'], l['n_departs_max'] * l['p_max_depart']) for l in self.lims.values()), default=0.0)
        agg = {r: self._aggregate(r) for r in self.members}
        pairs = sorted(((s, a, b) for a, nbrs in self.saddle.items() for b, s in nbrs.items() if str(a) < str(b)),
                       key=lambda p: (-p[0], str(p[1]), str(p[2])))
        uf = UnionFind(self.members)
        for saddle, a, b in pairs:
            ra, rb = uf.find(a), uf.find(b)
            if ra == rb:
                continue
            n_s, n_l, power, supply = (x + y for x, y in zip(agg[ra], agg[rb]))
            if n_s > max_sources or n_l > max_loads or abs(power) > self.margin * (supply if n_s else unit) + 1e-9:
                continue
            root = uf.union(ra, rb)
            other = rb if root == ra else ra
            self.merge(root, other)
            agg[root] = (n_s, n_l, power, supply)
            self.merges.append({'region': str(other), 'into': str(root), 'saddle': saddle, 'reason': ['persistence']})

    def merge_until_sufficient(self, max_loads: float = math.inf) -> None:
        """
        The youngest region that is not self-sufficient is merged with a
        neighbour, until none is left. Among the neighbours, in this order of
        preference: those whose union is self-sufficient and has at most
        `max_loads` loads, those within `max_loads`, then all of them. The
        highest saddle decides within the preferred group (the merge ToMATo
        would do first); without the preference, merges snowball into the
        densest region.
        """
        while len(self.members) > 1:
            bad = [r for r in self.members if self.problems(r)]
            if not bad:
                return
            r = min(bad, key=lambda x: (self.peak[x], str(x)))
            if not self.saddle[r]:   # disconnected graph
                return
            n_r = self.stats(r)['n_loads']
            within = [c for c in self.saddle[r] if n_r + self.stats(c)['n_loads'] <= max_loads]
            good = [c for c in within if not self.problems_of(self.members[r] | self.members[c])]
            c = max(good or within or list(self.saddle[r]), key=lambda x: (self.saddle[r][x], str(x)))
            self.merges.append({'region': str(r), 'into': str(c), 'saddle': self.saddle[r][c],
                                'reason': self.problems(r)})
            self.merge(c, r)


def decompose(graph, bandwidth: float, mass: str = "power", tau: float = 0.0, max_sources: float = 1,
              margin: float = 0.85, max_loads: float = math.inf, weight: str = "weight") -> dict:
    """
    Regions from the ToMATo atoms of the density (see module docstring).
    Returns {regions: {id: nodes}, stats: {id: dict}, n_atoms, tau, rho, dgm, peaks, merges,
    unresolved}. Regions are numbered by decreasing peak height.
    """
    G = graph.G
    rho, h = power_density(G, bandwidth, weight, mass)
    full = tomato(G, rho, min_density=0.0)
    atoms = tomato(G, rho, tau, min_density=0.0) if tau < math.inf else full
    regs = Regions(graph, fill_unlabelled(G, atoms['label']), rho, margin, weight)
    regs.grow_from_sources((min(0.7, margin), margin))
    regs.merge_along_persistence(max_sources, max_loads)
    regs.merge_until_sufficient(max_loads)
    return _result(regs, tau, h, full['dgm'], atoms['peaks'])


def voronoi_control(graph, max_sources: float = 1, margin: float = 0.85, bandwidth: float = None,
                    max_loads: float = math.inf, weight: str = "weight") -> dict:
    """Same merging rules from the plain Voronoi cells of the sources instead of the ToMATo atoms."""
    G = graph.G
    rho, h = power_density(G, bandwidth, weight)
    _, cells, _ = voronoi_cells_SB(G, sources_of(G), weight=weight, weighting='none')
    label = {v: s for s, ns in cells.items() for v in ns}
    regs = Regions(graph, label, rho, margin, weight)
    regs.merge_along_persistence(max_sources, max_loads)
    regs.merge_until_sufficient(max_loads)
    return _result(regs, None, h, None, [])


def _result(regs: Regions, tau, h, dgm, peaks) -> dict:
    order = sorted(regs.members, key=lambda r: (-regs.peak[r], str(r)))
    return {'regions': {i: regs.members[r] for i, r in enumerate(order)},
            'stats': {i: regs.stats(r) for i, r in enumerate(order)},
            'unresolved': {i: regs.problems(r) for i, r in enumerate(order) if regs.problems(r)},
            'n_atoms': regs.n_atoms, 'tau': tau, 'bandwidth': h, 'rho': regs.rho, 'dgm': dgm, 'peaks': peaks,
            'merges': regs.merges}


# ----------------------------------------------------------------------------
# Subproblems
# ----------------------------------------------------------------------------

def solve_subproblem(sub, n_iter: int = 60, ls_time: float = 20.0) -> dict:
    """
    Best of the tree builders with the variable-distance loop of partition.py,
    then local search (kept only if it improves the certificate).
    Returns {trees, cert, builder, init, time_build, time_ls}.
    """
    t0 = time.perf_counter()
    inits = standard_inits(sub)
    if len(sources_of(sub.G)) == 1:   # the weights have no effect with one source
        inits = {'none': inits['none']}
    best, name = None, None
    for builder_name, builder in BUILDERS.items():
        run = best_of(sub, inits, tree_builder=builder, n_iter=n_iter, grow=True)
        if best is None or solution_key(run['cert']) < solution_key(best['cert']):
            best, name = run, builder_name
    t1 = time.perf_counter()
    trees, cert = best['trees'], best['cert']
    new, _ = ls.local_search(sub, trees, ls.candidate_insert_edges(sub), ls.candidate_relocations(sub),
                             time_limit=ls_time)
    new_cert = check_solution(sub, new)
    if solution_key(new_cert) < solution_key(cert):
        trees, cert = new, new_cert
    return {'trees': trees, 'cert': cert, 'builder': name, 'init': best.get('init'),
            'time_build': t1 - t0, 'time_ls': time.perf_counter() - t1}


def solve_regions(graph, regions: dict, n_iter: int = 60, ls_time: float = 20.0, polish_time: float = 60.0) -> dict:
    """
    Each region solved alone, trees put together and judged on the full graph.
    Returns {trees, cert, polished_trees, polished_cert, regions: {id: sub-result without trees},
    time_sum, time_max, time_polish}. Polish = local search on the full graph (moves across regions).
    """
    trees, per_region = {}, {}
    for i, nodes in regions.items():
        sub = SubGraph(graph, nodes)
        if not sources_of(sub.G):   # only transit nodes: nothing to serve
            per_region[i] = {'time': 0.0, 'J': 0.0, 'violations': {}, 'builder': None}
            continue
        res = solve_subproblem(sub, n_iter, ls_time)
        trees.update(res['trees'])
        per_region[i] = {'time': res['time_build'] + res['time_ls'], 'J': res['cert']['J'],
                         'violations': res['cert']['violations'], 'builder': res['builder']}
    cert = check_solution(graph, trees)
    t0 = time.perf_counter()
    new, _ = ls.local_search(graph, trees, ls.candidate_insert_edges(graph), ls.candidate_relocations(graph),
                             time_limit=polish_time)
    polished_cert = check_solution(graph, new)
    if solution_key(polished_cert) >= solution_key(cert):
        new, polished_cert = trees, cert
    times = [r['time'] for r in per_region.values()]
    return {'trees': trees, 'cert': cert, 'polished_trees': new, 'polished_cert': polished_cert,
            'regions': per_region, 'time_sum': sum(times), 'time_max': max(times),
            'time_polish': time.perf_counter() - t0}


def export_subinstances(source_json, regions: dict, out_dir, elec_defaults: bool = True) -> list:
    """
    One JSON per region, in the format of data/instances (readable by
    src.instances.load_instance), nodes renumbered, 'original_id' = node of the
    full instance. Returns the paths.
    """
    with open(source_json, encoding='utf-8') as f:
        data = json.load(f)
    elec = data.get('electrotechnique', {})
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for i, nodes in regions.items():
        new_id = {v: k for k, v in enumerate(sorted(nodes))}
        out_nodes = []
        for v in sorted(nodes):
            n = {k: x for k, x in data['nodes'][v].items() if k not in ('neighbors', 'connected_edges')}
            n.update(original_id=v, id_noeud=new_id[v])
            n.setdefault('demande_hivers', n['demande'])
            if n['type'] == 'Source' and elec_defaults:
                n.setdefault('p_max_depart', elec.get('puissance_max_ligne'))
                n.setdefault('l_max_depart', elec.get('longueur_max_depart'))
            out_nodes.append(n)
        out_edges = [dict(e, start_node_idx=new_id[e['start_node_idx']], end_node_idx=new_id[e['end_node_idx']])
                     for e in data['edges'] if e['start_node_idx'] in new_id and e['end_node_idx'] in new_id]
        path = out_dir / f"region_{i:02d}.json"
        with open(path, 'w', encoding='utf-8') as f:
            json.dump({'nodes': out_nodes, 'edges': out_edges, 'electrotechnique': elec}, f)
        paths.append(path)
    return paths
