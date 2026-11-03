"""
Persistence-based clustering of the loads (ToMATo, Chazal, Guibas, Oudot,
Skraba 2013) on the graph.

1. Power density on the nodes of the graph:
     rho(v) = sum over the loads l of |demande(l)| * exp(-(d(v, l) / h)^2)
   with d the shortest-path distance and h a bandwidth (default: twice the
   median distance between a load and its nearest other load).
2. Superlevel-set H0 of rho, swept from the densest node: a node with no
   denser neighbour starts a cluster (a peak); a node joins the cluster of its
   densest neighbour (discrete gradient), and two clusters meeting at a node v
   are merged if the lower peak has prominence rho(peak) - rho(v) < tau.
3. tau = inf gives the full persistence diagram (peak, merge density); tau is
   then chosen in the largest gap of the prominences (see vectorize.largest_gap),
   so that only the significant clusters (towns, villages, valleys) are kept.

Clusters are used to initialise the source weights of the Voronoi cells
(partition.tomato_init_weights): a cluster should not be cut by a cell boundary.
"""
import math

import numpy as np

from .persistence import shortest_distances
from .union_find import UnionFind
from .vectorize import largest_gap

LOAD_LABELS = ('Conso', 'Prod')


def loads_of(G) -> list:
    return [n for n, lab in G.nodes('label') if lab in LOAD_LABELS]


def power_density(G, bandwidth: float = None, weight: str = "weight", mass: str = "power") -> tuple:
    """
    ({node: rho(node)}, bandwidth). `mass` = 'power' (|demande|) or 'count'
    (every load weighs 1).
    """
    loads = loads_of(G)
    if not loads:
        return {v: 0.0 for v in G}, bandwidth or 1.0
    if bandwidth is None:
        _, index, D = shortest_distances(G, loads, weight)
        L = D[:, [index[l] for l in loads]]
        np.fill_diagonal(L, np.inf)
        nn = L.min(axis=1)
        nn = nn[np.isfinite(nn)]
        bandwidth = 2 * float(np.median(nn)) if len(nn) else 1.0
        bandwidth = max(bandwidth, 1e-6)
    nodes, index, D = shortest_distances(G, loads, weight, limit=3 * bandwidth)
    w = np.array([abs(G.nodes[l]['conso_ete']) if mass == 'power' else 1.0 for l in loads])
    K = np.exp(-(D / bandwidth) ** 2)   # exp(-inf) = 0 beyond the cutoff
    rho = w @ K
    return {v: float(rho[index[v]]) for v in nodes}, bandwidth


def tomato(G, rho: dict, tau: float = math.inf, min_density: float = 1e-9) -> dict:
    """
    ToMATo with merging threshold `tau` on the nodes of density > min_density * max rho.

    Returns
    -------
    dict with
      label    : {node: cluster id (= its peak node), -1 for the nodes below the density threshold}
      dgm      : (k, 2) array (rho(peak), rho at merge), inf merge for the surviving clusters
      peaks    : list of the peaks of the surviving clusters
      merge_at : {peak of a merged cluster: rho at its merge}
    """
    top = max(rho.values()) if rho else 0.0
    nodes = [v for v in G if rho[v] > min_density * top]
    order = sorted(nodes, key=lambda v: (-rho[v], str(v)))
    rank = {v: i for i, v in enumerate(order)}
    uf = UnionFind()
    peak = {}      # {root: peak node}
    dgm, merge_at = [], {}
    for v in order:
        higher = [u for u in G[v] if u in rank and rank[u] < rank[v]]
        uf.add(v)
        if not higher:
            peak[v] = v
            continue
        rg = uf.find(min(higher, key=rank.get))   # cluster of the densest neighbour
        p = peak[rg]
        peak[uf.union(v, rg)] = p
        for u in higher:
            ru, rv = uf.find(u), uf.find(v)
            if ru == rv:
                continue
            pu, pv = peak[ru], peak[rv]
            low, high = (pu, pv) if (rho[pu], str(pu)) < (rho[pv], str(pv)) else (pv, pu)
            if rho[low] - rho[v] < tau:
                dgm.append((rho[low], rho[v]))
                merge_at[low] = rho[v]
                peak[uf.union(ru, rv)] = high
    label = {v: -1 for v in G}
    peaks = []
    for root, members in uf.groups().items():
        p = peak[root]
        peaks.append(p)
        dgm.append((rho[p], math.inf))
        for v in members:
            label[v] = p
    return {'label': label, 'dgm': np.array(dgm, float).reshape(-1, 2), 'peaks': peaks,
            'merge_at': merge_at}


def auto_tau(dgm: np.ndarray) -> float:
    """Merging threshold in the largest gap of the prominences of the full diagram."""
    fin = dgm[np.isfinite(dgm[:, 1])]
    prom = fin[:, 0] - fin[:, 1]
    thr, _ = largest_gap(prom)
    return thr


def cluster_loads(G, bandwidth: float = None, tau: float = None, weight: str = "weight",
                  mass: str = "power") -> dict:
    """
    ToMATo clusters of the loads, tau chosen automatically if None.

    Returns dict with label, dgm (full diagram, tau = inf), peaks, tau, bandwidth,
    rho, and clusters = {peak: {'loads': [...], 'power': net demand,
    'gross': sum of the positive demands, 'prominence'}}.
    """
    rho, h = power_density(G, bandwidth, weight, mass)
    full = tomato(G, rho)
    if tau is None:
        tau = auto_tau(full['dgm'])
    res = tomato(G, rho, tau) if math.isfinite(tau) else full
    clusters = {}
    for l in loads_of(G):
        p = res['label'][l]
        c = clusters.setdefault(p, {'loads': [], 'power': 0.0, 'gross': 0.0})
        c['loads'].append(l)
        d = G.nodes[l]['conso_ete']
        c['power'] += d
        c['gross'] += max(d, 0.0)
    for p, c in clusters.items():   # prominence in the full diagram (inf for the global peak)
        c['prominence'] = rho[p] - full['merge_at'].get(p, -math.inf) if p != -1 else 0.0
    return {**res, 'full_dgm': full['dgm'], 'tau': tau, 'bandwidth': h, 'rho': rho,
            'clusters': clusters}
