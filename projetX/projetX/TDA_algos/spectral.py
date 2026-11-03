"""
Spectral tools on the instance graph. Edges are conductances 1 / length
(electrical analogy: a line of length l has resistance l).

  spectral_embedding     : first eigenvectors of the normalised Laplacian
                           (random-walk scaling D^-1/2 U), and its eigenvalues
  spectral_edge_costs    : edge attribute length * (1 + beta * jump of the
                           embedding across the edge): edges across bottlenecks
                           of the graph (valleys, sparse corridors) cost more,
                           so Voronoi cells computed with it cut there. Steiner
                           trees and J still use the true lengths.
  anchored_spectral_cells: multi-way spectral clustering anchored on the sources
                           (each node goes to the source closest in the
                           row-normalised embedding), then connectivity repair
  effective_resistance   : R(s, v) for a few sources s and all v (grounded
                           Laplacian, sparse LU). R(s, v) <= d(s, v), with
                           equality when there is a single route: R / d measures
                           the redundancy of the routes to a load
  kron_reduction         : Schur complement of the Laplacian onto a set of nodes
"""
import heapq
import math

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh, splu

MIN_LENGTH = 1e-6
ATTR = 'w_spectral'


def conductance_matrix(G, weight: str = "weight", nodes=None):
    """(nodes, index, A): symmetric sparse matrix of the conductances 1 / length."""
    nodes = list(G) if nodes is None else list(nodes)
    index = {v: i for i, v in enumerate(nodes)}
    rows, cols, vals = [], [], []
    for u, v, w in G.edges(nodes, data=weight, default=1.0):
        if u in index and v in index:
            c = 1.0 / max(float(w), MIN_LENGTH)
            rows += [index[u], index[v]]
            cols += [index[v], index[u]]
            vals += [c, c]
    n = len(nodes)
    return nodes, index, coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()


def laplacian(A):
    deg = np.asarray(A.sum(axis=1)).ravel()
    return diags(deg) - A, deg


def spectral_embedding(G, k: int, weight: str = "weight", nodes=None) -> dict:
    """
    {nodes, index, eigenvalues: k+1 smallest of the normalised Laplacian,
    X: (n, k) embedding without the trivial eigenvector}.
    """
    nodes, index, A = conductance_matrix(G, weight, nodes)
    n = len(nodes)
    _, deg = laplacian(A)
    dm12 = diags(1.0 / np.sqrt(np.maximum(deg, 1e-12)))
    M = dm12 @ A @ dm12                     # eigenvalues 1 - lambda(L_norm), in [-1, 1]
    k_eff = max(1, min(k + 1, n - 2))
    if n <= 600:
        vals, vecs = np.linalg.eigh(M.toarray())
        vals, vecs = vals[::-1][:k_eff], vecs[:, ::-1][:, :k_eff]
    else:
        vals, vecs = eigsh(M, k=k_eff, which='LA', tol=1e-8)
        order = np.argsort(-vals)
        vals, vecs = vals[order], vecs[:, order]
    X = dm12 @ vecs[:, 1:]
    return {'nodes': nodes, 'index': index, 'eigenvalues': 1 - vals, 'X': np.asarray(X)}


def spectral_edge_costs(G, k: int = None, beta: float = 1.0, weight: str = "weight",
                        attr: str = ATTR, cap: float = 10.0) -> str:
    """
    Set the edge attribute `attr` = length * (1 + beta * min(jump / median jump, cap)),
    jump = ||X_u - X_v|| in the spectral embedding of dimension k (default: number of sources).
    Returns `attr` (to be passed as `weight` to the Voronoi functions).
    """
    if k is None:
        k = max(2, sum(lab == 'Source' for _, lab in G.nodes('label')))
    emb = spectral_embedding(G, k, weight)
    X, index = emb['X'], emb['index']
    edges = list(G.edges)
    jump = np.array([np.linalg.norm(X[index[u]] - X[index[v]]) for u, v in edges])
    scale = np.median(jump[jump > 0]) if (jump > 0).any() else 1.0
    for (u, v), j in zip(edges, jump):
        G[u][v][attr] = G[u][v][weight] * (1 + beta * min(j / scale, cap))
    return attr


def anchored_spectral_cells(G, k: int = None, weight: str = "weight") -> dict:
    """
    {source: set of nodes}: each node goes to the source nearest in the
    row-normalised spectral embedding (dimension k, default number of
    sources); then each cell keeps the connected component of its source and
    the other nodes are given to the nearest cell (Dijkstra on the lengths).
    """
    sources = [n for n, lab in G.nodes('label') if lab == 'Source']
    if k is None:
        k = max(2, len(sources))
    emb = spectral_embedding(G, k, weight)
    X, index = emb['X'], emb['index']
    Xn = X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    S = Xn[[index[s] for s in sources]]
    nearest = np.argmin(((Xn[:, None, :] - S[None, :, :]) ** 2).sum(-1), axis=1)
    owner = {v: sources[nearest[index[v]]] for v in G}
    for s in sources:
        owner[s] = s
    # connectivity repair: keep the component of the source inside its cell
    final = {}
    for s in sources:
        members = {v for v, o in owner.items() if o == s}
        seen, stack = {s}, [s]
        while stack:
            v = stack.pop()
            for u in G[v]:
                if u in members and u not in seen:
                    seen.add(u)
                    stack.append(u)
        for v in seen:
            final[v] = s
    dist = {v: 0.0 for v in final}
    heap = [(0.0, str(v), v) for v in final]
    heapq.heapify(heap)
    while heap:
        d, _, v = heapq.heappop(heap)
        if d > dist.get(v, math.inf):
            continue
        for u, data in G[v].items():
            nd = d + data[weight]
            if G.nodes[u]['label'] != 'Source' and nd < dist.get(u, math.inf):
                dist[u] = nd           # the seeds (distance 0) are never overwritten
                final[u] = final[v]
                heapq.heappush(heap, (nd, str(u), u))
    cells = {s: set() for s in sources}
    for v, s in final.items():
        cells[s].add(v)
    missing = set(G) - set(final)
    if missing:
        cells['unreachable'] = missing
    return cells


def effective_resistance(G, sources, targets=None, weight: str = "weight") -> dict:
    """
    {s: {v: R(s, v)}} for v in `targets` (all nodes by default), with
    R(u, v) = Gi[u, u] + Gi[v, v] - 2 Gi[u, v], Gi the inverse of the
    Laplacian grounded at one node. The graph must be connected.
    """
    nodes, index, A = conductance_matrix(G, weight)
    L, _ = laplacian(A)
    g = len(nodes) - 1                               # grounded node
    lu = splu(L[:g, :g].tocsc())
    targets = list(G) if targets is None else list(targets)

    def column(i):
        if i == g:
            return np.zeros(len(nodes))
        b = np.zeros(g)
        b[i] = 1.0
        return np.append(lu.solve(b), 0.0)

    diag = {}
    for v in set(targets) | set(sources):
        i = index[v]
        diag[v] = column(i)[i] if i != g else 0.0
    out = {}
    for s in sources:
        cs = column(index[s])
        out[s] = {v: float(diag[s] + diag[v] - 2 * cs[index[v]]) for v in targets}
    return out


def kron_reduction(G, keep, weight: str = "weight") -> tuple:
    """(keep, K): Schur complement L_kk - L_kr L_rr^-1 L_rk of the Laplacian onto `keep` (dense)."""
    nodes, index, A = conductance_matrix(G, weight)
    L, _ = laplacian(A)
    L = L.tocsr()
    keep = list(keep)
    k_idx = np.array([index[v] for v in keep])
    r_idx = np.array(sorted(set(range(len(nodes))) - set(k_idx.tolist())))
    Lkk = L[k_idx][:, k_idx].toarray()
    if not len(r_idx):
        return keep, Lkk
    Lkr = L[k_idx][:, r_idx]
    lu = splu(L[r_idx][:, r_idx].tocsc())
    Y = lu.solve(Lkr.T.toarray())
    return keep, Lkk - Lkr @ Y
