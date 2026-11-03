"""
Persistent homology of the filtrations used on the instance graphs.

A diagram is an (n, 2) numpy array of (birth, death) pairs, death = inf for
essential classes.

  h0_lower_star  : H0 of the sublevel sets of a vertex function f (an edge
                   enters at max(f(u), f(v))). h0_superlevel for the superlevel
                   sets (births are maxima, deaths are saddles, birth >= death).
  terminal_h0    : H0 of the terminals under the shortest-path metric of the
                   graph (Vietoris-Rips H0 = single linkage). Groups of terminals
                   can be born already merged (e.g. all the sources, contracted
                   into one root). Computed in O(m log n) with Mehlhorn's
                   theorem: the MST of the graph G1' (one candidate edge per
                   graph edge joining two Voronoi regions of the terminals) is a
                   MST of the metric closure. Death times = MST edge lengths.
  h1_planar      : H1 of the 2-complex made of the graph and the faces of its
                   straight-line drawing. The drawings are not planar (source
                   links cross Gabriel edges), so the drawing is planarised
                   first (a virtual vertex at each crossing). An edge enters at
                   its length, a face at max(half perimeter, its longest edge).
                   Pairs are found by duality: primal Kruskal finds the edges
                   that kill H0, the others create cycles, and a union-find on
                   the faces in decreasing order (elder rule) pairs each
                   cycle-creating edge with the face that fills its cycle.
                   Persistence = half perimeter of the hole - length of the edge
                   closing it: large for loops around holes of the network
                   (alternative routes), ~0 for the mesh triangles.
"""
import math

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, triu
from scipy.sparse.csgraph import dijkstra

from .union_find import UnionFind

INF = math.inf
EPS_LENGTH = 1e-12   # zero-length edges would be dropped by scipy sparse matrices


# ----------------------------------------------------------------------------
# Graph to sparse matrix
# ----------------------------------------------------------------------------

def graph_csr(G, weight: str = "weight", nodes=None):
    """
    (nodes, index, M): list of the nodes, {node: row}, symmetric CSR matrix of
    the edge weights of the subgraph induced by `nodes` (all of G by default).
    """
    nodes = list(G) if nodes is None else list(nodes)
    index = {v: i for i, v in enumerate(nodes)}
    rows, cols, vals = [], [], []
    for u, v, w in G.edges(nodes, data=weight, default=1.0):
        if u in index and v in index:
            w = max(float(w), EPS_LENGTH)
            rows += [index[u], index[v]]
            cols += [index[v], index[u]]
            vals += [w, w]
    n = len(nodes)
    M = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    M.sum_duplicates()
    return nodes, index, M


def shortest_distances(G, sources, weight: str = "weight", nodes=None, limit=np.inf):
    """
    (nodes, index, D): D[i, j] = distance from sources[i] to nodes[j] in the
    subgraph induced by `nodes` (inf if unreachable or farther than `limit`).
    """
    nodes, index, M = graph_csr(G, weight, nodes)
    D = dijkstra(M, directed=False, indices=[index[s] for s in sources], limit=limit)
    return nodes, index, np.atleast_2d(D)


# ----------------------------------------------------------------------------
# H0 of vertex functions
# ----------------------------------------------------------------------------

def h0_lower_star(G, f: dict, nodes=None) -> dict:
    """
    H0 of the sublevel sets of f on the subgraph induced by `nodes`.

    Returns {dgm: (k, 2) array of (birth, death), pairs: list of
    (vertex of birth, vertex of death or None)}. Pairs of zero persistence
    (a vertex joining an existing component) are not reported.
    """
    nodes = list(G) if nodes is None else list(nodes)
    inside = set(nodes)
    order = sorted(nodes, key=lambda v: (f[v], str(v)))
    uf = UnionFind()
    born = {}    # {root: vertex of birth (minimum of the component)}
    dgm, pairs = [], []
    for v in order:
        uf.add(v)
        born[v] = v
        for u in G[v]:
            if u not in inside or u not in uf:
                continue
            ru, rv = uf.find(u), uf.find(v)
            if ru == rv:
                continue
            bu, bv = born[ru], born[rv]
            old, young = (bu, bv) if (f[bu], str(bu)) <= (f[bv], str(bv)) else (bv, bu)
            if f[young] < f[v]:
                dgm.append((f[young], f[v]))
                pairs.append((young, v))
            root = uf.union(ru, rv)
            born[root] = old
    for root, members in uf.groups().items():
        dgm.append((f[born[root]], INF))
        pairs.append((born[root], None))
    return {'dgm': np.array(dgm, float).reshape(-1, 2), 'pairs': pairs}


def h0_superlevel(G, f: dict, nodes=None) -> dict:
    """H0 of the superlevel sets of f: pairs (birth = local max, death = merge value <= birth)."""
    res = h0_lower_star(G, {v: -x for v, x in f.items()}, nodes)
    res['dgm'] = -res['dgm']
    return res


# ----------------------------------------------------------------------------
# H0 of the terminals under the shortest-path metric
# ----------------------------------------------------------------------------

def terminal_h0(G, terminals, merged_groups=(), weight: str = "weight") -> dict:
    """
    Single-linkage H0 of `terminals` under the shortest-path metric of G.
    Each group of `merged_groups` (e.g. [sources]) is born as one point.

    Returns
    -------
    dict with
      deaths    : death times of the finite classes (= MST edges of the metric closure)
      total     : sum of the deaths (= weight of that MST)
      mst_edges : list of (label a, label b, length)
      n_labels  : number of points (terminals, a group counting as one)
      nn        : {label: distance to the nearest other label}
      label_of  : {terminal: label}
      region    : {node: label of its nearest terminal} (Voronoi regions of the terminals)
      dist      : {node: distance to its nearest terminal}
    """
    label_of, n_labels = {}, 0
    for group in merged_groups:
        for v in group:
            label_of[v] = n_labels
        n_labels += 1
    for t in terminals:
        if t not in label_of:
            label_of[t] = n_labels
            n_labels += 1
    nodes, index, M = graph_csr(G, weight)
    seeds = [index[v] for v in label_of]
    dist, _, src = dijkstra(M, directed=False, indices=seeds, min_only=True, return_predecessors=True)
    seed_label = {index[v]: lab for v, lab in label_of.items()}
    lab = np.array([seed_label.get(s, -1) for s in src])
    upper = triu(M, k=1).tocoo()
    r, c, w = upper.row, upper.col, upper.data
    keep = (lab[r] != lab[c]) & (lab[r] >= 0) & (lab[c] >= 0)
    cand = dist[r[keep]] + w[keep] + dist[c[keep]]
    la, lb = lab[r[keep]], lab[c[keep]]
    nn = np.full(n_labels, INF)
    np.minimum.at(nn, la, cand)
    np.minimum.at(nn, lb, cand)
    uf = UnionFind(range(n_labels))
    mst = []
    for k in np.argsort(cand, kind='stable'):
        if uf.union(int(la[k]), int(lb[k])) is not None:
            mst.append((int(la[k]), int(lb[k]), float(cand[k])))
            if len(mst) == n_labels - 1:
                break
    deaths = np.array([e[2] for e in mst])
    return {'deaths': deaths, 'total': float(deaths.sum()), 'mst_edges': mst, 'n_labels': n_labels,
            'nn': dict(enumerate(nn.tolist())), 'label_of': label_of,
            'region': {nodes[i]: int(lab[i]) for i in range(len(nodes))},
            'dist': {nodes[i]: float(dist[i]) for i in range(len(nodes))}}


def terminal_h0_diagram(res: dict) -> np.ndarray:
    """Diagram of terminal_h0: every point born at 0, plus one essential class per connected component."""
    n_ess = res['n_labels'] - len(res['deaths'])
    births = np.zeros(len(res['deaths']) + n_ess)
    deaths = np.concatenate([res['deaths'], np.full(n_ess, INF)])
    return np.column_stack([births, deaths])


# ----------------------------------------------------------------------------
# H1 of the planarised drawing
# ----------------------------------------------------------------------------

def _orient(ax, ay, bx, by, cx, cy):
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def edge_crossings(P: np.ndarray, edges: np.ndarray) -> list:
    """
    Proper crossings of the straight-line edges. `P` (n, 2) positions,
    `edges` (m, 2) node indices. Returns a list of (i, j, t_i, t_j, x, y): edge
    i and edge j cross at (x, y), at fraction t_i of edge i and t_j of edge j.
    Sweep on x so that only edges with overlapping x ranges are compared.
    """
    A, B = P[edges[:, 0]], P[edges[:, 1]]
    xmin, xmax = np.minimum(A[:, 0], B[:, 0]), np.maximum(A[:, 0], B[:, 0])
    ymin, ymax = np.minimum(A[:, 1], B[:, 1]), np.maximum(A[:, 1], B[:, 1])
    order = np.argsort(xmin)
    xs = xmin[order]
    out = []
    for pos, i in enumerate(order):
        end = np.searchsorted(xs, xmax[i], side='right')
        js = order[pos + 1:end]
        if len(js) == 0:
            continue
        js = js[(ymin[js] <= ymax[i]) & (ymax[js] >= ymin[i])]
        shared = ((edges[js, 0] == edges[i, 0]) | (edges[js, 0] == edges[i, 1]) |
                  (edges[js, 1] == edges[i, 0]) | (edges[js, 1] == edges[i, 1]))
        js = js[~shared]
        if len(js) == 0:
            continue
        a, b, c, d = A[i], B[i], A[js], B[js]
        o1 = _orient(a[0], a[1], b[0], b[1], c[:, 0], c[:, 1])
        o2 = _orient(a[0], a[1], b[0], b[1], d[:, 0], d[:, 1])
        o3 = _orient(c[:, 0], c[:, 1], d[:, 0], d[:, 1], a[0], a[1])
        o4 = _orient(c[:, 0], c[:, 1], d[:, 0], d[:, 1], b[0], b[1])
        hit = (o1 * o2 < 0) & (o3 * o4 < 0)
        for k in np.nonzero(hit)[0]:
            j = js[k]
            t_i = o3[k] / (o3[k] - o4[k])   # fraction along i where the line of j is crossed
            t_j = o1[k] / (o1[k] - o2[k])
            x, y = a + t_i * (b - a)
            out.append((int(i), int(j), float(t_i), float(t_j), float(x), float(y)))
    return out


def planarize(G, weight: str = "weight") -> dict:
    """
    Planarised drawing of G (positions in the node attribute 'pos').

    Returns dict with
      pos       : {vertex: (x, y)} (virtual vertices are ('x', k))
      sub       : list of sub-edges (a, b, parent edge index, value, real length)
      edges     : list of the parent edges (u, v) of G
      n_crossings
    The value of a sub-edge is the weight of its parent edge (filtration
    value), its real length is the fraction of that weight it covers.
    """
    nodes = list(G)
    idx = {v: i for i, v in enumerate(nodes)}
    P = np.array([G.nodes[v]['pos'] for v in nodes], float)
    edges = [(u, v) if str(u) <= str(v) else (v, u) for u, v in G.edges]
    E = np.array([(idx[u], idx[v]) for u, v in edges], int).reshape(-1, 2)
    weights = [float(G[u][v][weight]) for u, v in edges]
    crossings = edge_crossings(P, E)
    cuts = {i: [] for i in range(len(edges))}
    pos = {v: tuple(P[idx[v]]) for v in nodes}
    for k, (i, j, t_i, t_j, x, y) in enumerate(crossings):
        vx = ('x', k)
        pos[vx] = (x, y)
        cuts[i].append((t_i, vx))
        cuts[j].append((t_j, vx))
    sub = []
    for i, (u, v) in enumerate(edges):
        chain = [(0.0, u)] + sorted(cuts[i]) + [(1.0, v)]
        for (t0, a), (t1, b) in zip(chain, chain[1:]):
            sub.append((a, b, i, weights[i], (t1 - t0) * weights[i]))
    return {'pos': pos, 'sub': sub, 'edges': edges, 'n_crossings': len(crossings)}


def faces_of(pos: dict, sub: list) -> tuple:
    """
    Faces of a straight-line plane graph given by its sub-edges.
    Returns (face of each half-edge, number of faces); half-edge 2k goes from
    sub[k][0] to sub[k][1], half-edge 2k+1 the other way.
    """
    out = {}
    for k, (a, b, *_) in enumerate(sub):
        out.setdefault(a, []).append(2 * k)
        out.setdefault(b, []).append(2 * k + 1)
    dest = lambda h: sub[h // 2][1] if h % 2 == 0 else sub[h // 2][0]
    rot, rank = {}, {}
    for v, hs in out.items():
        x, y = pos[v]
        hs.sort(key=lambda h: math.atan2(pos[dest(h)][1] - y, pos[dest(h)][0] - x))
        rot[v] = hs
        for r, h in enumerate(hs):
            rank[h] = r
    face = [-1] * (2 * len(sub))
    n_faces = 0
    for h0 in range(2 * len(sub)):
        if face[h0] >= 0:
            continue
        h = h0
        while face[h] < 0:
            face[h] = n_faces
            b = dest(h)
            twin = h ^ 1
            hs = rot[b]
            h = hs[(rank[twin] - 1) % len(hs)]
        n_faces += 1
    return face, n_faces


def h1_planar(G, weight: str = "weight") -> dict:
    """
    H1 persistence of the planarised drawing of G (see module docstring).

    Returns dict with
      dgm              : (k, 2) array of (birth, death) of the finite classes
      edges            : parent edge (u, v) of each pair of dgm
      edge_persistence : {(u, v): largest persistence of a cycle closed by this edge}
      n_essential      : cycle-creating edges with no face to fill them (0 for a sphere)
      n_crossings, n_faces, euler (V - E + F, 2 for a connected plane graph)
    """
    pl = planarize(G, weight)
    sub, pos = pl['sub'], pl['pos']
    face, n_faces = faces_of(pos, sub)
    half_perim = np.zeros(n_faces)
    max_edge = np.zeros(n_faces)
    for k, (_, _, _, value, length) in enumerate(sub):
        for h in (2 * k, 2 * k + 1):
            half_perim[face[h]] += length / 2
            max_edge[face[h]] = max(max_edge[face[h]], value)
    face_value = np.maximum(half_perim, max_edge)
    # primal: edges that merge two components are negative (kill H0)
    order = sorted(range(len(sub)), key=lambda k: (sub[k][3], k))
    uf = UnionFind(pos)
    positive = [k for k in order if uf.union(sub[k][0], sub[k][1]) is None]
    # dual: faces merged in decreasing order of the positive edges (elder rule)
    duf = UnionFind(range(n_faces))
    oldest = {f: f for f in range(n_faces)}   # {root: face of largest value of the component}
    dgm, pair_edges, n_essential = [], [], 0
    for k in sorted(positive, key=lambda k: (-sub[k][3], -k)):
        fa, fb = duf.find(face[2 * k]), duf.find(face[2 * k + 1])
        if fa == fb:
            n_essential += 1
            continue
        oa, ob = oldest[fa], oldest[fb]
        young, old = (oa, ob) if face_value[oa] < face_value[ob] else (ob, oa)
        dgm.append((sub[k][3], float(face_value[young])))
        pair_edges.append(pl['edges'][sub[k][2]])
        oldest[duf.union(fa, fb)] = old
    edge_pers = {}
    for (b, d), e in zip(dgm, pair_edges):
        edge_pers[e] = max(edge_pers.get(e, 0.0), d - b)
    n_vertices = len(pos)
    return {'dgm': np.array(dgm, float).reshape(-1, 2), 'edges': pair_edges,
            'edge_persistence': edge_pers, 'n_essential': n_essential,
            'n_crossings': pl['n_crossings'], 'n_faces': n_faces,
            'euler': n_vertices - len(sub) + n_faces}
