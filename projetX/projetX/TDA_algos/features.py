"""
Instance descriptors and configuration selection.

simple_features : size, slack of the capacities, limits (no topology)
tda_features    : summaries (vectorize.summary) of
                    h0_term   H0 of the loads with the sources merged (single linkage)
                    h1        H1 of the planarised graph (holes / alternative routes)
                    h0_dist   superlevel H0 of the distance to the nearest source (corridors)
                    tomato    ToMATo diagram of the power density (load clusters)
                  plus the lower bound, the spectral gap at k = number of sources
                  and (optional) the mean effective resistance / distance to the
                  nearest source (redundancy of the routes).
                  Lengths are divided by l_max_depart, densities by their max,
                  so that instances of different sizes are comparable.
knn_select      : leave-one-group-out k-nearest-neighbour prediction of the
                  best configuration of each instance from a feature table.
"""
import math
import time

import numpy as np

from .bounds import lower_bound
from .clustering import cluster_loads
from .persistence import h0_superlevel, h1_planar, shortest_distances, terminal_h0_diagram
from .vectorize import summary

LOAD_LABELS = ('Conso', 'Prod')


def _limits(graph):
    from verify import source_limits
    sources = [n for n, lab in graph.G.nodes('label') if lab == 'Source']
    lims = [source_limits(graph, s) for s in sources]
    return sources, lims


def simple_features(graph) -> dict:
    G = graph.G
    sources, lims = _limits(graph)
    loads = [n for n, lab in G.nodes('label') if lab in LOAD_LABELS]
    d = np.array([G.nodes[l]['conso_ete'] for l in loads]) if loads else np.zeros(1)
    cap = sum(min(l['capacity'], l['n_departs_max'] * l['p_max_depart']) for l in lims)
    l_max = min(l['l_max_depart'] for l in lims)
    p_max = min(l['p_max_depart'] for l in lims)
    _, index, D = shortest_distances(G, sources)
    nearest = D[:, [index[l] for l in loads]].min(axis=0) if loads else np.zeros(1)
    return {'n_nodes': G.number_of_nodes(), 'n_edges': G.number_of_edges(), 'n_loads': len(loads),
            'n_sources': len(sources), 'loads_per_source': len(loads) / max(len(sources), 1),
            'frac_prod': float((d < 0).mean()), 'mean_degree': 2 * G.number_of_edges() / G.number_of_nodes(),
            'net_over_cap': abs(d.sum()) / cap if cap > 0 else math.inf,
            'gross_over_cap': d[d > 0].sum() / cap if cap > 0 else math.inf,
            'gross_over_pmax': d[d > 0].sum() / p_max if p_max > 0 else math.inf,
            'mean_dist_over_lmax': float(nearest.mean() / l_max) if math.isfinite(l_max) else 0.0,
            'max_dist_over_lmax': float(nearest.max() / l_max) if math.isfinite(l_max) else 0.0}


def tda_features(graph, spectral: bool = True, resistance: bool = False, timings: dict = None) -> tuple:
    """
    (features, diagrams): features dict (see module docstring) and the raw
    diagrams {'h0_term', 'h1', 'h0_dist', 'tomato'}. `timings` (dict), if
    given, receives the computing time of each component.
    """
    G = graph.G
    sources, lims = _limits(graph)
    l_max = min(l['l_max_depart'] for l in lims)
    scale = l_max if math.isfinite(l_max) and l_max > 0 else 1.0
    timings = {} if timings is None else timings
    feats, dgms = {}, {}

    t = time.time()
    lb = lower_bound(graph)
    dgms['h0_term'] = terminal_h0_diagram(lb['h0'])
    feats.update(summary(dgms['h0_term'] / scale, 'h0term_'))
    feats['lb_over_lmax'] = lb['lb'] / scale
    timings['h0_term'] = time.time() - t

    t = time.time()
    h1 = h1_planar(G)
    dgms['h1'] = h1['dgm']
    feats.update(summary(h1['dgm'] / scale, 'h1_'))
    feats['h1_crossings'] = h1['n_crossings']
    timings['h1'] = time.time() - t

    t = time.time()
    _, index, D = shortest_distances(G, sources)
    f = {v: float(D[:, index[v]].min()) for v in G}
    dgms['h0_dist'] = h0_superlevel(G, f)['dgm']
    feats.update(summary(dgms['h0_dist'] / scale, 'h0dist_'))
    timings['h0_dist'] = time.time() - t

    t = time.time()
    cl = cluster_loads(G)
    top = max(cl['rho'].values()) or 1.0
    dgms['tomato'] = cl['full_dgm']
    feats.update(summary(cl['full_dgm'] / top, 'tomato_'))
    feats['tomato_n_clusters'] = len([p for p in cl['clusters'] if p != -1])
    feats['tomato_clusters_per_source'] = feats['tomato_n_clusters'] / max(len(sources), 1)
    timings['tomato'] = time.time() - t

    if spectral:
        from .spectral import spectral_embedding
        t = time.time()
        k = max(2, len(sources))
        ev = spectral_embedding(G, k)['eigenvalues']
        feats['lambda2'] = float(ev[1]) if len(ev) > 1 else 0.0
        feats['eigengap_k'] = float(ev[-1] - ev[-2]) if len(ev) > 2 else 0.0
        timings['spectral'] = time.time() - t
    if resistance:
        from .spectral import effective_resistance
        t = time.time()
        loads = [n for n, lab in G.nodes('label') if lab in LOAD_LABELS]
        R = effective_resistance(G, sources, loads)
        ratios = []
        for l in loads:
            j = int(np.argmin(D[:, index[l]]))
            d = D[j, index[l]]
            if d > 0:
                ratios.append(R[sources[j]][l] / d)
        feats['resistance_ratio_mean'] = float(np.mean(ratios)) if ratios else 1.0
        timings['resistance'] = time.time() - t
    return feats, dgms


def knn_select(X: np.ndarray, y: list, groups: list, k: int = 3) -> dict:
    """
    Leave-one-group-out k-NN (features standardised on the training part).
    Returns {accuracy, predictions, majority_accuracy (always predicting the
    most frequent label of the training part)}.
    """
    X = np.asarray(X, float)
    X = np.where(np.isfinite(X), X, np.nan)
    col_mean = np.nanmean(X, axis=0)
    X = np.where(np.isnan(X), col_mean, X)
    y, groups = np.asarray(y), np.asarray(groups)
    pred, base = [], []
    for i in range(len(y)):
        train = groups != groups[i]
        if not train.any():
            pred.append(None)
            base.append(None)
            continue
        mu, sd = X[train].mean(0), X[train].std(0)
        sd[sd == 0] = 1
        Z = (X - mu) / sd
        d = np.linalg.norm(Z[train] - Z[i], axis=1)
        nn = y[train][np.argsort(d)[:k]]
        vals, counts = np.unique(nn, return_counts=True)
        pred.append(vals[np.argmax(counts)])
        vals, counts = np.unique(y[train], return_counts=True)
        base.append(vals[np.argmax(counts)])
    ok = [p == t for p, t in zip(pred, y) if p is not None]
    ok_base = [p == t for p, t in zip(base, y) if p is not None]
    return {'accuracy': float(np.mean(ok)) if ok else math.nan,
            'majority_accuracy': float(np.mean(ok_base)) if ok_base else math.nan,
            'predictions': pred}
