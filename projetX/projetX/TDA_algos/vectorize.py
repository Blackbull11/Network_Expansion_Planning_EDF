"""
Vectorisations of persistence diagrams and distances between diagrams.

  summary            : scalar statistics (number of bars, total persistence, entropy...)
  persistence_image  : persistence image (Adams et al. 2017) on a fixed grid
  betti_curve        : number of bars alive at each value of a grid
  wasserstein        : exact p-Wasserstein distance (Hungarian algorithm on the
                       diagonal-augmented cost matrix, top `max_points` bars)
  bottleneck         : exact bottleneck distance (binary search on the cost
                       values + Hopcroft-Karp perfect matching)
Essential classes (death = inf) are compared by their births, sorted; a
different number of essential classes gives an infinite distance.
"""
import math

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

INF = math.inf


def split(dgm: np.ndarray) -> tuple:
    """(finite pairs, births of the essential classes)."""
    dgm = np.asarray(dgm, float).reshape(-1, 2)
    fin = np.isfinite(dgm[:, 1])
    return dgm[fin], np.sort(dgm[~fin, 0])


def persistence(dgm: np.ndarray) -> np.ndarray:
    """|death - birth| of the finite pairs (superlevel diagrams have death <= birth)."""
    fin, _ = split(dgm)
    return np.abs(fin[:, 1] - fin[:, 0])


def summary(dgm: np.ndarray, prefix: str = "", eps: float = 1e-9) -> dict:
    """Scalar statistics of a diagram, keys prefixed by `prefix`."""
    fin, ess = split(dgm)
    p = np.abs(fin[:, 1] - fin[:, 0])
    p = p[p > eps]
    tot = p.sum()
    q = p / tot if tot > 0 else p
    out = {'n_bars': len(p), 'n_essential': len(ess), 'total_pers': tot,
           'max_pers': p.max() if len(p) else 0.0, 'mean_pers': p.mean() if len(p) else 0.0,
           'std_pers': p.std() if len(p) else 0.0,
           'q50_pers': np.quantile(p, 0.5) if len(p) else 0.0,
           'q90_pers': np.quantile(p, 0.9) if len(p) else 0.0,
           'entropy': float(-(q * np.log(q)).sum()) if len(p) else 0.0,
           'gap_ratio': largest_gap(p)[1]}
    return {prefix + k: float(v) for k, v in out.items()}


def largest_gap(p: np.ndarray) -> tuple:
    """
    (threshold, ratio): the persistence values sorted decreasingly, the largest
    ratio p[i] / p[i+1] between consecutive values. `threshold` is the
    geometric mean of the two values of that gap: bars above are the
    "significant" ones. (inf, 1) if fewer than 2 values.
    """
    p = np.sort(np.asarray(p, float)[np.asarray(p) > 0])[::-1]
    if len(p) < 2:
        return (0.0 if len(p) else INF), 1.0
    ratios = p[:-1] / p[1:]
    i = int(np.argmax(ratios))
    return float(math.sqrt(p[i] * p[i + 1])), float(ratios[i])


def persistence_image(dgm: np.ndarray, grid: int = 10, birth_range=None, pers_range=None,
                      sigma: float = None) -> np.ndarray:
    """
    Persistence image: Gaussian kernels at (birth, persistence) weighted by the
    persistence, integrated on a grid x grid lattice. Ranges default to the
    diagram's own; fix them to compare diagrams. Returns a flat vector.
    """
    fin, _ = split(dgm)
    b = np.minimum(fin[:, 0], fin[:, 1])
    p = np.abs(fin[:, 1] - fin[:, 0])
    if birth_range is None:
        birth_range = (b.min(), b.max()) if len(b) else (0.0, 1.0)
    if pers_range is None:
        pers_range = (0.0, p.max()) if len(p) else (0.0, 1.0)
    bx = np.linspace(*birth_range, grid)
    py = np.linspace(*pers_range, grid)
    if sigma is None:
        sigma = max((pers_range[1] - pers_range[0]) / grid, 1e-9)
    img = np.zeros((grid, grid))
    for bi, pi in zip(b, p):
        gx = np.exp(-((bx - bi) ** 2) / (2 * sigma ** 2))
        gy = np.exp(-((py - pi) ** 2) / (2 * sigma ** 2))
        img += pi * np.outer(gy, gx)
    return img.ravel()


def betti_curve(dgm: np.ndarray, ts: np.ndarray) -> np.ndarray:
    """Number of bars alive at each t of `ts` (birth <= t < death, sublevel convention)."""
    d = np.asarray(dgm, float).reshape(-1, 2)
    lo, hi = np.minimum(d[:, 0], d[:, 1]), np.maximum(d[:, 0], d[:, 1])
    return np.array([int(((lo <= t) & (t < hi)).sum()) for t in ts])


def _essential_distance(e1, e2, p):
    if len(e1) != len(e2):
        return INF
    if not len(e1):
        return 0.0
    diff = np.abs(e1 - e2)
    return float(diff.max()) if p == INF else float((diff ** p).sum())


def _top(fin, max_points):
    if max_points is None or len(fin) <= max_points:
        return fin
    p = np.abs(fin[:, 1] - fin[:, 0])
    return fin[np.argsort(-p)[:max_points]]


def wasserstein(d1, d2, p: float = 1.0, max_points: int = 500) -> float:
    """
    p-Wasserstein distance (L-inf ground metric). Only the `max_points` most
    persistent bars of each diagram are matched (the rest go to the diagonal).
    """
    f1, e1 = split(d1)
    f2, e2 = split(d2)
    ess = _essential_distance(e1, e2, p)
    if ess == INF:
        return INF
    # bars beyond max_points are matched to the diagonal
    rest = 0.0
    for f in (f1, f2):
        if max_points is not None and len(f) > max_points:
            pers = np.sort(np.abs(f[:, 1] - f[:, 0]))[::-1][max_points:]
            rest += float(((pers / 2) ** p).sum())
    f1, f2 = _top(f1, max_points), _top(f2, max_points)
    n, m = len(f1), len(f2)
    if n + m == 0:
        return (ess + rest) ** (1 / p)
    C = np.zeros((n + m, n + m))
    if n and m:
        C[:n, :m] = np.max(np.abs(f1[:, None, :] - f2[None, :, :]), axis=2) ** p
    diag1 = (np.abs(f1[:, 1] - f1[:, 0]) / 2) ** p
    diag2 = (np.abs(f2[:, 1] - f2[:, 0]) / 2) ** p
    C[:n, m:] = INF
    C[n:, :m] = INF
    if n:
        C[np.arange(n), m + np.arange(n)] = diag1
    if m:
        C[n + np.arange(m), np.arange(m)] = diag2
    C[n:, m:] = 0.0
    C[np.isinf(C)] = 1e18
    r, c = linear_sum_assignment(C)
    return float((C[r, c].sum() + ess + rest) ** (1 / p))


def bottleneck(d1, d2) -> float:
    """Bottleneck distance (L-inf ground metric)."""
    f1, e1 = split(d1)
    f2, e2 = split(d2)
    ess = _essential_distance(e1, e2, INF)
    n, m = len(f1), len(f2)
    if n + m == 0:
        return ess
    # left: points of d1 then diagonal copies of d2; right: points of d2 then diagonal copies of d1
    C = np.full((n + m, m + n), INF)
    if n and m:
        C[:n, :m] = np.max(np.abs(f1[:, None, :] - f2[None, :, :]), axis=2)
    if n:
        C[np.arange(n), m + np.arange(n)] = np.abs(f1[:, 1] - f1[:, 0]) / 2
    if m:
        C[n + np.arange(m), np.arange(m)] = np.abs(f2[:, 1] - f2[:, 0]) / 2
    C[n:, m:] = 0.0
    values = np.unique(C[np.isfinite(C)])
    lo, hi = 0, len(values) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        A = csr_matrix(C <= values[mid])
        if (maximum_bipartite_matching(A, perm_type='column') >= 0).all():
            hi = mid
        else:
            lo = mid + 1
    return float(max(values[lo], ess))
