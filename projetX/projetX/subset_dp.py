"""
Borne inférieure par programmation dynamique sur les sous-ensembles de charges.

  1. Dreyfus-Wagner sur le graphe privé des sources : dp[S][v] = plus court arbre
     reliant les charges S et le sommet v.
  2. Pour chaque source r : cost_r[S] = plus court départ de r desservant
     exactement S (une arête r-h puis un arbre sans autre source). S est permis
     si |P(S)| <= p_max_depart et cost_r[S] <= l_max_depart.
  3. F_r[M] = meilleure façon pour r de desservir M avec au plus
     min(n_departs_max, degré) départs, |P(M)| <= capacité de r.
  4. B[M] = meilleure affectation de M aux sources.

B[toutes les charges] ne relâche que la disjonction des départs entre eux :
c'est une borne inférieure de J, et l'optimum quand les départs trouvés sont
disjoints. Si B est infini, l'instance est infaisable.

Coût : environ 3^k opérations pour k charges, 2^k x |V| réels en mémoire. Au-delà
d'une vingtaine de charges, le calcul n'est plus raisonnable.
"""
import time

import numpy as np

from verify import source_limits

LOAD_LABELS = ('Conso', 'Prod')
BATCH = 1 << 16   # sous-ensembles traités à la fois, pour borner la mémoire


def all_submasks(mask: int) -> np.ndarray:
    """Tous les sous-ensembles de `mask`, vide et `mask` compris."""
    a = np.array([0], dtype=np.int32)
    while mask:
        bit = mask & -mask
        a = np.concatenate([a, a | bit])
        mask ^= bit
    return a


def low_submasks(mask: int) -> np.ndarray:
    """Sous-ensembles stricts de `mask` contenant son plus petit élément."""
    low = mask & -mask
    a = all_submasks(mask ^ low) | low
    return a[a != mask]


def shortest_distances(H) -> tuple:
    """(nœuds, indice, matrice des plus courts chemins) de H."""
    nodes = list(H)
    index = {v: i for i, v in enumerate(nodes)}
    D = np.full((len(nodes), len(nodes)), np.inf)
    np.fill_diagonal(D, 0)
    for u, v, w in H.edges.data('weight'):
        D[index[u], index[v]] = D[index[v], index[u]] = w
    for m in range(len(nodes)):
        D = np.minimum(D, D[:, m, None] + D[m, None, :])
    return nodes, index, D


def merge(dp: np.ndarray, mask: int) -> np.ndarray:
    """min sur S (S contient le plus petit élément de mask) de dp[S] + dp[mask \\ S], sommet par sommet."""
    subs = low_submasks(mask)
    best = np.full(dp.shape[1], np.inf)
    for start in range(0, len(subs), BATCH):
        part = subs[start:start + BATCH]
        best = np.minimum(best, (dp[part] + dp[mask ^ part]).min(axis=0))
    return best


def source_table(c: np.ndarray, allowed: np.ndarray, n_max: int, power: np.ndarray, capacity: float) -> np.ndarray:
    """F[M] : coût minimal pour une source de desservir M (voir le docstring du module)."""
    F = np.where(allowed, c, np.inf)
    F[0] = 0.0
    for _ in range(n_max - 1):
        previous = F.copy()
        for mask in range(1, len(F)):
            if mask & (mask - 1) == 0:
                continue
            subs = low_submasks(mask)
            subs = subs[allowed[subs]]
            if len(subs):
                F[mask] = min(F[mask], float((c[subs] + previous[mask ^ subs]).min()))
        if np.array_equal(F, previous):
            break
    F = np.where(np.abs(power) <= capacity + 1e-9, F, np.inf)
    F[0] = 0.0
    return F


def lower_bound(graph) -> dict:
    """{bound (inf si infaisable), n_loads, time} : borne inférieure de J pour `graph`."""
    t0 = time.perf_counter()
    G = graph.G
    sources = [n for n, lab in G.nodes('label') if lab == 'Source']
    loads = [n for n, lab in G.nodes('label') if lab in LOAD_LABELS]
    H = G.copy()
    H.remove_nodes_from(sources)
    nodes, index, D = shortest_distances(H)
    n, k = len(nodes), len(loads)
    full = (1 << k) - 1
    dp = np.full((full + 1, n), np.inf)
    for i, t in enumerate(loads):
        dp[1 << i] = D[index[t]]
    for mask in range(1, full + 1):
        if mask & (mask - 1):
            dp[mask] = (merge(dp, mask)[None, :] + D).min(axis=1)
    power = np.zeros(full + 1)
    for i, t in enumerate(loads):
        bit = 1 << i
        power[bit:2 * bit] = power[:bit] + G.nodes[t]['conso_ete']
    best = np.full(full + 1, np.inf)
    best[0] = 0.0
    for r in sources:
        neighbours = [h for h in G[r] if h in index]
        lim = source_limits(graph, r)
        c = (dp[:, [index[h] for h in neighbours]] + np.array([G[r][h]['weight'] for h in neighbours])).min(axis=1)
        c[0] = np.inf
        allowed = (np.abs(power) <= lim['p_max_depart'] + 1e-9) & (c <= lim['l_max_depart'] + 1e-9)
        allowed[0] = False
        F = source_table(c, allowed, int(min(lim['n_departs_max'], len(neighbours))), power, lim['capacity'])
        new = np.full(full + 1, np.inf)
        for mask in range(full + 1):
            subs = all_submasks(mask)
            new[mask] = float((F[subs] + best[mask ^ subs]).min())
        best = new
    return dict(bound=float(best[full]), n_loads=k, time=time.perf_counter() - t0)
