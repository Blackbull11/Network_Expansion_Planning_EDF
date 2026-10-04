"""
Voronoi cells with a variable distance: each source s has a free weight w_s
(edge costs from s are divided by w_s, as in voronoi.py), adjusted so that the
Voronoi + Steiner solution satisfies the constraints.

At each iteration:
  1. Voronoi cells with the current weights, Steiner tree in each cell;
  2. the solution is checked by verify.check_solution;
  3. each source gets a load ratio r_s = max over its constraints of value / limit:
       |P_source| / capacity, |P_source| / (n_departs_max * p_max_depart),
       n_départs / n_departs_max, L_départ / l_max_depart, |P_départ| / p_max_depart;
  4. the weight of each overloaded source (r_s > 1) is divided by r_s ** step:
     its cell shrinks and its boundary nodes go to the neighbouring cells.
     The step decreases geometrically so that the weights settle.
It stops as soon as the solution is feasible (no source is overloaded any more,
the weights would not change). Cells stay connected whatever the weights
(each node is reached by a path inside its cell), so Steiner always applies.

The best solution met is kept (verify.solution_key: feasible first, then the
smallest J). compute_voronoi_variable runs this from the 3 weightings of
voronoi.py (none, capacity, sqrt) and keeps the best result.
"""
from voronoi import WEIGHTINGS, compute_voronoi_cells, source_weight, voronoi_cells_SB
from steiner import steiner_trees_by_cells
from verify import check_solution, solution_key
from nx_graph import GraphData
from src.outils import init_stream_logger

logger = init_stream_logger()

MAX_SHRINK = 0.5   # a weight is divided by at most 2 per iteration


def load_ratio(src: dict) -> float:
    """r_s = max value / limit of the constraints of a source (see module docstring)."""
    ratios = [abs(src['power']) / src['capacity'],
              abs(src['power']) / (src['n_departs_max'] * src['p_max_depart']),
              src['n_departs'] / src['n_departs_max']]
    for d in src['departs']:
        ratios += [d['length'] / src['l_max_depart'], abs(d['power']) / src['p_max_depart']]
    return max(ratios)


def voronoi_variable_distance(graph: GraphData, init="capacity", n_iter: int = 60,
                              step: float = 0.5, decay: float = 0.95) -> dict:
    """
    Adjust the source weights from the weighting `init` (one of WEIGHTINGS).

    Parameters
    ----------
    n_iter : maximum number of iterations
    step, decay : at iteration k, w_s is divided by r_s ** (step * decay ** k)

    Returns
    -------
    best : dict
        {weights: {source: w}, cells, trees: {source: edges}, cert: certificate
        of the trees (verify.check_solution), init, iteration: iteration of the
        best solution, n_iter: iterations done}
    """
    G = graph.G
    sources = [n for n, label in G.nodes('label') if label == 'Source']
    w = {s: source_weight(G, s, init) for s in sources}
    best = None
    for k in range(n_iter):
        _, cells, _ = voronoi_cells_SB(G, sources, weighting=w)
        trees = steiner_trees_by_cells(G, cells, log=False)
        cert = check_solution(graph, trees)
        if best is None or solution_key(cert) < solution_key(best['cert']):
            best = {'weights': dict(w), 'cells': cells, 'trees': trees, 'cert': cert,
                    'init': init, 'iteration': k}
        if cert['feasible']:
            break
        ratios = {s: load_ratio(src) for s, src in cert['sources'].items()}
        overloaded = {s: r for s, r in ratios.items() if r > 1}
        if not overloaded:   # infeasible for structural reasons, weights cannot help
            break
        for s, r in overloaded.items():
            w[s] *= max(r ** -(step * decay ** k), MAX_SHRINK)
    best['n_iter'] = k + 1
    return best


def compute_voronoi_variable(graph: GraphData, inits=WEIGHTINGS, **kwargs) -> dict:
    """
    voronoi_variable_distance from each weighting of `inits`, keep the best
    solution, and store it in `graph` like compute_voronoi_cells +
    compute_steiner_trees (cells, split edges for plotting, steiner_edges),
    plus graph.voronoi_weights = {source: w}. Returns the best result.
    """
    runs = [voronoi_variable_distance(graph, init, **kwargs) for init in inits]
    for r in runs:
        logger.info("Distance variable depuis '%s' : %s après %d itérations (meilleure : itération %d), J = %.1f",
                    r['init'], "faisable" if r['cert']['feasible'] else "infaisable",
                    r['n_iter'], r['iteration'], r['cert']['J'])
    best = min(runs, key=lambda r: solution_key(r['cert']))
    compute_voronoi_cells(graph, weighting=best['weights'])
    graph.steiner_edges = best['trees']
    graph.voronoi_weights = best['weights']
    return best
