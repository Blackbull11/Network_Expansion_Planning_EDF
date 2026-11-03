"""
Run of the Voronoi + Steiner heuristic of main.py on the instances of
verified_instances/ (format with an "electrotechnique" block, see nx_graph.py).
The size, density, season and net load of each instance come from the table
verified_instances/ratios_charge_nette.csv.

On top of the evaluation of main.py, for each instance and distance:
  time_s          : computation time of the distance (Voronoi + Steiner, all the
                    iterations and starts for 'variable')
  max_flow        : largest |power flowing on an edge| of the trees (MW)
  lines_over      : number of edges carrying more than puissance_max_ligne
  max_pd_depart   : largest PD of a départ (MW.km)
  departs_over_pl : number of départs with PD > PL_max_noeud
puissance_max_ligne is checked by verify.py on the départ heads only (p_max_depart),
and PL_max_noeud is not checked: these are indicators, not constraints.

Results in output/: recap_<stamp>.md/.csv, details_<stamp>.txt and charts as
main.py, plus verified_<stamp>.csv (one row per instance and distance, with the
indicators above, written as the run goes) and solutions_<stamp>.json (trees).

Usage:
  python run_verified.py                          all the instances
  python run_verified.py Petit_N60_S2_0.20 --label test
                                                  instances whose name contains a pattern
  python run_verified.py --figures variable       figures of the variable distance only
                                                  (all: the 4 distances, none: no figure, the default)
"""
import argparse
import csv
import json
import time
from datetime import datetime
from pathlib import Path

import networkx as nx

from nx_graph import GraphData
from voronoi import compute_voronoi_cells
from voronoi_variable_distance import compute_voronoi_variable
from steiner import compute_steiner_trees, show_steiner_trees
from verify import TERMINAL_LABELS, TOL, check_solution, necessary_conditions, source_param
from output import OUTPUT_DIR, METHODS, WEIGHTING_NAMES, format_evaluation, write_summary, write_details

VERIFIED_DIR = Path(__file__).resolve().parent / "verified_instances"
RATIOS_FILE = VERIFIED_DIR / "ratios_charge_nette.csv"
FIGURES_DIR = OUTPUT_DIR / "figures"
SIZE_ORDER = ('Petit', 'Moyen', 'Grand', 'TrèsGrand')
SEASON_ORDER = ('Ete', 'Hiver')


def sort_key(row: dict) -> tuple:
    return (SIZE_ORDER.index(row['taille']), int(row['N']), int(row['S']), float(row['densite']),
            SEASON_ORDER.index(row['saison']))


def list_instances(folder: Path = VERIFIED_DIR) -> list:
    """
    [(name, path, row of ratios_charge_nette)] of the instances of `folder` that
    have a JSON file, in the order Petit -> TrèsGrand, then N, S, densité, saison.
    """
    with open(RATIOS_FILE, encoding="utf-8") as f:
        rows = sorted(csv.DictReader(f), key=sort_key)
    return [(row['nom'], folder / f"{row['nom']}.json", row) for row in rows
            if (folder / f"{row['nom']}.json").exists()]


def solve(graph: GraphData, method: str) -> dict:
    """Voronoi + Steiner for one of METHODS; for 'variable', the best run of compute_voronoi_variable."""
    if method == 'variable':
        return compute_voronoi_variable(graph)
    compute_voronoi_cells(graph, weighting=method)
    compute_steiner_trees(graph)
    return {}


def edge_flows(G: nx.Graph, source, tree_edges: list) -> dict:
    """{n: power flowing on the edge (parent of n, n)} of the tree rooted at `source` (as in verify.py)."""
    T = nx.Graph(tree_edges)
    if source not in T:
        return {}
    parent = nx.dfs_predecessors(T, source)
    flow = {}
    for n in nx.dfs_postorder_nodes(T, source):
        if n != source:
            flow[n] = G.nodes[n]['conso_ete'] + sum(flow[c] for c in T.neighbors(n) if c != parent[n])
    return flow


def indicators(graph: GraphData, ev: dict) -> dict:
    """max_flow, lines_over, max_pd_depart, departs_over_pl (see module docstring) of the valid trees."""
    max_flow, lines_over, max_pd, over_pl = 0.0, 0, 0.0, 0
    for s, src in ev['sources'].items():
        p_line, pl_max = src['p_max_depart'], source_param(graph, s, 'pl_max')
        flows = [abs(f) for f in edge_flows(graph.G, s, graph.steiner_edges[s]).values()]
        max_flow = max([max_flow] + flows)
        lines_over += sum(f > p_line + TOL for f in flows)
        pds = [d['pd'] for d in src['departs']]
        max_pd = max([max_pd] + pds)
        over_pl += sum(pd_ > pl_max + TOL for pd_ in pds)
    return {'max_flow': max_flow, 'lines_over': lines_over, 'max_pd_depart': max_pd,
            'departs_over_pl': over_pl}


CSV_COLUMNS = ['nom', 'taille', 'N', 'S', 'densite', 'saison', 'n_charges', 'n_neg',
               'conso_MW', 'prod_MW', 'charge_nette_MW', 'ratio_sources', 'prouvee_infaisable',
               'distance', 'J', 'PD', 'faisable', 'erreurs', 'viol_L', 'viol_P', 'viol_source', 'viol_n',
               'n_departs', 'time_s', 'iterations', 'init', 'max_flow', 'lines_over',
               'max_pd_depart', 'departs_over_pl']


def run(instances: list, label: str = None, figures: str = "none") -> list:
    """
    Voronoi + Steiner + evaluation of each instance (list of list_instances) and each distance.
    `figures`: 'all', 'variable' or 'none' (see module docstring).
    """
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S") + (f"_{label}" if label else "")
    csv_path = OUTPUT_DIR / f"verified_{stamp}.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerow(CSV_COLUMNS)
    results, details, solutions = [], [], {}
    for name, path, info in instances:
        G = GraphData(folderpath=str(path.parent), filename=path.name)
        print(G)
        loads = G.nodes[G.nodes['type'].isin(TERMINAL_LABELS)]
        proofs = necessary_conditions(G)
        if proofs:
            lines = [f"=== {name} : instance prouvée infaisable ==="] + [f"  {p}" for p in proofs]
            print("\n".join(lines))
            details += lines + [""]
        for weighting in METHODS:
            t0 = time.perf_counter()
            best = solve(G, weighting)
            elapsed = time.perf_counter() - t0
            ev = check_solution(G, G.steiner_edges)
            lines = format_evaluation(name, weighting, ev,
                                      G.voronoi_weights if weighting == "variable" else None)
            lines.insert(1, f"temps de calcul : {elapsed:.2f} s")
            print("\n".join(lines))
            details += lines + [""]
            results.append({'name': name, 'type': info['taille'],
                            'n_loads': len(loads), 'n_neg': int((loads['demande'] < 0).sum()),
                            'n_sources': len(ev['sources']), 'weighting': weighting, 'ev': ev})
            solutions.setdefault(name, {})[weighting] = {str(s): e for s, e in G.steiner_edges.items()}
            viol = ev['violations']
            row = {'nom': name, 'taille': info['taille'], 'N': info['N'],
                   'S': info['S'], 'densite': info['densite'], 'saison': info['saison'],
                   'n_charges': len(loads), 'n_neg': results[-1]['n_neg'],
                   'conso_MW': info['conso_MW'], 'prod_MW': info['prod_MW'],
                   'charge_nette_MW': info['charge_nette_MW'], 'ratio_sources': info['ratio_sources'],
                   'prouvee_infaisable': int(bool(proofs)), 'distance': weighting,
                   'J': round(ev['J'], 3), 'PD': round(ev['PD'], 3), 'faisable': int(ev['feasible']),
                   'erreurs': len(ev['errors']), 'viol_L': viol['L'][0], 'viol_P': viol['P'][0],
                   'viol_source': viol['source'][0], 'viol_n': viol['n'][0],
                   'n_departs': sum(src['n_departs'] for src in ev['sources'].values()),
                   'time_s': round(elapsed, 3), 'iterations': best.get('n_iter', ''),
                   'init': best.get('init', ''),
                   **{k: round(v, 3) for k, v in indicators(G, ev).items()}}
            with open(csv_path, "a", newline="", encoding="utf-8") as f:
                csv.writer(f, delimiter=";").writerow([row[c] for c in CSV_COLUMNS])
            if figures == "all" or (figures == "variable" and weighting == "variable"):
                show_steiner_trees(G, title=f"{name} - distance {WEIGHTING_NAMES[weighting]}"
                                            f" - J = {ev['J']:.1f}",
                                   filepath=FIGURES_DIR / f"{name}_{weighting}.png", show=False)
    (OUTPUT_DIR / f"solutions_{stamp}.json").write_text(json.dumps(solutions), encoding="utf-8")
    summary = write_summary(results, stamp)
    write_details(details, stamp)
    print(f"\nRécapitulatif écrit dans {summary}\nIndicateurs par instance dans {csv_path}")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Voronoi + Steiner sur les instances de verified_instances")
    parser.add_argument("motifs", nargs="*", help="ne garder que les instances dont le nom contient un de ces motifs")
    parser.add_argument("--label", help="suffixe du nom des rapports (recap_<date>_<label>.md ...)")
    parser.add_argument("--figures", choices=("all", "variable", "none"), default="none",
                        help="figures des 4 distances, de la distance variable seulement, ou aucune")
    args = parser.parse_args()
    instances = [x for x in list_instances() if not args.motifs or any(m in x[0] for m in args.motifs)]
    run(instances, label=args.label, figures=args.figures)
