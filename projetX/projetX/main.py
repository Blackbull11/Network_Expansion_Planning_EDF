"""
Run on the instances of data/instances (see data/generate_instances.py): for
each graph and each Voronoi distance (non pondérée, pondérée par la capacité,
pondérée par sa racine carrée, variable = poids ajustés aux contraintes, see
voronoi_variable_distance.py), Voronoi cells of the sources, Steiner tree in
each cell, then check of the solution by verify.check_solution.

Départ : subtree of the Steiner tree hanging from one neighbour of the source.
  L_départ = total length of the subtree (incl. the edge from the source)
  P_départ = sum of the 'demande' of its nodes (productions < 0 offset consumptions)
Constraints:
  L_départ <= l_max_depart, |P_départ| <= p_max_depart   (attributes of the source)
  |P_source| <= capacity of the source ('demande'), nb départs <= n_departs_max

Objective, to minimise (as in data/exemple_ete_4_PS/model.lp): J = L_tot, total
length of the trees (km). The constraints are hard: a solution violating one of
them, or whose trees are not a valid radial network covering all the loads
(see verify.py), is infeasible.
Also computed (not in J): PD = sum over the edges of the trees of
|power flowing on the edge| * length of the edge (MW.km).

Results: printed, and written in output/ (summary table recap_*.md / .csv and
per-départ details details_*.txt, charts in output/charts/). The Voronoi cells and Steiner trees are
saved as images in output/figures/<graph>_<distance>.png.

Usage:
  python main.py                 all the instances
  python main.py g05 montagne    only the instances whose name contains g05 or montagne
  python main.py g05 --show      also open the figures (one window per distance)
  python main.py g21_ g22_ --label grandes   reports named recap_<date>_grandes.md ...
  python main.py --petites       the 20 small instances g01 to g20
  python main.py --grandes       the 20 large instances g21 to g40
"""
import argparse
from datetime import datetime
from pathlib import Path

from nx_graph import GraphData
from voronoi import compute_voronoi_cells
from voronoi_variable_distance import compute_voronoi_variable
from steiner import compute_steiner_trees, show_steiner_trees
from verify import check_solution, necessary_conditions
from output import OUTPUT_DIR, METHODS, WEIGHTING_NAMES, format_evaluation, write_summary, write_details

# Resolve data/ relative to this file so the script works from any working directory
DATA_DIR = Path(__file__).resolve().parent / "data"
INSTANCES_DIR = DATA_DIR / "instances"
FIGURES_DIR = OUTPUT_DIR / "figures"

SAVE_FIGURES = True    # save the Steiner trees of each graph and distance in FIGURES_DIR


def evaluate(graph: GraphData) -> dict:
    """
    Certificate of the Steiner trees of graph.steiner_edges (one per Voronoi
    cell): structure, départs, constraints and objective J (see verify.py).
    """
    return check_solution(graph, graph.steiner_edges)


def solve(graph: GraphData, method: str) -> None:
    """Voronoi cells of the sources then Steiner tree in each cell, for one of METHODS."""
    if method == 'variable':
        compute_voronoi_variable(graph)
    else:
        compute_voronoi_cells(graph, weighting=method)
        compute_steiner_trees(graph)


def run(instance_files: list, show_plots: bool = False, label: str = None) -> list:
    """
    Voronoi + Steiner + evaluation for each instance and each distance.
    `label` is appended to the name of the reports (recap_<date>_<label>.md ...).
    """
    if SAVE_FIGURES:
        FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S") + (f"_{label}" if label else "")
    results, details = [], []
    for path in instance_files:
        G = GraphData(folderpath=str(path.parent), filename=path.name)
        print(G)
        loads = G.nodes[G.nodes['type'].isin(['Conso', 'Prod'])]
        proofs = necessary_conditions(G)
        if proofs:
            lines = [f"=== {path.stem} : instance prouvée infaisable ==="] + [f"  {p}" for p in proofs]
            print("\n".join(lines))
            details += lines + [""]
        for weighting in METHODS:
            solve(G, weighting)
            ev = evaluate(G)
            lines = format_evaluation(path.stem, weighting, ev,
                                      G.voronoi_weights if weighting == "variable" else None)
            print("\n".join(lines))
            details += lines + [""]
            results.append({'name': path.stem, 'type': path.stem.split('_')[1],
                            'n_loads': len(loads), 'n_neg': int((loads['demande'] < 0).sum()),
                            'n_sources': len(ev['sources']), 'weighting': weighting, 'ev': ev})
            if SAVE_FIGURES or show_plots:
                filepath = FIGURES_DIR / f"{path.stem}_{weighting}.png" if SAVE_FIGURES else None
                show_steiner_trees(G, title=f"{path.stem} - distance {WEIGHTING_NAMES[weighting]}"
                                            f" - J = {ev['J']:.1f}",
                                   filepath=filepath, show=show_plots)
    summary = write_summary(results, stamp)
    write_details(details, stamp)
    print(f"\nRécapitulatif écrit dans {summary}")
    if SAVE_FIGURES:
        print(f"Figures enregistrées dans {FIGURES_DIR}")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Voronoi + Steiner sur les instances de data/instances")
    parser.add_argument("motifs", nargs="*", help="ne garder que les instances dont le nom contient un de ces motifs")
    parser.add_argument("--show", action="store_true", help="ouvrir les figures")
    parser.add_argument("--label", help="suffixe du nom des rapports (recap_<date>_<label>.md ...)")
    parser.add_argument("--petites", action="store_true", help="les 20 petites instances g01 à g20")
    parser.add_argument("--grandes", action="store_true", help="les 20 grandes instances g21 à g40")
    args = parser.parse_args()
    files = [f for f in sorted(INSTANCES_DIR.glob("*.json"))
             if not args.motifs or any(m in f.stem for m in args.motifs)]
    if args.petites:
        files = [f for f in files if int(f.stem[1:3]) <= 20]
    if args.grandes:
        files = [f for f in files if int(f.stem[1:3]) > 20]
    run(files, show_plots=args.show, label=args.label)
