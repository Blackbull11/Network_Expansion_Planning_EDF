"""
PLNE (plne_reseau.build_model) sur les instances de verified_instances/,
des plus petites aux plus grandes, avec vérification par verify.check_solution.

Alignement du modèle sur verify.py (contraintes ajoutées ici, plne_reseau.py
n'est pas modifié) :
  - un feeder = un départ : un seul arc sort de la source par feeder ;
  - |P_départ| <= Pmax (le modèle n'a que P <= Pmax) ;
  - |P_source| <= capacité de la source ;
  - un feeder ne passe par aucune autre source ;
  - départs disjoints (option disjoint de build_model) ;
  - L_départ <= l_max_depart (Ltree_max) ; Lmax = l_max_depart en mode
    'shortest' (coupe valide : un chemin de l'arbre est plus court que l'arbre) ;
  - pas de budget PD (absent de verify.py).
Nombre de feeders par source = n_departs_max (borne --max-feeders possible).

Solution initiale : la solution de l'heuristique Voronoi + Steiner à distance
variable (voronoi_variable_distance), si elle est faisable, est donnée à
HiGHS (valeurs de y, z, x, w et du flot f). Sans elle, HiGHS ne trouve souvent
aucune solution entière en 10 min sur ces instances.

Chaque instance est résolue dans un sous-processus tué après
1.5 x time_limit + 120 s : HiGHS ne respecte pas toujours sa limite de temps.

Résultats : output/plne_verified/<instance>.json (statut, J, borne inférieure
de HiGHS, J de la solution initiale, temps, certificat) et
output/plne_verified/summary.csv, mis à jour après chaque instance.

Usage :
  python run_plne_verified.py --max-nodes 60            les instances de 60 nœuds au plus
  python run_plne_verified.py --n 10 --time-limit 600   les 10 plus petites
  python run_plne_verified.py --no-warm-start ...       sans solution initiale
"""
import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import pulp

from plne_reseau import build_model
from src.instances import load_instance, instance_size, verified_files
from verify import check_solution, departs_of_tree, source_limits

OUT_DIR = Path(__file__).resolve().parent / "output" / "plne_verified"
NO_PD = 1e12   # budget PD inactif


def instance_from_graph(graph, max_feeders=None) -> dict:
    """Dictionnaire `inst` de plne_reseau.build_model à partir d'un GraphData."""
    G = graph.G
    sources = [n for n, lab in G.nodes('label') if lab == 'Source']
    lims = {s: source_limits(graph, s) for s in sources}
    feeders, src_of = [], {}
    for s in sources:
        n_k = int(lims[s]['n_departs_max'])
        if max_feeders:
            n_k = min(n_k, max_feeders)
        for c in range(n_k):
            k = f"s{s}_f{c}"
            feeders.append(k)
            src_of[k] = s
    p_max = {lims[s]['p_max_depart'] for s in sources}
    l_max = {lims[s]['l_max_depart'] for s in sources}
    if len(p_max) > 1 or len(l_max) > 1:
        raise ValueError("build_model attend des p_max_depart et l_max_depart communs à toutes les sources")
    G2 = G.__class__()
    G2.add_nodes_from(G.nodes)
    G2.add_edges_from((u, v, {'L': d}) for u, v, d in G.edges.data('weight'))
    l_max = l_max.pop()
    return dict(G=G2, pos=dict(G.nodes('pos')), sources=sources,
                prosumers=[n for n, lab in G.nodes('label') if lab in ('Conso', 'Prod')],
                transship=[n for n, lab in G.nodes('label') if lab == 'Transit'],
                D={n: G.nodes[n]['conso_ete'] for n in G if G.nodes[n]['label'] in ('Conso', 'Prod')},
                feeders=feeders, src_of=src_of, capacity={s: lims[s]['capacity'] for s in sources},
                Lmax=l_max, Ltree_max=l_max, Pmax=p_max.pop(), PD=NO_PD)


def align_with_verify(m, var, inst) -> None:
    """Contraintes de verify.py absentes de build_model (voir le docstring du module)."""
    G, D, s, x, y, w = inst["G"], inst["D"], inst["src_of"], var["x"], var["y"], var["w"]
    for k in inst["feeders"]:
        r = s[k]
        m += pulp.lpSum(x[(r, j), k] for j in G.neighbors(r)) <= 1, f"one_arc_{k}"
        m += pulp.lpSum(D[v] * y[v, k] for v in inst["prosumers"]) >= -inst["Pmax"], f"cap_neg_{k}"
        for other in inst["sources"]:
            if other != r:
                m += w[other, k] == 0, f"no_src_{other}_{k}"
    for src in inst["sources"]:
        P = pulp.lpSum(D[v] * y[v, k] for k in inst["feeders"] if s[k] == src for v in inst["prosumers"])
        m += P <= inst["capacity"][src], f"capsrc_{src}"
        m += P >= -inst["capacity"][src], f"capsrc_neg_{src}"


def flow_vars(m, inst) -> dict:
    """{(arc, feeder): variable f}: build_model ne les renvoie pas, on les retrouve par leur nom f_<n>."""
    E = [tuple(sorted(e)) for e in inst["G"].edges]
    A = [(i, j) for (i, j) in E] + [(j, i) for (i, j) in E]
    keys = [(a, k) for a in A for k in inst["feeders"]]
    by_name = {v.name: v for v in m.variables()}
    return {key: by_name[f"f_{n}"] for n, key in enumerate(keys) if f"f_{n}" in by_name}


def warm_start_values(graph, inst, var, f, trees: dict) -> dict:
    """
    {variable: valeur} de la solution `trees` = {source: arêtes} (faisable pour
    verify.py) : chaque départ va dans un feeder de sa source, par longueur
    décroissante (brisure de symétrie de build_model). None si une source a
    plus de départs que de feeders.
    """
    import networkx as nx
    G = graph.G
    val = {v: 0.0 for v in list(var["y"].values()) + list(var["z"].values()) +
           list(var["x"].values()) + list(var["w"].values()) + list(f.values())}
    by_src = {r: [k for k in inst["feeders"] if inst["src_of"][k] == r] for r in inst["sources"]}
    for k in inst["feeders"]:
        val[var["w"][inst["src_of"][k], k]] = 1.0
    for r, edges in trees.items():
        departs = sorted(departs_of_tree(G, r, edges), key=lambda d: -d['length'])
        if len(departs) > len(by_src[r]):
            return None
        T = nx.Graph(edges)
        T.add_node(r)
        parent = nx.dfs_predecessors(T, r)
        size = {}
        for n in nx.dfs_postorder_nodes(T, r):
            size[n] = 1 + sum(size[c] for c in T.neighbors(n) if parent.get(c) == n)
        for d, k in zip(departs, by_src[r]):
            nodes = nx.node_connected_component(T.subgraph(set(T) - {r}), d['head'])
            for n in nodes:
                p = parent[n]
                val[var["w"][n, k]] = 1.0
                val[var["z"][tuple(sorted((p, n))), k]] = 1.0
                val[var["x"][(p, n), k]] = 1.0
                val[f[(p, n), k]] = float(size[n])
                if n in inst["D"]:
                    val[var["y"][n, k]] = 1.0
    return val


class HiGHSWarm(pulp.HiGHS):
    """pulp.HiGHS qui donne une solution initiale {variable: valeur} au solveur avant de lancer la résolution."""

    def __init__(self, start=None, **kwargs):
        super().__init__(**kwargs)
        self.start = start

    def callSolver(self, lp):
        if self.start:
            import highspy
            sol = highspy.HighsSolution()
            sol.col_value = [float(self.start.get(v, 0.0)) for v in lp.variables()]
            sol.value_valid = True
            lp.solverModel.setSolution(sol)
        lp.solverModel.run()


def heuristic_trees(graph) -> tuple:
    """(arbres, J) de la distance variable (voronoi_variable_distance), (None, None) si infaisable."""
    from voronoi_variable_distance import voronoi_variable_distance
    from voronoi import WEIGHTINGS
    from verify import solution_key
    runs = [voronoi_variable_distance(graph, init) for init in WEIGHTINGS]
    best = min(runs, key=lambda r: solution_key(r['cert']))
    return (best['trees'], best['cert']['J']) if best['cert']['feasible'] else (None, None)


def solve_graph(graph, time_limit: float, max_feeders=None, threads=None, warm: bool = True) -> tuple:
    """(résultats, arbres {source: arêtes} ou None) de la PLNE alignée sur verify.py pour `graph`."""
    inst = instance_from_graph(graph, max_feeders)
    t0 = time.time()
    m, var = build_model(inst, lmax_mode="shortest", disjoint=True)
    align_with_verify(m, var, inst)
    t_build = time.time() - t0
    res = dict(n_loads=len(inst["prosumers"]), n_feeders=len(inst["feeders"]),
               n_vars=len(m.variables()), n_cons=len(m.constraints), time_build=t_build,
               time_limit=time_limit, J_warm=None, time_warm=None)
    start = None
    if warm:
        t0 = time.time()
        trees0, res['J_warm'] = heuristic_trees(graph)
        if trees0 is not None:
            start = warm_start_values(graph, inst, var, flow_vars(m, inst), trees0)
        res['time_warm'] = time.time() - t0
    solver = HiGHSWarm(start=start, msg=False, timeLimit=time_limit, threads=threads)
    t0 = time.time()
    m.solve(solver)
    res.update(status=pulp.LpStatus[m.status], time_solve=time.time() - t0,
               obj=None, bound=None, gap=None, feasible=None, J_verify=None)
    try:   # pulp dit 'Optimal' même à la limite de temps : le statut de HiGHS fait foi
        res['bound'] = m.solverModel.getInfo().mip_dual_bound
        res['highs_status'] = m.solverModel.modelStatusToString(m.solverModel.getModelStatus())
    except Exception:
        pass
    if m.sol_status not in (pulp.LpSolutionOptimal, pulp.LpSolutionIntegerFeasible):
        return res, None
    res['obj'] = pulp.value(m.objective)
    if res['bound'] is not None and res['obj']:
        res['gap'] = (res['obj'] - res['bound']) / res['obj']
    trees = {r: set() for r in inst["sources"]}
    for k in inst["feeders"]:
        trees[inst["src_of"][k]].update(e for e in var["E"] if var["z"][e, k].value() > 0.5)
    return res, {r: sorted(es) for r, es in trees.items()}


def solve_instance(path: Path, time_limit: float, max_feeders=None, threads=None, warm: bool = True) -> dict:
    graph = load_instance(path)
    res, trees = solve_graph(graph, time_limit, max_feeders, threads, warm)
    n_nodes, n_edges, n_src = instance_size(path)
    res = dict(instance=path.stem, n_nodes=n_nodes, n_edges=n_edges, n_sources=n_src, **res)
    if trees is not None:
        cert = check_solution(graph, trees)
        res.update(feasible=cert['feasible'], J_verify=cert['J'], errors=cert['errors'],
                   violations=cert['violations'], trees={str(r): es for r, es in trees.items()})
    return res


SUMMARY_COLS = ('instance', 'n_nodes', 'n_edges', 'n_sources', 'n_loads', 'n_feeders', 'n_vars',
                'status', 'highs_status', 'obj', 'bound', 'gap', 'J_warm', 'feasible', 'J_verify', 'time_solve', 'time_limit')


def run_single(a) -> None:
    """Mode --single : une instance, résultat écrit dans OUT_DIR/<instance>.json."""
    f = Path(a.single)
    try:
        res = solve_instance(f, a.time_limit, a.max_feeders, a.threads, not a.no_warm_start)
    except Exception as exc:   # une instance qui plante ne doit pas arrêter la campagne
        res = dict(instance=f.stem, status=f"error: {exc!r}")
    with open(OUT_DIR / f"{f.stem}.json", "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, default=str)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--max-nodes", type=int, help="seulement les instances d'au plus ce nombre de nœuds")
    p.add_argument("--n", type=int, help="seulement les n plus petites instances")
    p.add_argument("--prefix", default="", help="Petit, Moyen, Grand ou TrèsGrand")
    p.add_argument("--time-limit", type=float, default=600)
    p.add_argument("--max-feeders", type=int, help="nombre maximal de feeders par source dans le modèle")
    p.add_argument("--threads", type=int)
    p.add_argument("--no-warm-start", action="store_true", help="ne pas donner la solution de l'heuristique")
    p.add_argument("--single", help=argparse.SUPPRESS)
    a = p.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if a.single:
        return run_single(a)
    files = verified_files(a.prefix)
    if a.max_nodes:
        files = [f for f in files if instance_size(f)[0] <= a.max_nodes]
    if a.n:
        files = files[:a.n]
    print(f"{len(files)} instances, limite {a.time_limit:.0f} s chacune", flush=True)
    rows = []
    for i, f in enumerate(files, 1):
        print(f"[{i}/{len(files)}] {f.stem} {instance_size(f)}", flush=True)
        cmd = [sys.executable, __file__, "--single", str(f), "--time-limit", str(a.time_limit)]
        cmd += ["--max-feeders", str(a.max_feeders)] if a.max_feeders else []
        cmd += ["--threads", str(a.threads)] if a.threads else []
        cmd += ["--no-warm-start"] if a.no_warm_start else []
        out = OUT_DIR / f"{f.stem}.json"
        out.unlink(missing_ok=True)
        t0 = time.time()
        try:
            subprocess.run(cmd, timeout=1.5 * a.time_limit + 120, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        except subprocess.TimeoutExpired:
            pass
        if out.exists():
            with open(out, encoding="utf-8") as fh:
                res = json.load(fh)
        else:
            res = dict(instance=f.stem, status="tué (temps dépassé)", time_solve=time.time() - t0)
            with open(out, "w", encoding="utf-8") as fh:
                json.dump(res, fh, indent=1)
        print(f"    statut={res.get('highs_status', res.get('status'))} J={res.get('obj')} borne={res.get('bound')} "
              f"J_heuristique={res.get('J_warm')} faisable(verify)={res.get('feasible')} "
              f"t={res.get('time_solve') or 0:.0f}s", flush=True)
        rows.append(res)
        with open(OUT_DIR / "summary.csv", "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=SUMMARY_COLS, extrasaction='ignore')
            wr.writeheader()
            wr.writerows(rows)


if __name__ == "__main__":
    main()
