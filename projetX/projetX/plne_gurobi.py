"""
PLNE de plne_reseau.py et run_plne_verified.py, résolue avec l'API Python de
Gurobi (gurobipy) au lieu de PuLP + HiGHS. Le modèle est le même : les résultats
sont comparables à ceux de HiGHS.

Deux commandes :
  demo      modèle de plne_reseau.py sur une instance aléatoire
            (options --lmax-mode tree|shortest et --disjoint) ;
  verified  modèle aligné sur verify.py (voir run_plne_verified.py) sur
            data/verified_instances, des plus petites aux plus grandes, avec
            solution initiale de l'heuristique et vérification par verify.py.

Résultats de `verified` : output/plne_gurobi/<instance>.json (mêmes champs que
output/plne_verified), <instance>.png (solution, une couleur par arbre de
départ), <instance>.log (journal de Gurobi) et summary.csv, mis à jour après
chaque instance. `demo` écrit aussi une image dans output/plne_gurobi/. Gurobi
respecte sa limite de temps : pas besoin du sous-processus de run_plne_verified.py.

Licence : le gurobipy de pip contient une licence limitée à 2000 variables et
2000 contraintes. Elle suffit pour `demo --n 10 --seed 5`, pas pour
verified_instances (plus de 7000 variables sur la plus petite), qui demande une
licence complète.

Usage :
  python plne_gurobi.py demo --n 10 --seed 5 --lmax-mode shortest --disjoint
  python plne_gurobi.py demo --n 30 --lmax-mode shortest --disjoint
  python plne_gurobi.py verified --max-nodes 60 --time-limit 600
  python plne_gurobi.py verified --prefix Petit_N60_S1_0.20_Ete --log
  python plne_gurobi.py verified --prefix Petit_N60 --tighten --dp-stop --no-warm-start

--tighten renforce le modèle sans perdre de solution optimale. --dp-stop calcule
la borne inférieure de subset_dp.py (jusqu'à 20 charges) et arrête Gurobi dès
qu'il trouve une solution de cette valeur : une solution faisable à une borne
inférieure est optimale (statut OPTIMAL_BY_DP).
"""
import argparse
import csv
import json
import logging
import sys
import time
from pathlib import Path

import gurobipy as gp
from gurobipy import GRB
import networkx as nx

from plne_reseau import check, make_instance, plot
from run_plne_verified import NO_PD, SUMMARY_COLS, heuristic_trees, instance_from_graph, warm_start_values
from src.instances import SIZE_CLASSES, instance_size, is_instance, load_instance
from verify import check_solution

HERE = Path(__file__).resolve().parent
INSTANCE_DIR = HERE / "data" / "verified_instances"
OUT_DIR = HERE / "output" / "plne_gurobi"
COLS = tuple(c for c in SUMMARY_COLS if c != 'highs_status')
STATUS_NAMES = {getattr(GRB.Status, n): n for n in dir(GRB.Status) if n.isupper()}
LICENSE_ERRORS = (GRB.Error.NO_LICENSE, GRB.Error.SIZE_LIMIT_EXCEEDED)
LICENSE_HELP = ("Licence Gurobi absente ou insuffisante pour ce modèle. La licence incluse dans pip est limitée à "
                "2000 variables et 2000 contraintes : installer une licence complète (gurobi.lic).")
logging.getLogger("gurobipy").propagate = False   # sinon le logger racine de src.outils réaffiche tout le journal
DP_TOL = 1e-4   # km : les longueurs ont 4 décimales
DP_MAX_LOADS = 20   # au-delà, subset_dp coûte trop cher (3^charges)
FEEDER_COLORS = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
MAX_ROWS = 10   # lignes du tableau des départs sur l'image


def add_vars(m, name, keys, vtype=GRB.CONTINUOUS, ub=GRB.INFINITY) -> dict:
    """{clé: variable >= 0} nommées <name>_<n>, comme plne_reseau.new_vars."""
    return {key: m.addVar(ub=ub, vtype=vtype, name=f"{name}_{n}") for n, key in enumerate(keys)}


def build_model(inst, lmax_mode="tree", disjoint=False) -> tuple:
    """(modèle Gurobi, variables) : traduction de plne_reseau.build_model, voir ce module pour les notations."""
    G, D = inst["G"], inst["D"]
    V = list(G.nodes)
    E = [tuple(sorted(e)) for e in G.edges]
    L = {e: G.edges[e]["L"] for e in E}
    A = [(i, j) for (i, j) in E] + [(j, i) for (i, j) in E]
    Vc, K, s = inst["prosumers"], inst["feeders"], inst["src_of"]
    nV = len(V)

    sp = {src: nx.single_source_dijkstra_path_length(G, src, weight="L") for src in inst["sources"]}
    delta = {(v, k): sp[s[k]][v] for v in Vc for k in K}

    m = gp.Model("network_planning")
    y = add_vars(m, "y", [(v, k) for v in Vc for k in K], GRB.BINARY)
    z = add_vars(m, "z", [(e, k) for e in E for k in K], GRB.BINARY)
    x = add_vars(m, "x", [(a, k) for a in A for k in K], GRB.BINARY)
    w = add_vars(m, "w", [(i, k) for i in V for k in K], GRB.BINARY)
    f = add_vars(m, "f", [(a, k) for a in A for k in K])
    length = {k: gp.quicksum(L[e] * z[e, k] for e in E) for k in K}

    # (2.8) objectif
    m.setObjective(gp.quicksum(length[k] for k in K), GRB.MINIMIZE)

    # (2.3) affectation
    for v in Vc:
        m.addConstr(gp.quicksum(y[v, k] for k in K) == 1, name=f"assign_{v}")

    for k in K:
        r = s[k]
        # (2.4) capacité du feeder
        m.addConstr(gp.quicksum(D[v] * y[v, k] for v in Vc) <= inst["Pmax"], name=f"cap_{k}")
        # (2.6) longueur max de l'arbre
        m.addConstr(length[k] <= inst["Ltree_max"], name=f"Ltree_{k}")
        # (2.7) budget puissance-distance (NO_PD : pas de contrainte, évite un second membre de 1e12)
        if inst["PD"] < NO_PD:
            m.addConstr(gp.quicksum(D[v] * delta[v, k] * y[v, k] for v in Vc) <= inst["PD"], name=f"PD_{k}")

        # (2.5) structure d'arbre enraciné en s(k)
        m.addConstr(w[r, k] == 1)
        for (i, j) in E:
            m.addConstr(x[(i, j), k] + x[(j, i), k] == z[(i, j), k])
        for i in V:
            inflow = gp.quicksum(x[(j, i), k] for j in G.neighbors(i))
            if i == r:
                m.addConstr(inflow == 0)
            else:
                m.addConstr(inflow == w[i, k])
                m.addConstr(gp.quicksum(f[(j, i), k] for j in G.neighbors(i))
                            - gp.quicksum(f[(i, j), k] for j in G.neighbors(i)) == w[i, k])
        for a in A:
            m.addConstr(f[a, k] <= (nV - 1) * x[a, k])
            m.addConstr(x[a, k] <= w[a[0], k])
            m.addConstr(x[a, k] <= w[a[1], k])
        m.addConstr(gp.quicksum(z[e, k] for e in E) == gp.quicksum(w[i, k] for i in V) - 1)
        for v in Vc:
            m.addConstr(y[v, k] <= w[v, k])

    # (2.1) L_max
    d = {}
    if lmax_mode == "shortest":
        for v in Vc:
            for k in K:
                if delta[v, k] > inst["Lmax"]:
                    m.addConstr(y[v, k] == 0)
    else:   # distance mesurée le long de l'arbre construit
        d = add_vars(m, "d", [(i, k) for i in V for k in K], ub=inst["Ltree_max"])
        for k in K:
            m.addConstr(d[s[k], k] == 0)
            for (i, j) in A:
                Lij = L[tuple(sorted((i, j)))]
                m.addConstr(d[j, k] >= d[i, k] + Lij - (inst["Ltree_max"] + Lij) * (1 - x[(i, j), k]))
            for v in Vc:
                m.addConstr(d[v, k] <= inst["Lmax"] + inst["Ltree_max"] * (1 - y[v, k]))
                if delta[v, k] > inst["Lmax"]:
                    m.addConstr(y[v, k] == 0)

    # brisure de symétrie entre les feeders d'une même source
    for src in inst["sources"]:
        Ks = [k for k in K if s[k] == src]
        for k1, k2 in zip(Ks, Ks[1:]):
            m.addConstr(length[k1] >= length[k2])

    if disjoint:
        for i in V:
            if i not in inst["sources"]:
                m.addConstr(gp.quicksum(w[i, k] for k in K) <= 1)

    return m, dict(y=y, z=z, x=x, w=w, f=f, d=d, E=E, L=L, delta=delta)


def align_with_verify(m, var, inst) -> None:
    """Contraintes de verify.py absentes de build_model (voir run_plne_verified.align_with_verify)."""
    G, D, s, x, y, w = inst["G"], inst["D"], inst["src_of"], var["x"], var["y"], var["w"]
    for k in inst["feeders"]:
        r = s[k]
        m.addConstr(gp.quicksum(x[(r, j), k] for j in G.neighbors(r)) <= 1, name=f"one_arc_{k}")
        m.addConstr(gp.quicksum(D[v] * y[v, k] for v in inst["prosumers"]) >= -inst["Pmax"], name=f"cap_neg_{k}")
        for other in inst["sources"]:
            if other != r:
                m.addConstr(w[other, k] == 0, name=f"no_src_{other}_{k}")
    for src in inst["sources"]:
        P = gp.quicksum(D[v] * y[v, k] for k in inst["feeders"] if s[k] == src for v in inst["prosumers"])
        m.addConstr(P <= inst["capacity"][src], name=f"capsrc_{src}")
        m.addConstr(P >= -inst["capacity"][src], name=f"capsrc_neg_{src}")


def tighten_model(m, var, inst) -> None:
    """
    Inégalités valides qui gardent toutes les solutions optimales (limiter les
    feeders au degré des sources se fait avant build_model, voir limit_feeders) :
      - f >= x : un arc ouvert porte au moins l'unité consommée par son extrémité ;
      - un nœud de transit d'un feeder a un fils (sinon on l'enlève, plus court) ;
      - une unité de flot des sources vers chaque charge sur les arcs du réseau,
        sommés sur les feeders : relaxation forte de l'arborescence de Steiner.
    """
    G, x, w, f, K = inst["G"], var["x"], var["w"], var["f"], inst["feeders"]
    sources = set(inst["sources"])
    arcs = [(i, j) for (i, j) in var["E"]] + [(j, i) for (i, j) in var["E"]]
    for k in K:
        for a in arcs:
            m.addConstr(f[a, k] >= x[a, k])
        for i in inst["transship"]:
            m.addConstr(gp.quicksum(x[(i, j), k] for j in G.neighbors(i)) >= w[i, k])
    for t in inst["prosumers"]:
        usable = [a for a in arcs if a[1] not in sources and a[0] != t]
        g = add_vars(m, f"g{t}", usable, ub=1.0)
        for a in usable:
            m.addConstr(g[a] <= gp.quicksum(x[a, k] for k in K))
        m.addConstr(gp.quicksum(g[a] for a in usable if a[0] in sources) == 1)
        for v in G.nodes:
            if v not in sources:
                m.addConstr(gp.quicksum(g[(v, j)] for j in G.neighbors(v) if (v, j) in g)
                            - gp.quicksum(g[(j, v)] for j in G.neighbors(v) if (j, v) in g) == (-1 if v == t else 0))


def limit_feeders(inst) -> None:
    """Garde au plus deg(s) feeders par source s : chaque feeder non vide a son propre voisin de s."""
    G, sources = inst["G"], set(inst["sources"])
    degree = {s: sum(j not in sources for j in G.neighbors(s)) for s in sources}
    used, keep = {}, []
    for k in inst["feeders"]:
        s = inst["src_of"][k]
        used[s] = used.get(s, 0) + 1
        if used[s] <= degree[s]:
            keep.append(k)
    inst["feeders"] = keep


def certify_with_bound(res: dict, cert: dict) -> None:
    """Passe res['status'] à OPTIMAL_BY_DP si la solution vérifiée atteint la borne de la programmation dynamique."""
    lb = res.get('dp_bound')
    if lb is not None and cert['feasible'] and cert['J'] <= lb + DP_TOL:
        res.update(status='OPTIMAL_BY_DP', bound=lb, gap=0.0)


def optimize(m, time_limit: float, threads=None, log_file=None, console=False, params=None) -> dict:
    """Lance Gurobi ; {status, time_solve, obj, bound, gap}, obj et gap à None sans solution entière."""
    m.Params.TimeLimit = time_limit
    m.Params.LogToConsole = int(console)
    if threads:
        m.Params.Threads = threads
    if log_file:
        m.Params.LogFile = str(log_file)
    for name, value in (params or {}).items():
        m.setParam(name, value)
    t0 = time.time()
    m.optimize()
    res = dict(status=STATUS_NAMES.get(m.Status, str(m.Status)), time_solve=time.time() - t0,
               obj=None, bound=None, gap=None)
    try:   # pas de borne finie si le modèle est infaisable
        if abs(m.ObjBound) < GRB.INFINITY:
            res['bound'] = m.ObjBound
    except gp.GurobiError:
        pass
    if m.SolCount > 0:
        res['obj'], res['gap'] = m.ObjVal, m.MIPGap
    return res


def solve(inst, lmax_mode="tree", disjoint=False, time_limit=300, console=False) -> tuple:
    """Équivalent de plne_reseau.solve : (résultats, {assign, trees} ou None)."""
    m, var = build_model(inst, lmax_mode, disjoint)
    m.update()
    res = dict(optimize(m, time_limit, console=console), n_vars=m.NumVars, n_cons=m.NumConstrs)
    if m.SolCount == 0:
        return res, None
    K, Vc, E = inst["feeders"], inst["prosumers"], var["E"]
    assign = {v: k for v in Vc for k in K if var["y"][v, k].X > 0.5}
    trees = {k: [e for e in E if var["z"][e, k].X > 0.5] for k in K}
    return res, dict(assign=assign, trees=trees)


def solve_graph(graph, time_limit: float, max_feeders=None, threads=None, warm: bool = True,
                log_file=None, console=False, tighten: bool = False, dp_bound=None, params=None) -> tuple:
    """
    Équivalent de run_plne_verified.solve_graph : (résultats, arbres {source: arêtes} ou None).
    tighten : modèle renforcé (tighten_model). dp_bound : borne inférieure de subset_dp ; Gurobi
    s'arrête dès qu'il trouve une solution de cette valeur, ce qui prouve son optimalité.
    """
    inst = instance_from_graph(graph, max_feeders)
    if tighten:
        limit_feeders(inst)
    t0 = time.time()
    m, var = build_model(inst, lmax_mode="shortest", disjoint=True)
    align_with_verify(m, var, inst)
    if tighten:
        tighten_model(m, var, inst)
    m.update()
    res = dict(n_loads=len(inst["prosumers"]), n_feeders=len(inst["feeders"]), n_vars=m.NumVars,
               n_cons=m.NumConstrs, time_build=time.time() - t0, time_limit=time_limit,
               J_warm=None, time_warm=None, feasible=None, J_verify=None)
    if warm:
        t0 = time.time()
        trees0, res['J_warm'] = heuristic_trees(graph)
        start = warm_start_values(graph, inst, var, var["f"], trees0) if trees0 is not None else None
        for v, val in (start or {}).items():
            v.Start = val
        res['time_warm'] = time.time() - t0
    if dp_bound is not None:
        params = dict(params or {}, BestObjStop=dp_bound + DP_TOL / 2)
        res['dp_bound'] = dp_bound
    res.update(optimize(m, time_limit, threads, log_file, console, params))
    if m.SolCount == 0:
        return res, None
    trees = {r: set() for r in inst["sources"]}
    for k in inst["feeders"]:
        trees[inst["src_of"][k]].update(e for e in var["E"] if var["z"][e, k].X > 0.5)
    return res, {r: sorted(es) for r, es in trees.items()}


def plot_solution(graph, trees: dict, cert: dict, res: dict, out: Path) -> None:
    """Écrit dans `out` l'image de la solution : une couleur par arbre de départ, chiffres de verify.py à droite."""
    from math import sqrt

    import matplotlib.patheffects as path_effects
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    import plot_demands as PD

    G = graph.G
    departs, member = [], {}   # member : {nœud : rang de son départ}
    for s, src in cert['sources'].items():
        T = nx.Graph(trees[s])
        if s in T:
            T.remove_node(s)
        for d in sorted(src['departs'], key=lambda d: -d['length']):
            member.update((n, len(departs)) for n in nx.node_connected_component(T, d['head']))
            departs.append(d)
    n_d = len(departs)
    colors = (FEEDER_COLORS[:n_d] if n_d <= len(FEEDER_COLORS)   # au-delà : teintes éloignées d'un départ au suivant
              else [plt.get_cmap('hsv')((i * 0.618) % 1) for i in range(n_d)])
    fr = lambda v, digits=2: f"{v:.{digits}f}".replace('.', ',')
    signed = lambda v: f"{v:+g}".replace('-', '−')
    pos = dict(G.nodes('pos'))
    nodes = [dict(id=n, geometry=pos[n], demande=G.nodes[n]['conso_ete'], type=lab) for n, lab in G.nodes('label')]
    sources = [n for n in nodes if n['type'] == 'Source']
    loads = sorted((n for n in nodes if n['type'] in ('Conso', 'Prod')), key=lambda n: -abs(n['demande']))
    transit = [n for n in nodes if n['type'] == 'Transit']
    used = [n for n in transit if n['id'] in member]
    rows = min(n_d, MAX_ROWS)

    k = max(1.0, sqrt(len(nodes) / PD.REF_NODES))
    fig_h = (9.0 + 0.38 * max(0, rows - 3)) * k
    f = plt.figure(figsize=(13.0 * k, fig_h))
    f.patch.set_facecolor(PD.SURFACE)
    ax = f.add_axes((0.015, 0.02, 0.665, 0.86))
    ax.add_collection(LineCollection([(pos[u], pos[v]) for u, v in G.edges], colors=PD.GRID, linewidths=0.7, zorder=1))
    segments = [[] for _ in departs]
    for s, edges in trees.items():
        for u, v in edges:
            segments[member[v if u == s else u]].append((pos[u], pos[v]))
    for seg, color in zip(segments, colors):
        ax.add_collection(LineCollection(seg, colors=[color], linewidths=3.2, capstyle='round', zorder=2))
    idle = [n['geometry'] for n in transit if n['id'] not in member]
    if idle:
        ax.scatter(*zip(*idle), s=7, color=PD.EDGE, linewidths=0, zorder=2)
    if used:
        ax.scatter(*zip(*(n['geometry'] for n in used)), s=46, c=[colors[member[n['id']]] for n in used],
                   edgecolors=PD.SURFACE, linewidths=1.5, zorder=3)
    for n in loads:
        x, y = n['geometry']
        color = colors[member[n['id']]]
        if n['demande'] == 0:
            ax.scatter([x], [y], s=70, facecolors=PD.SURFACE, edgecolors=color, linewidths=2, zorder=4)
        else:
            marker, _, area = PD.marker_of(n)
            ax.scatter([x], [y], s=area, marker=marker, color=color, edgecolors=PD.SURFACE, linewidths=1.5, zorder=4)
    ax.scatter(*zip(*(n['geometry'] for n in sources)), s=210, marker='s', color=PD.INK, edgecolors=PD.SURFACE,
               linewidths=1.5, zorder=5)
    ax.set_aspect('equal')
    ax.margins(0.05)
    ax.autoscale_view()
    ax.axis('off')
    PD.place_labels(f, ax, loads, [n['geometry'] for n in sources], [n['geometry'] for n in used])
    halo = [path_effects.withStroke(linewidth=3, foreground=PD.SURFACE)]
    for n in sources:
        ax.annotate("source", n['geometry'], xytext=(0, -17), textcoords='offset points', ha='center', va='top',
                    fontsize=10 * k, color=PD.INK, fontweight='bold', path_effects=halo, zorder=6)

    optimal = res['status'] == 'OPTIMAL'
    gap = res.get('gap')
    f.text(0.02, 0.955, f"{'Solution optimale' if optimal else 'Meilleure solution'} Gurobi : réseau de {fr(cert['J'])} km",
           fontsize=16 * k, color=PD.INK, fontweight='bold', va='center')
    f.text(0.02, 0.915, f"{res['instance']}  ·  {sum(len(e) for e in trees.values())} lignes construites sur "
                        f"{G.number_of_edges()} possibles  ·  {n_d} départs  ·  écart à l'optimum {fr(100 * (gap or 0), 1)} %  ·  "
                        "valeur écrite = demande (MW)", fontsize=10.5 * k, color=PD.INK_2, va='center')

    p = f.add_axes((0.70, 0.02, 0.285, 0.86))
    top = 100.0 * fig_h / (9.0 * k)
    p.set_xlim(0, 1)
    p.set_ylim(0, top)
    p.axis('off')
    y = top - 3

    def heading(title):
        nonlocal y
        p.text(0, y, title, fontsize=12.5 * k, color=PD.INK, fontweight='bold', va='center')
        y -= 3.2
        p.plot([0, 1], [y + 0.6, y + 0.6], color=PD.GRID, lw=k)
        y -= 2.2

    def line(label, value):
        nonlocal y
        p.text(0, y, label, fontsize=10.5 * k, color=PD.INK, va='center')
        p.text(1, y, value, fontsize=10.5 * k, color=PD.INK, ha='right', va='center', fontweight='bold')
        y -= 3.8

    heading("Départs")
    for x, title in ((0.31, "Longueur"), (0.63, "Puissance"), (0.86, "Charges")):
        p.text(x, y, title, fontsize=9.5 * k, color=PD.INK_2, ha='right', va='center')
    y -= 3.4
    for i, d in enumerate(departs[:rows]):
        n_loads = sum(1 for n in loads if member[n['id']] == i)
        p.plot([0.01, 0.06], [y, y], color=colors[i], lw=4 * k, solid_capstyle='round')
        for x, text in ((0.31, f"{fr(d['length'], 1)} km"), (0.63, f"{signed(d['power'])} MW"), (0.86, str(n_loads))):
            p.text(x, y, text, fontsize=10.5 * k, color=PD.INK, ha='right', va='center')
        y -= 4.2
    if n_d > rows:
        p.text(0, y, f"+ {n_d - rows} autres départs", fontsize=9.5 * k, color=PD.INK_2, va='center')
        y -= 3.4
    limits = next(iter(cert['sources'].values()))
    n_max = ' / '.join(f"{v:g}" for v in sorted({s['n_departs_max'] for s in cert['sources'].values()}))
    p.text(0, y, f"Limites : {limits['l_max_depart']:g} km et {limits['p_max_depart']:g} MW par départ, "
                 f"{n_max} départs par source", fontsize=9.5 * k, color=PD.INK_2, va='center')
    y -= 3
    broken = ', '.join(f"{c} ({v[0]})" for c, v in cert['violations'].items() if v[0])
    p.text(0, y, f"Non respectées : {broken}" if broken else "Toutes respectées (verify.py)", fontsize=9.5 * k,
           color=PD.INK_2, va='center')
    y -= 6

    heading("Longueur totale J")
    if res.get('J_warm') is not None:
        line("Heuristique Voronoi + Steiner", f"{fr(res['J_warm'])} km")
    line("Gurobi (optimum prouvé)" if optimal else f"Gurobi ({res['status']})", f"{fr(cert['J'])} km")
    if not optimal and res.get('bound') is not None:
        line("Borne inférieure", f"{fr(res['bound'])} km")
    if res.get('J_warm') is not None:
        saved = cert['J'] - res['J_warm']
        line("Gain", f"{fr(saved).replace('-', '−')} km ({fr(100 * saved / res['J_warm'], 1).replace('-', '−')} %)")
    line("Temps de résolution", f"{res['time_solve']:.0f} s")
    y -= 2.2

    heading("Bilan de puissance aux sources")
    conso = sum(n['demande'] for n in loads if n['type'] == 'Conso')
    prod = -sum(n['demande'] for n in loads if n['type'] == 'Prod')
    line("Consommation", f"{conso:g} MW")
    line("Production", f"{prod:g} MW")
    line("Net aux sources", f"{signed(conso - prod)} MW")
    line("Capacité des sources", f"{sum(s['capacity'] for s in cert['sources'].values()):g} MW")

    handles = [Line2D([], [], marker='s', ls='', mfc=PD.INK, mec=PD.INK, ms=11 * k, label="source")]
    if conso:
        handles.append(Line2D([], [], marker='o', ls='', mfc=PD.INK_2, mec=PD.SURFACE, ms=sqrt(PD.area(2)), label="consommateur (aire ∝ |MW|)"))
    if prod:
        handles.append(Line2D([], [], marker='^', ls='', mfc=PD.INK_2, mec=PD.SURFACE, ms=sqrt(PD.area(2, '^')), label="producteur (aire ∝ |MW|)"))
    if any(n['demande'] == 0 for n in loads):
        handles.append(Line2D([], [], marker='o', ls='', mfc=PD.SURFACE, mec=PD.INK_2, mew=2, ms=8, label="charge de 0 MW"))
    handles += [Line2D([], [], marker='o', ls='', mfc=PD.INK_2, mec=PD.SURFACE, ms=7, label="nœud de transit utilisé"),
                Line2D([], [], color=PD.INK_2, lw=3.2, label="ligne construite (couleur = départ)"),
                Line2D([], [], color=PD.EDGE, lw=1, label="ligne possible non construite")]
    f.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.698, 0.015), frameon=False, fontsize=10 * k,
             labelcolor=PD.INK, labelspacing=0.95, handlelength=2.2, handletextpad=1.0, borderaxespad=0)
    f.savefig(out, dpi=max(100, 150 / k), facecolor=PD.SURFACE)
    plt.close(f)


def solve_instance(path: Path, a) -> dict:
    """Résout le fichier `path` ; dictionnaire écrit dans OUT_DIR/<instance>.json."""
    graph = load_instance(path)
    dp_bound, t_dp = None, None
    n_loads = sum(lab in ('Conso', 'Prod') for _, lab in graph.G.nodes('label'))
    if a.dp_stop and n_loads <= DP_MAX_LOADS:
        import subset_dp
        dp = subset_dp.lower_bound(graph)
        dp_bound, t_dp = dp['bound'], dp['time']
        print(f"    borne subset_dp = {dp_bound:.4f} km ({t_dp:.1f} s)", flush=True)
    if dp_bound == float('inf'):   # relaxation infaisable : l'instance l'est aussi
        return dict(instance=path.stem, status='INFEASIBLE_BY_DP', time_dp=t_dp)
    res, trees = solve_graph(graph, a.time_limit, a.max_feeders, a.threads, not a.no_warm_start,
                             OUT_DIR / f"{path.stem.replace('è', 'e')}.log", a.log,   # Gurobi déforme le è
                             tighten=a.tighten, dp_bound=dp_bound)
    n_nodes, n_edges, n_src = instance_size(path)
    res = dict(instance=path.stem, n_nodes=n_nodes, n_edges=n_edges, n_sources=n_src, time_dp=t_dp, **res)
    if trees is not None:
        cert = check_solution(graph, trees)
        res.update(feasible=cert['feasible'], J_verify=cert['J'], errors=cert['errors'],
                   violations=cert['violations'], trees={str(r): es for r, es in trees.items()})
        certify_with_bound(res, cert)
        try:
            plot_solution(graph, trees, cert, res, OUT_DIR / f"{path.stem}.png")
        except Exception as exc:   # une image ratée ne doit pas faire perdre le résultat
            res['plot_error'] = repr(exc)
            print(f"    image non écrite : {exc!r}", flush=True)
    return res


def run_verified(a) -> None:
    files = [f for f in INSTANCE_DIR.glob(f"{a.prefix}*.json")
             if f.stem.split('_')[0] in SIZE_CLASSES and is_instance(f)]
    files.sort(key=lambda f: (instance_size(f), f.name))
    if a.max_nodes:
        files = [f for f in files if instance_size(f)[0] <= a.max_nodes]
    if a.n:
        files = files[:a.n]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"{len(files)} instances, limite {a.time_limit:.0f} s chacune", flush=True)
    rows = []
    for i, f in enumerate(files, 1):
        print(f"[{i}/{len(files)}] {f.stem} {instance_size(f)}", flush=True)
        try:
            res = solve_instance(f, a)
        except gp.GurobiError as exc:
            if exc.errno in LICENSE_ERRORS:   # même erreur sur toutes les instances suivantes
                raise
            res = dict(instance=f.stem, status=f"error: {exc!r}")
        except Exception as exc:   # une instance qui plante ne doit pas arrêter la campagne
            res = dict(instance=f.stem, status=f"error: {exc!r}")
        with open(OUT_DIR / f"{f.stem}.json", "w", encoding="utf-8") as fh:
            json.dump(res, fh, indent=1, default=str)
        print(f"    statut={res['status']} J={res.get('obj')} borne={res.get('bound')} gap={res.get('gap')} "
              f"J_heuristique={res.get('J_warm')} faisable(verify)={res.get('feasible')} "
              f"t={res.get('time_solve') or 0:.0f}s", flush=True)
        rows.append(res)
        with open(OUT_DIR / "summary.csv", "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=COLS, extrasaction='ignore')
            wr.writeheader()
            wr.writerows(rows)


def run_demo(a) -> None:
    inst = make_instance(n=a.n, seed=a.seed)
    print(f"|V|={len(inst['G'])}  |E|={inst['G'].number_of_edges()}  sources={inst['sources']}  "
          f"|Vc|={len(inst['prosumers'])}  |C|={len(inst['feeders'])}")
    res, sol = solve(inst, a.lmax_mode, a.disjoint, a.time_limit, a.log)
    print(f"statut={res['status']}  objectif={res['obj']}  borne={res['bound']}  temps={res['time_solve']:.1f}s  "
          f"vars={res['n_vars']}  contraintes={res['n_cons']}")
    if sol is None:
        return
    check(inst, sol)
    for k in inst["feeders"]:
        pv = sorted(v for v, kk in sol["assign"].items() if kk == k)
        ln = sum(inst["G"].edges[e]["L"] for e in sol["trees"][k])
        print(f"  {k:8s} charge={sum(inst['D'][v] for v in pv):3d}  longueur={ln:6.2f}  prosumers={pv}")
    out = Path(a.plot) if a.plot else OUT_DIR / f"demo_n{a.n}_seed{a.seed}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    plot(inst, sol, out, f"Solution Gurobi : longueur totale = {res['obj']:.2f}")
    print("figure :", out)


def main():
    p = argparse.ArgumentParser(description="PLNE du réseau de distribution avec Gurobi")
    sub = p.add_subparsers(dest="command", required=True)
    demo = sub.add_parser("demo", help="instance aléatoire de plne_reseau.py")
    demo.add_argument("--seed", type=int, default=3)
    demo.add_argument("--n", type=int, default=30)
    demo.add_argument("--lmax-mode", choices=["tree", "shortest"], default="tree")
    demo.add_argument("--disjoint", action="store_true")
    demo.add_argument("--time-limit", type=float, default=300)
    demo.add_argument("--plot", help="chemin de l'image (défaut : output/plne_gurobi/demo_n<n>_seed<seed>.png)")
    demo.add_argument("--log", action="store_true", help="afficher le journal de Gurobi")
    demo.set_defaults(run=run_demo)
    ver = sub.add_parser("verified", help="instances de data/verified_instances")
    ver.add_argument("--max-nodes", type=int, help="seulement les instances d'au plus ce nombre de nœuds")
    ver.add_argument("--n", type=int, help="seulement les n plus petites instances")
    ver.add_argument("--prefix", default="", help="début du nom : Petit, Moyen, Grand, TrèsGrand ou une instance")
    ver.add_argument("--time-limit", type=float, default=600)
    ver.add_argument("--max-feeders", type=int, help="nombre maximal de feeders par source dans le modèle")
    ver.add_argument("--threads", type=int)
    ver.add_argument("--no-warm-start", action="store_true", help="ne pas donner la solution de l'heuristique")
    ver.add_argument("--log", action="store_true", help="afficher le journal de Gurobi")
    ver.add_argument("--tighten", action="store_true", help="modèle renforcé (inégalités valides, voir tighten_model)")
    ver.add_argument("--dp-stop", action="store_true",
                     help=f"arrêter Gurobi dès la borne de subset_dp (jusqu'à {DP_MAX_LOADS} charges)")
    ver.set_defaults(run=run_verified)
    a = p.parse_args()
    try:
        a.run(a)
    except gp.GurobiError as exc:
        if exc.errno in LICENSE_ERRORS:
            sys.exit(f"{exc}\n{LICENSE_HELP}")
        raise


if __name__ == "__main__":
    main()
