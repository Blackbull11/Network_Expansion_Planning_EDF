"""
PLNE : planification simplifiee d'un reseau de distribution
(foret de Steiner enracinee et capacitee).

Modele (notations du sujet) :
  y[v,k]   = 1 si le prosumer v est affecte au feeder k
  z[e,k]   = 1 si une ligne est construite sur l'arete e pour le feeder k
Variables auxiliaires pour imposer la structure d'arbre :
  x[i,j,k] = 1 si l'arc (i->j) est utilise dans l'arborescence du feeder k
  w[i,k]   = 1 si le noeud i appartient a l'arbre du feeder k
  f[i,j,k] >= 0 flot (une unite envoyee par la racine vers chaque noeud de l'arbre)
  d[i,k]   >= 0 distance le long de l'arbre entre s(k) et i (option L_max "tree")

Resolution : PuLP + HiGHS.
Usage : python3 plne_reseau.py [--seed 3] [--n 30] [--lmax-mode tree|shortest] [--disjoint]
"""
import argparse
import math
import random
import time

import networkx as nx
import pulp


# ----------------------------------------------------------------------------
# 1. Generation d'une instance
# ----------------------------------------------------------------------------
def make_instance(n=30, n_sources=3, feeders_per_source=2, frac_prosumers=0.65,
                  k_nn=4, seed=3):
    rng = random.Random(seed)
    pos = {i: (rng.random() * 10, rng.random() * 10) for i in range(n)}

    # graphe des lignes possibles : k plus proches voisins (symetrise)
    G = nx.Graph()
    G.add_nodes_from(pos)
    for i in pos:
        neigh = sorted((j for j in pos if j != i),
                       key=lambda j: math.dist(pos[i], pos[j]))[:k_nn]
        for j in neigh:
            G.add_edge(i, j, L=round(math.dist(pos[i], pos[j]), 3))
    # on rend le graphe connexe si besoin
    comps = list(nx.connected_components(G))
    while len(comps) > 1:
        a, b = comps[0], comps[1]
        i, j = min(((i, j) for i in a for j in b),
                   key=lambda p: math.dist(pos[p[0]], pos[p[1]]))
        G.add_edge(i, j, L=round(math.dist(pos[i], pos[j]), 3))
        comps = list(nx.connected_components(G))

    nodes = list(pos)
    # sources : noeuds bien espaces (farthest point sampling)
    sources = [rng.choice(nodes)]
    while len(sources) < n_sources:
        sources.append(max((v for v in nodes if v not in sources),
                           key=lambda v: min(math.dist(pos[v], pos[s]) for s in sources)))
    others = [v for v in nodes if v not in sources]
    rng.shuffle(others)
    n_pros = int(frac_prosumers * len(others))
    prosumers = sorted(others[:n_pros])
    transship = sorted(others[n_pros:])

    # demandes : surtout positives, quelques producteurs (D<0)
    D = {v: (rng.randint(-3, -1) if rng.random() < 0.15 else rng.randint(1, 5))
         for v in prosumers}

    feeders = []
    src_of = {}
    for s in sources:
        for c in range(feeders_per_source):
            k = f"s{s}_f{c}"
            feeders.append(k)
            src_of[k] = s

    total = sum(max(d, 0) for d in D.values())
    params = dict(
        Lmax=9.0,                                  # distance max prosumer-source
        Pmax=math.ceil(1.6 * total / len(feeders)),  # capacite d'un feeder
        Ltree_max=22.0,                            # longueur max d'un arbre
        PD=60.0,                                   # budget puissance x distance
    )
    return dict(G=G, pos=pos, sources=sources, prosumers=prosumers,
                transship=transship, D=D, feeders=feeders, src_of=src_of, **params)


def new_vars(m, name, keys, lowBound=None, upBound=None, cat="Continuous"):
    """Cree un dictionnaire de variables, compatible PuLP 2.x et 4.x."""
    out = {}
    for n, key in enumerate(keys):
        vname = f"{name}_{n}"
        if hasattr(m, "add_variable"):           # PuLP >= 4
            out[key] = m.add_variable(vname, lowBound=lowBound, upBound=upBound, cat=cat)
        else:                                     # PuLP 2.x
            out[key] = pulp.LpVariable(vname, lowBound=lowBound, upBound=upBound, cat=cat)
    return out


# ----------------------------------------------------------------------------
# 2. Modele PLNE
# ----------------------------------------------------------------------------
def build_model(inst, lmax_mode="tree", disjoint=False):
    G, D = inst["G"], inst["D"]
    V = list(G.nodes)
    E = [tuple(sorted(e)) for e in G.edges]
    L = {e: G.edges[e]["L"] for e in E}
    A = [(i, j) for (i, j) in E] + [(j, i) for (i, j) in E]   # arcs
    Vc, K, s = inst["prosumers"], inst["feeders"], inst["src_of"]
    nV = len(V)

    # distances plus court chemin delta[v,k] = L(s(k), v)
    sp = {src: nx.single_source_dijkstra_path_length(G, src, weight="L")
          for src in inst["sources"]}
    delta = {(v, k): sp[s[k]][v] for v in Vc for k in K}

    m = pulp.LpProblem("network_planning", pulp.LpMinimize)
    y = new_vars(m, "y", [(v, k) for v in Vc for k in K], cat="Binary")
    z = new_vars(m, "z", [(e, k) for e in E for k in K], cat="Binary")
    x = new_vars(m, "x", [(a, k) for a in A for k in K], cat="Binary")
    w = new_vars(m, "w", [(i, k) for i in V for k in K], cat="Binary")
    f = new_vars(m, "f", [(a, k) for a in A for k in K], lowBound=0)

    # (2.8) objectif
    m += pulp.lpSum(L[e] * z[e, k] for e in E for k in K)

    # (2.3) affectation
    for v in Vc:
        m += pulp.lpSum(y[v, k] for k in K) == 1, f"assign_{v}"

    for k in K:
        r = s[k]
        # (2.4) capacite du feeder
        m += pulp.lpSum(D[v] * y[v, k] for v in Vc) <= inst["Pmax"], f"cap_{k}"
        # (2.6) longueur max de l'arbre
        m += pulp.lpSum(L[e] * z[e, k] for e in E) <= inst["Ltree_max"], f"Ltree_{k}"
        # (2.7) budget puissance-distance
        m += pulp.lpSum(D[v] * delta[v, k] * y[v, k] for v in Vc) <= inst["PD"], f"PD_{k}"

        # (2.5) structure d'arbre enracine en s(k)
        m += w[r, k] == 1
        for e in E:
            i, j = e
            m += x[(i, j), k] + x[(j, i), k] == z[e, k]          # orientation
        for i in V:
            inc = [(j, i) for j in G.neighbors(i)]
            out = [(i, j) for j in G.neighbors(i)]
            if i == r:
                m += pulp.lpSum(x[a, k] for a in inc) == 0         # racine : pas d'arc entrant
            else:
                m += pulp.lpSum(x[a, k] for a in inc) == w[i, k]   # 1 parent si dans l'arbre
                # conservation : chaque noeud de l'arbre consomme 1 unite de flot
                m += (pulp.lpSum(f[a, k] for a in inc)
                      - pulp.lpSum(f[a, k] for a in out)) == w[i, k]
        for a in A:
            m += f[a, k] <= (nV - 1) * x[a, k]                     # flot seulement sur arcs ouverts
            m += x[a, k] <= w[a[0], k]                             # arc => extremites dans l'arbre
            m += x[a, k] <= w[a[1], k]
        # nb aretes = nb noeuds - 1  (avec connexite par le flot => arbre)
        m += pulp.lpSum(z[e, k] for e in E) == pulp.lpSum(w[i, k] for i in V) - 1
        # prosumer affecte => dans l'arbre
        for v in Vc:
            m += y[v, k] <= w[v, k]

    # (2.1) L_max
    if lmax_mode == "shortest":
        for v in Vc:
            for k in K:
                if delta[v, k] > inst["Lmax"]:
                    m += y[v, k] == 0
    else:  # distance mesuree le long de l'arbre construit
        # d <= Ltree_max (un chemin dans l'arbre est plus court que l'arbre) :
        # cela permet un big-M serre M_ij = Ltree_max + L_ij
        d = new_vars(m, "d", [(i, k) for i in V for k in K],
                     lowBound=0, upBound=inst["Ltree_max"])
        for k in K:
            m += d[s[k], k] == 0
            for (i, j) in A:
                Lij = L[tuple(sorted((i, j)))]
                m += d[j, k] >= d[i, k] + Lij - (inst["Ltree_max"] + Lij) * (1 - x[(i, j), k])
            for v in Vc:
                m += d[v, k] <= inst["Lmax"] + inst["Ltree_max"] * (1 - y[v, k])
                # elimination a priori (coupe valide)
                if delta[v, k] > inst["Lmax"]:
                    m += y[v, k] == 0

    # brisure de symetrie : les feeders d'une meme source sont interchangeables,
    # on impose longueur(feeder c) >= longueur(feeder c+1)
    for src in inst["sources"]:
        Ks = [k for k in K if s[k] == src]
        for k1, k2 in zip(Ks, Ks[1:]):
            m += (pulp.lpSum(L[e] * z[e, k1] for e in E)
                  >= pulp.lpSum(L[e] * z[e, k2] for e in E))

    # option : un noeud non-source appartient a au plus un feeder
    if disjoint:
        for i in V:
            if i not in inst["sources"]:
                m += pulp.lpSum(w[i, k] for k in K) <= 1

    return m, dict(y=y, z=z, x=x, w=w, E=E, L=L, delta=delta)


# ----------------------------------------------------------------------------
# 3. Resolution, verification, trace
# ----------------------------------------------------------------------------
def solve(inst, lmax_mode="tree", disjoint=False, time_limit=300):
    m, var = build_model(inst, lmax_mode, disjoint)
    t0 = time.time()
    m.solve(pulp.HiGHS(msg=False, timeLimit=time_limit))
    t = time.time() - t0
    status = pulp.LpStatus[m.status]
    res = dict(status=status, time=t, obj=pulp.value(m.objective),
               n_vars=len(m.variables()), n_cons=len(m.constraints))
    if status != "Optimal" and m.sol_status not in (1, 2):
        return res, None
    K, Vc, E = inst["feeders"], inst["prosumers"], var["E"]
    assign = {v: k for v in Vc for k in K if var["y"][v, k].value() > 0.5}
    trees = {k: [e for e in E if var["z"][e, k].value() > 0.5] for k in K}
    return res, dict(assign=assign, trees=trees)


def check(inst, sol):
    """Verification independante de la solution."""
    G, D, s = inst["G"], inst["D"], inst["src_of"]
    for k, edges in sol["trees"].items():
        load = sum(D[v] for v, kk in sol["assign"].items() if kk == k)
        assert load <= inst["Pmax"] + 1e-6, f"capacite {k}"
        if not edges:
            assert all(kk != k for kk in sol["assign"].values())
            continue
        T = nx.Graph()
        T.add_edges_from(edges)
        assert nx.is_tree(T) and s[k] in T, f"arbre {k}"
        length = sum(G.edges[e]["L"] for e in edges)
        assert length <= inst["Ltree_max"] + 1e-6, f"longueur {k}"
        for v, kk in sol["assign"].items():
            if kk == k:
                assert v in T, f"{v} hors de l'arbre {k}"
    return True


def plot(inst, sol, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    G, pos, s = inst["G"], inst["pos"], inst["src_of"]
    palette = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd",
               "#8c564b", "#e377c2", "#17becf", "#bcbd22", "#7f7f7f"]
    col = {k: palette[i % len(palette)] for i, k in enumerate(inst["feeders"])}

    fig, ax = plt.subplots(figsize=(8, 8))
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#dddddd", width=1)
    for k, edges in sol["trees"].items():
        if edges:
            nx.draw_networkx_edges(G, pos, edgelist=edges, ax=ax,
                                   edge_color=col[k], width=3, label=k)
    nx.draw_networkx_nodes(G, pos, nodelist=inst["transship"], ax=ax,
                           node_color="white", edgecolors="grey", node_size=60)
    for v in inst["prosumers"]:
        k = sol["assign"][v]
        marker = "^" if inst["D"][v] < 0 else "o"
        ax.scatter(*pos[v], s=140, c=col[k], marker=marker,
                   edgecolors="black", zorder=3)
        ax.annotate(f"{inst['D'][v]}", pos[v], xytext=(5, 5),
                    textcoords="offset points", fontsize=8)
    ax.scatter([pos[r][0] for r in inst["sources"]], [pos[r][1] for r in inst["sources"]],
               s=400, marker="s", c="black", zorder=4, label="sources")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=3)
    p.add_argument("--n", type=int, default=30)
    p.add_argument("--lmax-mode", choices=["tree", "shortest"], default="tree")
    p.add_argument("--disjoint", action="store_true")
    p.add_argument("--plot", default="solution.png")
    a = p.parse_args()

    inst = make_instance(n=a.n, seed=a.seed)
    print(f"|V|={len(inst['G'])}  |E|={inst['G'].number_of_edges()}  "
          f"sources={inst['sources']}  |Vc|={len(inst['prosumers'])}  "
          f"|C|={len(inst['feeders'])}")
    print(f"Pmax={inst['Pmax']}  Lmax={inst['Lmax']}  "
          f"Ltree_max={inst['Ltree_max']}  PD={inst['PD']}")
    res, sol = solve(inst, a.lmax_mode, a.disjoint)
    print(f"statut={res['status']}  objectif={res['obj']}  temps={res['time']:.1f}s  "
          f"vars={res['n_vars']}  contraintes={res['n_cons']}")
    if sol:
        check(inst, sol)
        for k in inst["feeders"]:
            pv = sorted(v for v, kk in sol["assign"].items() if kk == k)
            ln = sum(inst["G"].edges[e]["L"] for e in sol["trees"][k])
            ld = sum(inst["D"][v] for v in pv)
            print(f"  {k:8s} charge={ld:3d}  longueur={ln:6.2f}  prosumers={pv}")
        plot(inst, sol, a.plot,
             f"Solution optimale : longueur totale = {res['obj']:.2f}")
        print("figure :", a.plot)
