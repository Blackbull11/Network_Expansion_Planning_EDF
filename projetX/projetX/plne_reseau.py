#!/usr/bin/env python3
"""
PLNE - Simple power distribution network expansion planning (section 2 du sujet).

Lit une instance JSON (nodes / edges / arcs / parametres), construit la PLNE
"foret de Steiner enracinee avec capacites" et la resout.

Usage minimal :
    python plne_reseau.py g15_epars_c100_s8.json

Exemples :
    python plne_reseau.py inst.json --time-limit 600 --gap 0.01
    python plne_reseau.py inst.json --ltree 120 --pd 40 --lmax 35 --lmax-exact
    python plne_reseau.py inst.json --solver gurobi --lp modele.lp
    # modele reduit a 1 depart par source, puis reprise sur le modele complet
    python plne_reseau.py inst.json --departs-max 1 --sortie sol1.json
    python plne_reseau.py inst.json --depart sol1.json

Dependances : pulp (API classique, versions 2.8 a 3.x : la 4.0 a change d'API),
highspy (solveur HiGHS), networkx ; matplotlib pour --plot.
    pip install "pulp>=2.8,<4" highspy networkx matplotlib

Notations (identiques a la formulation papier) :
    y[v,k]  binaire  : prosumer v affecte au depart (feeder) k
    z[e,k]  binaire  : arete e construite pour le depart k
    x[a,k]  binaire  : arc a utilise par le depart k   (z_(u,v)^k dans le sujet)
    f[a,k]  continue : flot sur l'arc a pour le depart k
"""
import argparse
import json
import sys
import time
from pathlib import Path

import networkx as nx
import pulp

if int(pulp.__version__.split(".")[0]) >= 4:
    sys.exit('Ce script utilise l\'API classique de PuLP : pip install "pulp>=2.8,<4"')


# --------------------------------------------------------------------------- #
# 1. Lecture de l'instance
# --------------------------------------------------------------------------- #
def lire_instance(chemin, scenario="demande", departs_max=None):
    with open(chemin) as fh:
        data = json.load(fh)
    if "noeuds" in data:                          # format data/JCgraph/json
        import jcgraph
        data = jcgraph.convertir(data, Path(chemin).stem)

    inst = {"nom": data.get("parametres", {}).get("nom", chemin)}
    noeuds = {n["id_noeud"]: n for n in data["nodes"]}

    inst["Vs"] = sorted(i for i, n in noeuds.items() if n["type"] == "Source")
    inst["Vt"] = sorted(i for i, n in noeuds.items() if n["type"] == "Transit")
    # Prosumers = tout le reste (Conso : D > 0, Prod : D < 0)
    inst["Vc"] = sorted(i for i, n in noeuds.items()
                        if n["type"] not in ("Source", "Transit"))
    inst["D"] = {v: float(noeuds[v][scenario]) for v in inst["Vc"]}
    inst["pos"] = {i: tuple(n["geometry"]) for i, n in noeuds.items()}
    inst["type"] = {i: n["type"] for i, n in noeuds.items()}

    # Aretes
    inst["E"] = [e["edge_idx"] for e in data["edges"]]
    inst["L"] = {e["edge_idx"]: float(e["length_km"]) for e in data["edges"]}
    inst["ext"] = {e["edge_idx"]: (e["start_node_idx"], e["end_node_idx"])
                   for e in data["edges"]}
    paire_vers_arete = {}
    for e, (u, v) in inst["ext"].items():
        cle = frozenset((u, v))
        if cle in paire_vers_arete:
            sys.exit(f"Aretes paralleles entre {u} et {v} : non gere.")
        paire_vers_arete[cle] = e

    # Arcs (deux orientations par arete)
    inst["A"] = [a["arc_idx"] for a in data["arcs"]]
    inst["arc"] = {a["arc_idx"]: (a["start_node_idx"], a["end_node_idx"])
                   for a in data["arcs"]}
    inst["arete_de"] = {a: paire_vers_arete[frozenset(uv)]
                        for a, uv in inst["arc"].items()}
    inst["inverse"] = {a["arc_idx"]: a["reversed_arc_idx"] for a in data["arcs"]}
    inst["entrants"] = {i: [] for i in noeuds}
    inst["sortants"] = {i: [] for i in noeuds}
    for a, (u, v) in inst["arc"].items():
        inst["sortants"][u].append(a)
        inst["entrants"][v].append(a)

    # Departs : C = union des C_s, n_departs_max departs par source
    inst["C"], inst["src"], inst["Pmax"], inst["Lmax_depart"] = [], {}, {}, {}
    inst["departs_de"] = {s: [] for s in inst["Vs"]}
    for s in inst["Vs"]:
        nb = int(noeuds[s]["n_departs_max"])
        for j in range(min(nb, departs_max) if departs_max else nb):
            k = (s, j)
            inst["C"].append(k)
            inst["src"][k] = s
            inst["Pmax"][k] = float(noeuds[s]["p_max_depart"])
            # longueur totale max du depart (l_max_depart de verify.py), None si absente
            lmd = noeuds[s].get("l_max_depart")
            inst["Lmax_depart"][k] = float(lmd) if lmd is not None else None
            inst["departs_de"][s].append(k)
    inst["cap_source"] = {s: float(noeuds[s][scenario]) for s in inst["Vs"]}

    # delta[v, s] : plus court chemin source -> noeud (section 2.8)
    G = nx.Graph()
    G.add_nodes_from(noeuds)
    for e, (u, v) in inst["ext"].items():
        G.add_edge(u, v, weight=inst["L"][e], idx=e)
    inst["G"] = G
    inst["delta"] = {s: nx.single_source_dijkstra_path_length(G, s)
                     for s in inst["Vs"]}
    return inst


# --------------------------------------------------------------------------- #
# 2. Construction du modele
# --------------------------------------------------------------------------- #
def borne_ltree(ltree, k):
    """L^tree_max du depart k : ltree est un nombre, None ou {depart: valeur}."""
    return ltree.get(k) if isinstance(ltree, dict) else ltree


def construire_modele(inst, lmax=None, ltree=None, pd=None, lmax_exact=False,
                      pmax_sym=False, cap_source=False, disjoint=False,
                      renfort=True, un_arc_source=False):
    Vs, Vc, Vt = inst["Vs"], inst["Vc"], inst["Vt"]
    E, A, C = inst["E"], inst["A"], inst["C"]
    L, D, src = inst["L"], inst["D"], inst["src"]
    ent, sor = inst["entrants"], inst["sortants"]
    arc, arete_de, delta = inst["arc"], inst["arete_de"], inst["delta"]
    n = len(Vc)                                   # big-M du couplage f / x
    S = pulp.lpSum

    m = pulp.LpProblem("expansion_reseau", pulp.LpMinimize)
    kid = {k: f"{k[0]}_{k[1]}" for k in C}        # suffixe de nom par depart

    y = {(v, k): pulp.LpVariable(f"y_{v}_{kid[k]}", cat="Binary")
         for v in Vc for k in C}
    z = {(e, k): pulp.LpVariable(f"z_{e}_{kid[k]}", cat="Binary")
         for e in E for k in C}
    x = {(a, k): pulp.LpVariable(f"x_{a}_{kid[k]}", cat="Binary")
         for a in A for k in C}
    f = {(a, k): pulp.LpVariable(f"f_{a}_{kid[k]}", lowBound=0, upBound=n)
         for a in A for k in C}

    # Objectif (2.9)
    m += S(L[e] * z[e, k] for e in E for k in C), "longueur_totale"

    # (1) Affectation (2.3)
    for v in Vc:
        m += S(y[v, k] for k in C) == 1, f"affect_{v}"

    for k in C:
        s, K = src[k], kid[k]

        # (2) Capacite du depart (2.4)
        charge = S(D[v] * y[v, k] for v in Vc)
        m += charge <= inst["Pmax"][k], f"pmax_{K}"
        if pmax_sym:                              # option : borne aussi en injection
            m += charge >= -inst["Pmax"][k], f"pmax_inf_{K}"

        # (3) Lmax : fixation par le plus court chemin (condition necessaire)
        if lmax is not None:
            for v in Vc:
                if delta[s].get(v, float("inf")) > lmax:
                    y[v, k].upBound = 0

        # (4) Couplage arcs / arete, (5) couplage flot / arc
        for e in E:
            a1 = next(a for a in sor[inst["ext"][e][0]] if arete_de[a] == e)
            m += x[a1, k] + x[inst["inverse"][a1], k] <= z[e, k], f"cpl_{e}_{K}"
        for a in A:
            m += f[a, k] <= n * x[a, k], f"bigM_{a}_{K}"
            # un arc n'est "utilise" que s'il porte au moins une unite
            m += x[a, k] <= f[a, k], f"supp_{a}_{K}"

        # (6) Emission a la source du depart
        m += S(f[a, k] for a in sor[s]) == S(y[v, k] for v in Vc), f"emis_{K}"

        # (7) Prosumers : une unite consommee SSI v est affecte a k
        for v in Vc:
            m += (S(f[a, k] for a in ent[v]) - S(f[a, k] for a in sor[v])
                  == y[v, k]), f"cons_c_{v}_{K}"
        # (8) Transit : conservation
        for v in Vt:
            m += (S(f[a, k] for a in ent[v]) - S(f[a, k] for a in sor[v])
                  == 0), f"cons_t_{v}_{K}"

        # (9) Aucun arc entrant dans une source ; (10) seule s_k emet pour k
        for s2 in Vs:
            for a in ent[s2]:
                x[a, k].upBound = 0
                f[a, k].upBound = 0
            if s2 != s:
                for a in sor[s2]:
                    x[a, k].upBound = 0
                    f[a, k].upBound = 0

        # (11) Au plus un arc entrant par noeud non source => arborescence
        for v in Vc + Vt:
            m += S(x[a, k] for a in ent[v]) <= 1, f"in1_{v}_{K}"

        # (12) Taille maximale de l'arbre (2.7)
        if borne_ltree(ltree, k) is not None:
            m += S(L[e] * z[e, k] for e in E) <= borne_ltree(ltree, k), f"ltree_{K}"

        # Option : un depart = une seule ligne quittant la source (comme verify.py,
        # ou un depart est le sous-arbre accroche a un voisin de la source)
        if un_arc_source:
            m += S(x[a, k] for a in sor[s]) <= 1, f"un_arc_{K}"

        # (13) Puissance-distance (2.8)
        if pd is not None:
            m += S(D[v] * delta[s][v] * y[v, k] for v in Vc
                   if v in delta[s]) <= pd, f"pd_{K}"

        # (14)-(15) Inegalites valides sur les noeuds de transit (2.6)
        if renfort:
            for u in Vt:
                m += (S(x[a, k] for a in ent[u])
                      <= S(x[a, k] for a in sor[u])), f"tr_out_{u}_{K}"
                for a in sor[u]:
                    w = arc[a][1]
                    m += x[a, k] <= S(x[b, k] for b in ent[u]
                                      if arc[b][0] != w), f"tr_in_{a}_{K}"
            # un prosumer affecte a k a exactement un arc entrant pour k
            for v in Vc:
                m += y[v, k] <= S(x[a, k] for a in ent[v]), f"y_in_{v}_{K}"

        # Lmax exact : distance dans l'arbre (potentiels + big-M)
        if lmax is not None and lmax_exact:
            d = {v: pulp.LpVariable(f"d_{v}_{K}", lowBound=0, upBound=lmax)
                 for v in Vc + Vt}
            for a in A:
                u, w = arc[a]
                if w in Vs or (u in Vs and u != s):
                    continue
                la = L[arete_de[a]]
                du = 0 if u == s else d[u]
                m += d[w] >= du + la - (lmax + la) * (1 - x[a, k]), f"dist_{a}_{K}"

    # Options hors sujet -----------------------------------------------------
    if cap_source:      # charge nette totale d'une source <= sa capacite
        for s in Vs:
            m += S(D[v] * y[v, k] for v in Vc for k in inst["departs_de"][s]) \
                <= inst["cap_source"][s], f"cap_src_{s}"
    if disjoint:        # departs disjoints en aretes et en noeuds
        for e in E:
            m += S(z[e, k] for k in C) <= 1, f"disj_e_{e}"
        for v in Vc + Vt:
            m += S(x[a, k] for k in C for a in ent[v]) <= 1, f"disj_v_{v}"
        for v in Vc:    # pas de transit par un prosumer d'un autre depart
            for k in C:
                m += S(x[a, k] for a in ent[v]) <= y[v, k], f"notr_{v}_{kid[k]}"

    # Cassage de symetrie : les departs d'une meme source sont interchangeables
    for s in Vs:
        ks = inst["departs_de"][s]
        for k1, k2 in zip(ks, ks[1:]):
            m += (S(y[v, k1] for v in Vc) >= S(y[v, k2] for v in Vc)), \
                f"sym_{kid[k1]}_{kid[k2]}"

    return m, {"y": y, "z": z, "x": x, "f": f}


# --------------------------------------------------------------------------- #
# 3. Heuristique gloutonne (solution initiale fournie au solveur)
# --------------------------------------------------------------------------- #
def heuristique(inst, lmax=None, ltree=None, pd=None, pmax_sym=False,
                cap_source=False):
    """Plus proche insertion avec capacites : a chaque iteration on raccorde le
    prosumer le moins cher a raccorder a un arbre existant (plus court chemin
    depuis n'importe quel noeud de l'arbre). Renvoie une solution ou None."""
    G, D, L, delta = inst["G"], inst["D"], inst["L"], inst["delta"]
    C, src = inst["C"], inst["src"]
    # graphe prive des autres sources, par source
    H = {s: G.subgraph([v for v in G if v == s or v not in inst["Vs"]])
         for s in inst["Vs"]}
    arbre = {k: nx.Graph() for k in C}
    for k in C:
        arbre[k].add_node(src[k])
    pros = {k: [] for k in C}
    charge = {k: 0.0 for k in C}
    pdist = {k: 0.0 for k in C}
    lg = {k: 0.0 for k in C}
    ch_src = {s: 0.0 for s in inst["Vs"]}
    cache = {}
    restants = set(inst["Vc"])

    def admissible(v, k, cout):
        s = src[k]
        if charge[k] + D[v] > inst["Pmax"][k] + 1e-9:
            return False
        if cap_source and ch_src[s] + D[v] > inst["cap_source"][s] + 1e-9:
            return False
        if lmax is not None and delta[s][v] > lmax:
            return False
        lt = borne_ltree(ltree, k)
        if lt is not None and lg[k] + cout > lt + 1e-9:
            return False
        if pd is not None and pdist[k] + D[v] * delta[s][v] > pd + 1e-9:
            return False
        return True

    while restants:
        meilleur = None
        for k in C:
            if k not in cache:
                cache[k] = nx.multi_source_dijkstra(H[src[k]], set(arbre[k].nodes))
            dist, chemins = cache[k]
            for v in restants:
                if v in dist and (meilleur is None or dist[v] < meilleur[0]) \
                        and admissible(v, k, dist[v]):
                    meilleur = (dist[v], v, k, chemins[v])
        if meilleur is None:
            return None                       # glouton bloque
        cout, v, k, chemin = meilleur
        nx.add_path(arbre[k], chemin)
        pros[k].append(v)
        charge[k] += D[v]
        ch_src[src[k]] += D[v]
        pdist[k] += D[v] * delta[src[k]][v]
        lg[k] += cout
        restants.discard(v)
        cache.pop(k)

    if pmax_sym and any(charge[k] < -inst["Pmax"][k] - 1e-9 for k in C):
        return None

    sol = {"instance": inst["nom"], "statut": "Heuristique", "departs": []}
    for k in C:
        if not pros[k]:
            continue
        # amelioration : ACM sur les noeuds de l'arbre puis elagage des feuilles
        T = nx.minimum_spanning_tree(H[src[k]].subgraph(arbre[k].nodes), weight="weight")
        garder = set(pros[k]) | {src[k]}
        feuilles = [u for u in T if T.degree(u) == 1 and u not in garder]
        while feuilles:
            u = feuilles.pop()
            (w,) = T.neighbors(u)
            T.remove_node(u)
            if T.degree(w) == 1 and w not in garder:
                feuilles.append(w)
        aretes = sorted(G[u][w]["idx"] for u, w in T.edges)
        sol["departs"].append({
            "source": k[0], "depart": k[1], "prosumers": sorted(pros[k]),
            "aretes": aretes,
            "longueur_km": round(sum(L[e] for e in aretes), 4),
            "charge_MW": round(charge[k], 4)})
    sol["objectif_km"] = sum(d["longueur_km"] for d in sol["departs"])
    return sol


def injecter_depart(inst, var, sol):
    """Traduit une solution (affectations + aretes) en valeurs de y, z, x, f."""
    for famille in var.values():
        for v in famille.values():
            v.setInitialValue(0)
    arc_de = {uv: a for a, uv in inst["arc"].items()}
    for dep in sol["departs"]:
        k, s = (dep["source"], dep["depart"]), dep["source"]
        T = nx.Graph([inst["ext"][e] for e in dep["aretes"]])
        T.add_node(s)
        affectes = set(dep["prosumers"])
        for v in affectes:
            var["y"][v, k].setInitialValue(1)
        for e in dep["aretes"]:
            var["z"][e, k].setInitialValue(1)
        parent = dict(nx.bfs_predecessors(T, s))
        nb = {u: int(u in affectes) for u in T}       # prosumers dans le sous-arbre
        for u in reversed(list(nx.bfs_tree(T, s))):
            if u != s:
                nb[parent[u]] += nb[u]
                a = arc_de[parent[u], u]
                var["x"][a, k].setInitialValue(1)
                var["f"][a, k].setInitialValue(nb[u])


# --------------------------------------------------------------------------- #
# 4. Resolution
# --------------------------------------------------------------------------- #
def resoudre(m, solver="highs", time_limit=600, gap=0.01, threads=None, msg=True,
             depart=False):
    """Renvoie (duree, borne inferieure ou None)."""
    opts = dict(timeLimit=time_limit, gapRel=gap, msg=msg)
    if threads:
        opts["threads"] = threads
    t0 = time.time()
    borne = None
    if solver == "highs":
        # PuLP ne transmet pas de solution initiale a HiGHS : on le fait a la main
        import highspy
        slv = pulp.HiGHS(**opts)
        slv.createAndConfigureSolver(m)
        slv.buildSolverModel(m)
        if depart:
            ini = highspy.HighsSolution()
            ini.col_value = [v.varValue or 0.0 for v in m.variables()]
            m.solverModel.setSolution(ini)
        slv.callSolver(m)
        statut, statut_sol = slv.findSolutionValues(m)
        m.assignStatus(statut, statut_sol)
        borne = m.solverModel.getInfo().mip_dual_bound
    else:
        solveurs = {"cbc": pulp.PULP_CBC_CMD, "gurobi": pulp.GUROBI,
                    "cplex": pulp.CPLEX_PY, "scip": pulp.SCIP_PY}
        m.solve(solveurs[solver](warmStart=depart, **opts))
        if solver == "gurobi":
            borne = m.solverModel.ObjBound
    return time.time() - t0, borne


def extraire_solution(inst, m, var):
    y, z = var["y"], var["z"]
    sol = {"instance": inst["nom"],
           "statut": pulp.LpStatus[m.status],
           "objectif_km": pulp.value(m.objective),
           "departs": []}
    for k in inst["C"]:
        prosumers = [v for v in inst["Vc"] if (y[v, k].value() or 0) > 0.5]
        aretes = [e for e in inst["E"] if (z[e, k].value() or 0) > 0.5]
        if not prosumers and not aretes:
            continue
        sol["departs"].append({
            "source": k[0], "depart": k[1],
            "prosumers": prosumers,
            "aretes": aretes,
            "longueur_km": round(sum(inst["L"][e] for e in aretes), 4),
            "charge_MW": round(sum(inst["D"][v] for v in prosumers), 4),
        })
    return sol


# --------------------------------------------------------------------------- #
# 5. Verification independante du modele (ne lit que y et z)
# --------------------------------------------------------------------------- #
def verifier(inst, sol, lmax=None, ltree=None, pd=None, tol=1e-6):
    erreurs, vus = [], {}
    for dep in sol["departs"]:
        s, nom = dep["source"], f"depart ({dep['source']},{dep['depart']})"
        T = nx.Graph()
        T.add_node(s)
        for e in dep["aretes"]:
            u, v = inst["ext"][e]
            T.add_edge(u, v, weight=inst["L"][e])
        if not nx.is_connected(T):
            erreurs.append(f"{nom} : aretes non connexes")
            continue
        if T.number_of_edges() != T.number_of_nodes() - 1:
            erreurs.append(f"{nom} : contient un cycle")
        autres = [t for t in inst["Vs"] if t != s and t in T]
        if autres:
            erreurs.append(f"{nom} : passe par une autre source {autres}")
        dist = nx.single_source_dijkstra_path_length(T, s)
        dep["distance_max_km"] = round(
            max((dist[v] for v in dep["prosumers"] if v in dist), default=0.0), 4)
        for v in dep["prosumers"]:
            if v not in T:
                erreurs.append(f"{nom} : prosumer {v} hors de l'arbre")
            if v in vus:
                erreurs.append(f"prosumer {v} affecte deux fois")
            vus[v] = nom
            if lmax is not None and dist.get(v, 0) > lmax + tol:
                dep.setdefault("depassements_lmax", []).append(v)
        charge = sum(inst["D"][v] for v in dep["prosumers"])
        if charge > inst["Pmax"][(s, dep["depart"])] + tol:
            erreurs.append(f"{nom} : charge {charge:.3f} > Pmax")
        lg = sum(inst["L"][e] for e in dep["aretes"])
        lt = borne_ltree(ltree, (s, dep["depart"]))
        if lt is not None and lg > lt + tol:
            erreurs.append(f"{nom} : longueur {lg:.2f} > Ltree_max")
        pdist = sum(inst["D"][v] * inst["delta"][s][v] for v in dep["prosumers"])
        dep["puissance_distance"] = round(pdist, 4)
        if pd is not None and pdist > pd + tol:
            erreurs.append(f"{nom} : puissance-distance {pdist:.2f} > PD")
    manquants = set(inst["Vc"]) - set(vus)
    if manquants:
        erreurs.append(f"prosumers non affectes : {sorted(manquants)}")
    total = sum(d["longueur_km"] for d in sol["departs"])
    if abs(total - sol["objectif_km"]) > 1e-3:
        erreurs.append(f"objectif {sol['objectif_km']:.4f} != somme {total:.4f}")
    return erreurs


# --------------------------------------------------------------------------- #
# 6. Trace
# --------------------------------------------------------------------------- #
def tracer(inst, sol, chemin):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pos = inst["pos"]
    fig, ax = plt.subplots(figsize=(11, 11))
    for e, (u, v) in inst["ext"].items():
        ax.plot(*zip(pos[u], pos[v]), color="#d9d9d9", lw=0.6, zorder=1)
    cmap = plt.get_cmap("tab10")
    styles = ["-", "--", ":"]
    rang = {s: i for i, s in enumerate(inst["Vs"])}
    for dep in sol["departs"]:
        c = cmap(rang[dep["source"]] % 10)
        for e in dep["aretes"]:
            u, v = inst["ext"][e]
            ax.plot(*zip(pos[u], pos[v]), color=c, lw=2.2,
                    ls=styles[dep["depart"] % 3], zorder=2)
        for v in dep["prosumers"]:
            prod = inst["D"][v] < 0
            ax.scatter(*pos[v], s=45, color=c, marker="^" if prod else "o",
                       edgecolor="black", linewidth=0.5, zorder=3)
    for v in inst["Vt"]:
        ax.scatter(*pos[v], s=6, color="#999999", zorder=1)
    for s in inst["Vs"]:
        ax.scatter(*pos[s], s=260, marker="s", color=cmap(rang[s] % 10),
                   edgecolor="black", linewidth=1.2, zorder=4)
        ax.annotate(str(s), pos[s], ha="center", va="center", color="white",
                    fontsize=9, fontweight="bold", zorder=5)
    ax.set_title(f"{inst['nom']} - {sol['objectif_km']:.2f} km ({sol['statut']})\n"
                 "carre = source, rond = conso, triangle = prod, "
                 "style de trait = depart")
    ax.set_aspect("equal")
    ax.set_xlabel("km")
    ax.set_ylabel("km")
    fig.tight_layout()
    fig.savefig(chemin, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("instance")
    p.add_argument("--scenario", default="demande",
                   choices=["demande", "demande_hivers"])
    p.add_argument("--departs-max", type=int,
                   help="limite le nombre de departs par source (reduit le modele)")
    p.add_argument("--lmax", type=float, help="L_max (km)")
    p.add_argument("--lmax-exact", action="store_true",
                   help="L_max sur la distance dans l'arbre (sinon plus court chemin)")
    p.add_argument("--ltree", type=float, help="L^tree_max (km)")
    p.add_argument("--pd", type=float, help="budget puissance-distance (MW.km)")
    p.add_argument("--pmax-sym", action="store_true",
                   help="impose aussi -Pmax <= charge du depart")
    p.add_argument("--cap-source", action="store_true",
                   help="charge nette d'une source <= sa capacite (hors sujet)")
    p.add_argument("--disjoint", action="store_true",
                   help="departs disjoints en aretes et en noeuds (hors sujet)")
    p.add_argument("--sans-renfort", action="store_true",
                   help="retire les inegalites valides de 2.6")
    p.add_argument("--solver", default="highs",
                   choices=["highs", "cbc", "gurobi", "cplex", "scip"])
    p.add_argument("--depart", help="solution.json d'un run precedent, "
                   "utilisee comme solution initiale")
    p.add_argument("--sans-heuristique", action="store_true",
                   help="ne fournit pas de solution initiale au solveur")
    p.add_argument("--time-limit", type=float, default=600)
    p.add_argument("--gap", type=float, default=0.01)
    p.add_argument("--threads", type=int)
    p.add_argument("--silencieux", action="store_true")
    p.add_argument("--lp", help="ecrit le modele au format .lp")
    p.add_argument("--sortie", default="solution.json")
    p.add_argument("--plot", help="fichier image de la solution")
    a = p.parse_args()

    inst = lire_instance(a.instance, a.scenario, a.departs_max)
    print(f"Instance {inst['nom']} : {len(inst['Vs'])} sources, "
          f"{len(inst['Vc'])} prosumers, {len(inst['Vt'])} transits, "
          f"{len(inst['E'])} aretes, {len(inst['C'])} departs")

    m, var = construire_modele(
        inst, lmax=a.lmax, ltree=a.ltree, pd=a.pd, lmax_exact=a.lmax_exact,
        pmax_sym=a.pmax_sym, cap_source=a.cap_source, disjoint=a.disjoint,
        renfort=not a.sans_renfort)
    print(f"Modele : {m.numVariables()} variables, {m.numConstraints()} contraintes")
    if a.lp:
        m.writeLP(a.lp)

    depart = None
    if a.depart:                              # solution d'un run precedent
        with open(a.depart) as fh:
            depart = json.load(fh)
        depart["departs"] = [d for d in depart["departs"]
                             if (d["source"], d["depart"]) in inst["src"]]
        depart["objectif_km"] = sum(d["longueur_km"] for d in depart["departs"])
        origine = f"Solution initiale {a.depart}"
    elif not a.sans_heuristique and not a.disjoint:
        depart = heuristique(inst, a.lmax, a.ltree, a.pd, a.pmax_sym, a.cap_source)
        origine = "Heuristique gloutonne"
    if depart is None:
        pass
    elif not verifier(inst, depart, a.lmax if a.lmax_exact else None, a.ltree, a.pd):
        injecter_depart(inst, var, depart)
        print(f"{origine} : {depart['objectif_km']:.2f} km")
    else:
        depart = None
        print(f"{origine} : non realisable, ignoree")

    duree, borne = resoudre(m, a.solver, a.time_limit, a.gap, a.threads,
                            msg=not a.silencieux, depart=depart is not None)
    if m.sol_status not in (pulp.LpSolutionOptimal,
                            pulp.LpSolutionIntegerFeasible):
        sys.exit(f"Pas de solution realisable trouvee (statut "
                 f"{pulp.LpStatus[m.status]}).")

    sol = extraire_solution(inst, m, var)
    sol["duree_s"] = round(duree, 1)
    if depart:
        sol["solution_initiale_km"] = round(depart["objectif_km"], 4)
    if borne is not None and borne > -1e29:
        sol["borne_inf_km"] = round(borne, 4)
        sol["gap"] = round(max(0.0, sol["objectif_km"] - borne) / sol["objectif_km"], 4)
        sol["statut"] = ("Optimal (gap <= tolerance)" if sol["gap"] <= a.gap + 1e-9
                         else "Realisable, optimalite non prouvee")
    erreurs = verifier(inst, sol, a.lmax, a.ltree, a.pd)
    sol["verification"] = erreurs or "OK"

    print(f"\nStatut {sol['statut']} - {sol['objectif_km']:.2f} km "
          f"en {duree:.0f} s - {len(sol['departs'])} departs utilises")
    if "gap" in sol:
        print(f"Borne inferieure {sol['borne_inf_km']:.2f} km - gap {100*sol['gap']:.1f} %")
    print(f"{'source':>6} {'dep':>3} {'#pros':>5} {'charge MW':>9} "
          f"{'long km':>8} {'dist max':>8} {'P.dist':>8}")
    for d in sol["departs"]:
        print(f"{d['source']:>6} {d['depart']:>3} {len(d['prosumers']):>5} "
              f"{d['charge_MW']:>9.2f} {d['longueur_km']:>8.2f} "
              f"{d.get('distance_max_km', 0):>8.2f} "
              f"{d.get('puissance_distance', 0):>8.2f}")
    print("Verification :", "OK" if not erreurs else "\n  - " + "\n  - ".join(erreurs))

    with open(a.sortie, "w") as fh:
        json.dump(sol, fh, indent=1)
    if a.plot:
        tracer(inst, sol, a.plot)


if __name__ == "__main__":
    main()
