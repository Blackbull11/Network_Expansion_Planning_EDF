"""
Génération des 40 instances de test de data/instances/ (même format que
graphe_maille.json, avec en plus, sur les sources, les attributs
p_max_depart et l_max_depart comme dans exemple_monosource) :
g01-g20 de 20 à 100 charges, g21-g40 de 200 à 500 charges.

Unités : distances en km, puissances en MW. 'demande' > 0 pour une
consommation, < 0 pour une production (type 'Prod'). Pour une source,
'demande' est sa capacité (P_max de la source).

Quatre types de répartition des charges :
  - ville    : charges concentrées autour d'un ou plusieurs centres urbains ;
  - montagne : charges le long d'une longue route de vallée (et de quelques
               vallées secondaires), graphe presque linéaire ;
  - epars    : charges dispersées uniformément sur une grande zone rurale ;
  - mixte    : ville(s) + campagne éparse autour.

Usage : python data/generate_instances.py [motifs...]
  (réécrit data/instances/*.json et leurs images data/instances/graphics/*.png,
  seulement celles dont le nom contient un des motifs s'il y en a)
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

OUT_DIR = Path(__file__).resolve().parent / "instances"
GRAPHICS_DIR = OUT_DIR / "graphics"


# --------------------------------------------------------------------------
# Graphe
# --------------------------------------------------------------------------

def gabriel_edges(P: np.ndarray) -> set:
    """Arêtes (i, j), i < j, du graphe de Gabriel des points P (planaire, connexe)."""
    D2 = ((P[:, None, :] - P[None, :, :]) ** 2).sum(-1)
    n = len(P)
    edges = set()
    for i in range(n):
        # (i, j) est une arête si aucun k n'est dans le disque de diamètre [i, j]
        s = D2[i][None, :] + D2[i + 1:]          # s[j-i-1, k] = D2[i,k] + D2[j,k]
        s[:, i] = np.inf
        s[np.arange(n - i - 1), np.arange(i + 1, n)] = np.inf
        ok = s.min(axis=1) >= D2[i, i + 1:] - 1e-9
        edges.update((i, j) for j in np.nonzero(ok)[0] + i + 1)
    return edges


def connect_components(P: np.ndarray, edges: set) -> set:
    """Ajoute les plus courtes arêtes entre composantes jusqu'à ce que le graphe soit connexe."""
    n = len(P)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for i, j in edges:
        parent[find(i)] = find(j)
    D = np.sqrt(((P[:, None, :] - P[None, :, :]) ** 2).sum(-1))
    while len({find(i) for i in range(n)}) > 1:
        comp = np.array([find(i) for i in range(n)])
        mask = comp[:, None] != comp[None, :]
        i, j = np.unravel_index(np.where(mask, D, np.inf).argmin(), D.shape)
        edges.add((min(i, j), max(i, j)))
        parent[find(i)] = find(j)
    return edges


def add_source_links(P: np.ndarray, edges: set, sources, k: int) -> set:
    """Relie chaque source à ses k plus proches voisins (plusieurs départs possibles)."""
    for s in sources:
        d = np.sqrt(((P - P[s]) ** 2).sum(-1))
        d[list(sources)] = np.inf
        for j in np.argsort(d)[:k]:
            edges.add((min(s, j), max(s, j)))
    return edges


# --------------------------------------------------------------------------
# Positions
# --------------------------------------------------------------------------

def clusters(rng, n, centers, sigmas):
    """n points répartis sur des nuages gaussiens (centres, écarts-types en km)."""
    idx = rng.integers(len(centers), size=n)
    return np.array(centers)[idx] + rng.normal(size=(n, 2)) * np.array(sigmas)[idx, None]


def random_road(rng, start, heading, length, step, wiggle):
    """Route sinueuse : marche aléatoire à pas `step` km sur `length` km."""
    pts = [np.array(start, float)]
    for _ in range(int(length / step)):
        heading += rng.normal(0, wiggle)
        pts.append(pts[-1] + step * rng.uniform(0.7, 1.3) * np.array([np.cos(heading), np.sin(heading)]))
    return np.array(pts)


# --------------------------------------------------------------------------
# Charges et sources
# --------------------------------------------------------------------------

def draw_loads(rng, n, n_neg, conso, prod):
    """
    n charges : n - n_neg consommations log-normales (médiane, sigma, max)
    et n_neg productions uniformes dans [prod[0], prod[1]] (valeurs < 0).
    """
    med, sig, cap = conso
    c = np.minimum(med * np.exp(sig * rng.normal(size=n - n_neg)), cap)
    p = rng.uniform(prod[0], prod[1], size=n_neg)
    loads = np.concatenate([c, p])
    rng.shuffle(loads)
    return np.round(loads, 2)


# --------------------------------------------------------------------------
# Scénarios
# --------------------------------------------------------------------------

def source_array(src_pos):
    """Positions des sources, None si elles sont à placer (voir kmeans_sites)."""
    return None if src_pos is None else np.array(src_pos, float)


def scenario_ville(rng, n, centers, sigmas, n_transit, src_pos=None):
    loads = clusters(rng, n, centers, sigmas)
    transit = clusters(rng, n_transit, centers, [1.4 * s for s in sigmas])
    return source_array(src_pos), loads, transit


def scenario_epars(rng, n, size, n_transit, src_pos=None):
    loads = rng.uniform(0, size, size=(n, 2))
    transit = rng.uniform(0, size, size=(n_transit, 2))
    return source_array(src_pos), loads, transit


def scenario_mixte(rng, n, frac_city, centers, sigmas, size, n_transit, src_pos=None):
    n_city = int(frac_city * n)
    loads = np.vstack([clusters(rng, n_city, centers, sigmas),
                       rng.uniform(0, size, size=(n - n_city, 2))])
    transit = np.vstack([clusters(rng, n_transit // 2, centers, [1.4 * s for s in sigmas]),
                         rng.uniform(0, size, size=(n_transit - n_transit // 2, 2))])
    return source_array(src_pos), loads, transit


def kmeans_sites(P: np.ndarray, k: int, iters: int = 30) -> np.ndarray:
    """
    k positions de sources au centre des charges P : k-moyennes (Lloyd),
    initialisées par les points les plus éloignés (déterministe).
    """
    idx = [int(np.argmin(((P - P.mean(0)) ** 2).sum(1)))]
    d = ((P - P[idx[0]]) ** 2).sum(1)
    for _ in range(k - 1):
        idx.append(int(d.argmax()))
        d = np.minimum(d, ((P - P[idx[-1]]) ** 2).sum(1))
    C = P[idx].copy()
    for _ in range(iters):
        label = ((P[:, None, :] - C[None, :, :]) ** 2).sum(-1).argmin(1)
        C = np.array([P[label == c].mean(0) if (label == c).any() else C[c] for c in range(k)])
    return C


def sized_capacities(rng, S, L, loads, margin, p_max_depart):
    """
    Capacité de chaque source = |charge nette des charges dont elle est la plus
    proche (à vol d'oiseau)| x un facteur tiré dans `margin`, au moins p_max_depart / 2.
    Les cellules de Voronoi sur le graphe ne coïncident pas avec ces zones :
    certaines sources peuvent donc être surchargées.
    """
    nearest = ((L[:, None, :] - S[None, :, :]) ** 2).sum(-1).argmin(1)
    net = np.array([loads[nearest == s].sum() for s in range(len(S))])
    return np.round(np.maximum(np.abs(net), p_max_depart / 2) * rng.uniform(*margin, len(S)), 1)


def scenario_montagne(rng, n, road_length, n_branches, branch_length, n_transit, width, src_pos=None,
                      n_src=2):
    """
    Couloir de vallée : route principale sinueuse + `n_branches` vallées secondaires.
    Les nœuds sont tirés dans une bande de demi-largeur ~`width` km autour des
    routes (fond de vallée, versants) : charges groupées en villages (80 %) ou
    isolées (fermes, refuges), nœuds de transit répartis le long du couloir.
    Les sources sont sur la route principale, régulièrement espacées.
    Le graphe (graphe de Gabriel limité au couloir, voir make_instance) est un
    ruban maillé : plusieurs chemins possibles, mais pas de raccourci à travers
    la montagne.
    """
    step = 0.5
    main = random_road(rng, (0, 0), 0.0, road_length, step, 0.12)
    roads = [main]
    heading = np.arctan2(*np.diff(main, axis=0)[:, ::-1].T)
    for _ in range(n_branches):
        i = rng.integers(len(main) // 10, len(main) - len(main) // 10)
        h = heading[min(i, len(heading) - 1)] + rng.choice([-1, 1]) * rng.uniform(1.0, 2.0)
        roads.append(random_road(rng, main[i], h, branch_length, step, 0.15)[1:])
    C = np.vstack(roads)                                  # points des routes
    road_id = np.concatenate([[k] * len(r) for k, r in enumerate(roads)])
    tangent = np.vstack([np.gradient(r, axis=0) for r in roads])
    normal = np.stack([-tangent[:, 1], tangent[:, 0]], axis=1)
    normal /= np.linalg.norm(normal, axis=1)[:, None]

    def along(idx, shift):
        """indices décalés de `shift` points le long de la même route"""
        j = np.clip(idx + shift, 0, len(C) - 1)
        return np.where(road_id[j] == road_id[idx], j, idx)

    def corridor(idx, sigma):
        return C[idx] + normal[idx] * rng.normal(0, sigma, len(idx))[:, None]

    # charges : villages (grappes de ~1-2 km) + charges isolées
    n_village = max(3, n // 7)
    villages = rng.integers(len(C), size=n_village)
    n_grouped = int(0.8 * n)
    idx = along(rng.choice(villages, n_grouped), rng.integers(-3, 4, n_grouped))
    loads = np.vstack([corridor(idx, 0.35 * width),
                       corridor(rng.integers(len(C), size=n - n_grouped), width)])
    transit = corridor(rng.integers(len(C), size=n_transit), 0.7 * width)
    src = main[[int(round(x)) for x in np.linspace(0.08, 0.92, n_src) * (len(main) - 1)]]
    return src, loads, transit


# --------------------------------------------------------------------------
# Assemblage et écriture au format graphe_maille.json
# --------------------------------------------------------------------------

def build_graph(rng, positions, edges, types, demandes, src_params, detour):
    n = len(positions)
    nodes = []
    counter = {}
    for i in range(n):
        t = types[i]
        nodes.append({"id_noeud": i, "type": t, "demande": float(demandes[i]),
                      "demande_hivers": float(demandes[i]),
                      "id_selon_type": counter.setdefault(t, 0), "n_departs_max": None,
                      "geometry": [round(float(positions[i][0]), 2), round(float(positions[i][1]), 2)],
                      "neighbors": [], "connected_edges": [],
                      "incomming_arcs": [], "outgoing_arcs": []})
        counter[t] += 1
    for s, prm in src_params.items():
        nodes[s].update(prm)
    edges = sorted((int(u), int(v)) for u, v in edges)
    m = len(edges)
    edge_list, arcs = [], []
    for e_idx, (u, v) in enumerate(edges):
        length = float(np.hypot(*(positions[u] - positions[v])) * rng.uniform(*detour))
        length = max(round(length, 2), 0.05)
        geom = [nodes[u]["geometry"], nodes[v]["geometry"]]
        edge_list.append({"edge_idx": e_idx, "length_km": length, "start_node_idx": u,
                          "end_node_idx": v, "nb_lignes_existantes": 0, "geometry": geom})
        arcs.append({"arc_idx": e_idx, "length_km": length, "start_node_idx": u,
                     "end_node_idx": v, "nb_lignes_existantes": 0, "reversed_arc_idx": e_idx + m})
        nodes[u]["neighbors"].append(v)
        nodes[v]["neighbors"].append(u)
        nodes[u]["connected_edges"].append(e_idx)
        nodes[v]["connected_edges"].append(e_idx)
        nodes[u]["outgoing_arcs"].append(e_idx)
        nodes[v]["incomming_arcs"].append(e_idx)
        nodes[v]["outgoing_arcs"].append(e_idx + m)
        nodes[u]["incomming_arcs"].append(e_idx + m)
    for e_idx, (u, v) in enumerate(edges):
        arcs.append({"arc_idx": e_idx + m, "length_km": edge_list[e_idx]["length_km"],
                     "start_node_idx": v, "end_node_idx": u, "nb_lignes_existantes": 0,
                     "reversed_arc_idx": e_idx})
    return {"nodes": nodes, "edges": edge_list, "arcs": arcs}


def plot_instance(data: dict, filepath: Path) -> None:
    """
    Image du graphe : sources en bleu, charges positives en rouge, charges
    négatives (productions) en vert, nœuds de transit en gris. La taille des
    charges est proportionnelle à |demande|.
    """
    nodes, prm = data["nodes"], data["parametres"]
    f, ax = plt.subplots(figsize=(8, 8))
    for e in data["edges"]:
        (x0, y0), (x1, y1) = e["geometry"]
        ax.plot([x0, x1], [y0, y1], color="0.75", linewidth=0.6, zorder=1)
    groups = {"Transit": [], "Conso": [], "Prod": [], "Source": []}
    for n in nodes:
        groups[n["type"]].append(n)
    xy = lambda ns: np.array([n["geometry"] for n in ns]).reshape(-1, 2)
    size = lambda ns: 15 + 25 * np.abs([n["demande"] for n in ns])
    ax.scatter(*xy(groups["Transit"]).T, s=6, color="grey", zorder=2)
    ax.scatter(*xy(groups["Conso"]).T, s=size(groups["Conso"]), color="tab:red",
               edgecolor="k", linewidth=0.3, zorder=3)
    ax.scatter(*xy(groups["Prod"]).T, s=size(groups["Prod"]), color="tab:green",
               edgecolor="k", linewidth=0.3, zorder=3)
    ax.scatter(*xy(groups["Source"]).T, s=120, marker="s", color="tab:blue",
               edgecolor="k", linewidth=0.6, zorder=4)
    sum_p = lambda ns: sum(n["demande"] for n in ns)
    legend = [
        Line2D([], [], marker="s", ls="", color="tab:blue", markersize=9,
               label=f"Poste source ({len(groups['Source'])}, {sum_p(groups['Source']):.1f} MW)"),
        Line2D([], [], marker="o", ls="", color="tab:red", markersize=7,
               label=f"Charge positive ({len(groups['Conso'])}, {sum_p(groups['Conso']):.1f} MW)"),
        Line2D([], [], marker="o", ls="", color="tab:green", markersize=7,
               label=f"Charge négative ({len(groups['Prod'])}, {sum_p(groups['Prod']):.1f} MW)"),
        Line2D([], [], marker="o", ls="", color="grey", markersize=3,
               label=f"Noeud de transit ({len(groups['Transit'])})"),
        Line2D([], [], color="0.75", label=f"Arête ({len(data['edges'])})"),
    ]
    ax.legend(handles=legend, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8,
              title="taille des charges ∝ |P|", title_fontsize=7)
    ax.set_title(f"{prm['nom']} ({prm['type']})", fontsize=11)
    ax.set_xlabel("km")
    ax.set_ylabel("km")
    ax.set_aspect("equal")
    f.savefig(filepath, dpi=150, bbox_inches="tight")
    plt.close(f)


def make_instance(name, seed, kind, n, n_neg, conso, prod, capacities, p_max_depart,
                  l_max_depart, n_departs_max, detour=(1.05, 1.3), src_links=6, n_src=None,
                  margin=None, **geo):
    """
    Construit une instance et l'écrit dans OUT_DIR/name.json.

    capacities, p_max_depart, l_max_depart, n_departs_max : une valeur par source.
    Grandes instances : capacities=None, n_src sources placées par kmeans_sites
    (sauf montagne : sur la route), capacités tirées par sized_capacities(margin),
    p_max_depart et l_max_depart identiques pour toutes les sources et
    n_departs_max = (min, max) : ceil(capacité / p_max_depart) + 1 borné à [min, max].
    """
    rng = np.random.default_rng(seed)
    if capacities is not None:
        n_src = len(capacities)
    scen = {"ville": scenario_ville, "epars": scenario_epars, "mixte": scenario_mixte,
            "montagne": scenario_montagne}[kind]
    if kind == "montagne":
        S, L, T = scen(rng, n, n_src=n_src, **geo)
        detour = (1.2, 1.6)   # routes de montagne sinueuses
        src_links = min(src_links, 4)
    else:
        S, L, T = scen(rng, n, **geo)
        if S is None:
            S = kmeans_sites(L, n_src)
    # sources (0..n_src-1), puis charges, puis transit
    P = np.vstack([S, L, T])
    src = list(range(n_src))
    edges = gabriel_edges(P)
    if kind == "montagne":
        # pas d'arête à travers la montagne (entre deux branches ou deux lacets éloignés)
        max_edge = 2.5 * geo["width"]
        edges = {(i, j) for i, j in edges if np.hypot(*(P[i] - P[j])) <= max_edge}
    edges = connect_components(P, edges)
    edges = add_source_links(P, edges, src, src_links)
    loads = draw_loads(rng, n, n_neg, conso, prod)
    if capacities is None:
        capacities = sized_capacities(rng, S, L, loads, margin, p_max_depart)
        n_departs_max = np.clip(np.ceil(capacities / p_max_depart) + 1, *n_departs_max)
        p_max_depart, l_max_depart = rep(p_max_depart, n_src), rep(l_max_depart, n_src)
    types =["Source"] * n_src + ["Conso" if x >= 0 else "Prod" for x in loads] + \
            ["Transit"] * (len(P) - n_src - n)
    demandes = np.concatenate([capacities, loads, np.zeros(len(P) - n_src - n)])
    src_params = {s: {"n_departs_max": int(n_departs_max[s]),
                      "p_max_depart": float(p_max_depart[s]),
                      "l_max_depart": float(l_max_depart[s])} for s in range(n_src)}
    data = build_graph(rng, P, edges, types, demandes, src_params, detour)
    data["parametres"] = {"nom": name, "type": kind, "seed": seed, "n_charges": n,
                          "n_charges_negatives": n_neg, "n_sources": n_src,
                          "unites": {"longueur": "km", "puissance": "MW"},
                          "conso_totale": round(float(loads[loads > 0].sum()), 2),
                          "prod_totale": round(float(-loads[loads < 0].sum()), 2),
                          "capacite_sources": [float(c) for c in capacities]}
    OUT_DIR.mkdir(exist_ok=True)
    with open(OUT_DIR / f"{name}.json", "w") as f:
        json.dump(data, f, indent=1)
    GRAPHICS_DIR.mkdir(exist_ok=True)
    plot_instance(data, GRAPHICS_DIR / f"{name}.png")
    print(f"{name:28s} {len(P):4d} noeuds {len(edges):4d} arêtes, conso {data['parametres']['conso_totale']:7.1f} MW,"
          f" prod {data['parametres']['prod_totale']:6.1f} MW, capacité sources {sum(capacities):7.1f} MW")


def rep(x, k):
    return [x] * k


# Ordres de grandeur (MW, km) :
#   ville    : charges 0.3-3 MW (quartiers, immeubles), PV/cogénération -0.3 à -1.5 MW,
#              départs courts (L_max 15-40 km) mais chargés (P_max 7-12 MW)
#   montagne : villages 0.05-0.8 MW, quelques gros clients (station de ski) 2-3 MW,
#              petite hydro -0.5 à -3 MW, départs longs (L_max 30-57 km), P_max 4-7.2 MW
#   epars    : hameaux / fermes 0.05-0.7 MW, parcs PV/éoliens -0.5 à -3 MW,
#              L_max 80-100 km, P_max 3-4 MW
# Les seuils (L_max, P_max, capacités) ont été ajustés après un premier run pour que
# les contraintes soient respectées sur environ la moitié des graphes seulement.
# Multiplier toutes les capacités d'un graphe par un même facteur ne change pas
# les cellules de Voronoi (pondérées ou non), seulement la contrainte P_max source.
INSTANCES = [
    # --- ville ---------------------------------------------------------------
    dict(name="g01_ville_c20_s2", seed=1, kind="ville", n=20, n_neg=5,
         conso=(0.9, 0.5, 3.0), prod=(-1.2, -0.3), capacities=[12.6, 6.3],
         p_max_depart=rep(7, 2), l_max_depart=rep(18, 2), n_departs_max=[4, 3],
         centers=[(5, 5)], sigmas=[2.0], n_transit=15, src_pos=[(1, 2), (8, 8)]),
    dict(name="g02_ville_c50_s3", seed=2, kind="ville", n=50, n_neg=10,
         conso=(0.9, 0.6, 3.0), prod=(-1.5, -0.3), capacities=[25, 10, 8],
         p_max_depart=rep(10, 3), l_max_depart=rep(30, 3), n_departs_max=[6, 4, 4],
         centers=[(10, 10)], sigmas=[3.0], n_transit=35, src_pos=[(10, 10), (3, 4), (17, 14)]),
    dict(name="g03_ville_c100_s4", seed=3, kind="ville", n=100, n_neg=12,
         conso=(0.7, 0.6, 3.0), prod=(-1.5, -0.3), capacities=[30, 20, 15, 10],
         p_max_depart=rep(12, 4), l_max_depart=rep(40, 4), n_departs_max=[6, 5, 4, 4],
         centers=[(12, 12), (20, 8)], sigmas=[3.5, 2.0], n_transit=60,
         src_pos=[(12, 12), (20, 8), (5, 18), (6, 4)]),
    dict(name="g04_ville_c100_s6", seed=4, kind="ville", n=100, n_neg=10,
         conso=(0.8, 0.6, 3.0), prod=(-1.5, -0.3), capacities=[25, 25, 15, 12, 10, 8],
         p_max_depart=rep(12, 6), l_max_depart=rep(40, 6), n_departs_max=[6, 6, 5, 4, 4, 3],
         centers=[(8, 8), (22, 18)], sigmas=[3.0, 3.0], n_transit=60,
         src_pos=[(8, 8), (22, 18), (15, 13), (3, 15), (27, 8), (14, 24)]),
    dict(name="g05_ville_c50_s2_desequilibre", seed=5, kind="ville", n=50, n_neg=8,
         conso=(0.8, 0.5, 2.5), prod=(-1.2, -0.3), capacities=[40, 6],
         p_max_depart=[10, 8], l_max_depart=[20, 15], n_departs_max=[8, 3],
         centers=[(10, 10)], sigmas=[3.0], n_transit=30, src_pos=[(2, 3), (10, 11)]),
    dict(name="g06_ville_c20_s3", seed=6, kind="ville", n=20, n_neg=6,
         conso=(1.2, 0.5, 3.0), prod=(-1.5, -0.5), capacities=[10, 8, 5],
         p_max_depart=rep(8, 3), l_max_depart=rep(15, 3), n_departs_max=[3, 3, 2],
         centers=[(4, 4)], sigmas=[1.5], n_transit=12, src_pos=[(4, 4), (0, 7), (8, 1)]),
    # --- montagne ------------------------------------------------------------
    dict(name="g07_montagne_c20_s2", seed=7, kind="montagne", n=20, n_neg=5,
         conso=(0.15, 0.7, 0.8), prod=(-2.5, -0.5), capacities=[6, 4],
         p_max_depart=rep(4, 2), l_max_depart=rep(30, 2), n_departs_max=[3, 3],
         road_length=35, n_branches=1, branch_length=10, n_transit=50, width=1.5),
    dict(name="g08_montagne_c50_s3", seed=8, kind="montagne", n=50, n_neg=10,
         conso=(0.2, 0.7, 2.5), prod=(-3.0, -0.5), capacities=[10, 6, 5],
         p_max_depart=rep(7.2, 3), l_max_depart=rep(42, 3), n_departs_max=[4, 3, 3],
         road_length=60, n_branches=3, branch_length=12, n_transit=100, width=1.8),
    dict(name="g09_montagne_c100_s4", seed=9, kind="montagne", n=100, n_neg=12,
         conso=(0.15, 0.8, 3.0), prod=(-3.0, -0.5), capacities=[7.2, 4.8, 4.8, 3],
         p_max_depart=rep(5, 4), l_max_depart=rep(57, 4), n_departs_max=[4, 4, 3, 3],
         road_length=90, n_branches=5, branch_length=15, n_transit=170, width=1.8),
    dict(name="g10_montagne_c50_s5", seed=10, kind="montagne", n=50, n_neg=9,
         conso=(0.2, 0.6, 2.0), prod=(-2.0, -0.5), capacities=[6, 5, 5, 4, 3],
         p_max_depart=rep(4, 5), l_max_depart=rep(30, 5), n_departs_max=[3, 3, 3, 3, 2],
         road_length=70, n_branches=4, branch_length=10, n_transit=100, width=1.6),
    dict(name="g11_montagne_c100_s8", seed=11, kind="montagne", n=100, n_neg=11,
         conso=(0.2, 0.7, 3.0), prod=(-3.0, -0.5), capacities=[8, 6, 6, 5, 5, 4, 4, 3],
         p_max_depart=rep(6.9, 8), l_max_depart=rep(35, 8), n_departs_max=rep(3, 8),
         road_length=110, n_branches=6, branch_length=12, n_transit=180, width=1.8),
    # --- épars ---------------------------------------------------------------
    dict(name="g12_epars_c20_s2", seed=12, kind="epars", n=20, n_neg=6,
         conso=(0.2, 0.6, 0.6), prod=(-2.0, -0.5), capacities=[4, 3],
         p_max_depart=rep(3, 2), l_max_depart=rep(90, 2), n_departs_max=[3, 3],
         size=40, n_transit=15, src_pos=[(10, 12), (30, 28)]),
    dict(name="g13_epars_c50_s4", seed=13, kind="epars", n=50, n_neg=10,
         conso=(0.2, 0.6, 0.6), prod=(-3.0, -0.5), capacities=[7.2, 6, 4.8, 3.6],
         p_max_depart=rep(4, 4), l_max_depart=rep(80, 4), n_departs_max=[4, 4, 3, 3],
         size=60, n_transit=40, src_pos=[(15, 15), (45, 15), (15, 45), (45, 45)]),
    dict(name="g14_epars_c100_s5", seed=14, kind="epars", n=100, n_neg=12,
         conso=(0.2, 0.6, 0.6), prod=(-3.0, -0.5), capacities=[10, 6, 5, 4, 3],
         p_max_depart=rep(4, 5), l_max_depart=rep(100, 5), n_departs_max=[5, 4, 4, 3, 3],
         size=80, n_transit=70, src_pos=[(40, 40), (15, 15), (65, 15), (15, 65), (65, 65)]),
    dict(name="g15_epars_c100_s8", seed=15, kind="epars", n=100, n_neg=10,
         conso=(0.25, 0.6, 0.7), prod=(-2.5, -0.5), capacities=[7.2, 6, 6, 4.8, 4.8, 3.6, 3.6, 2.4],
         p_max_depart=rep(3.5, 8), l_max_depart=rep(100, 8), n_departs_max=rep(3, 8),
         size=90, n_transit=70,
         src_pos=[(15, 15), (45, 15), (75, 15), (15, 75), (45, 75), (75, 75), (25, 45), (65, 45)]),
    dict(name="g16_epars_c50_s3", seed=16, kind="epars", n=50, n_neg=9,
         conso=(0.2, 0.7, 0.7), prod=(-2.5, -0.5), capacities=[7, 3, 2],
         p_max_depart=rep(3, 3), l_max_depart=rep(80, 3), n_departs_max=[4, 3, 3],
         size=70, n_transit=40, src_pos=[(35, 35), (10, 60), (60, 10)]),
    # --- mixte ---------------------------------------------------------------
    dict(name="g17_mixte_c20_s3", seed=17, kind="mixte", n=20, n_neg=6,
         conso=(0.4, 0.8, 2.5), prod=(-2.0, -0.5), capacities=[6, 3, 2],
         p_max_depart=rep(6, 3), l_max_depart=rep(30, 3), n_departs_max=[3, 3, 2],
         frac_city=0.5, centers=[(15, 15)], sigmas=[2.0], size=40, n_transit=15,
         src_pos=[(15, 15), (5, 35), (35, 5)]),
    dict(name="g18_mixte_c50_s3", seed=18, kind="mixte", n=50, n_neg=10,
         conso=(0.4, 0.8, 3.0), prod=(-3.0, -0.5), capacities=[15, 6, 4],
         p_max_depart=rep(5, 3), l_max_depart=rep(100, 3), n_departs_max=[5, 3, 3],
         frac_city=0.5, centers=[(25, 25)], sigmas=[3.0], size=60, n_transit=35,
         src_pos=[(22, 22), (50, 50), (8, 50)]),
    dict(name="g19_mixte_c100_s5", seed=19, kind="mixte", n=100, n_neg=12,
         conso=(0.4, 0.8, 3.0), prod=(-3.0, -0.5), capacities=[20, 12, 6, 5, 4],
         p_max_depart=rep(8, 5), l_max_depart=rep(60, 5), n_departs_max=[6, 5, 4, 3, 3],
         frac_city=0.55, centers=[(20, 20), (55, 50)], sigmas=[3.0, 2.5], size=70, n_transit=70,
         src_pos=[(20, 20), (55, 50), (60, 10), (10, 60), (38, 36)]),
    dict(name="g20_mixte_c100_s7", seed=20, kind="mixte", n=100, n_neg=11,
         conso=(0.35, 0.8, 3.0), prod=(-3.0, -0.5), capacities=[15, 12, 8, 5, 4, 3, 3],
         p_max_depart=rep(8, 7), l_max_depart=rep(60, 7), n_departs_max=[5, 5, 4, 3, 3, 3, 3],
         frac_city=0.5, centers=[(20, 25), (60, 30)], sigmas=[2.5, 3.0], size=80, n_transit=70,
         src_pos=[(20, 25), (60, 30), (40, 28), (10, 70), (70, 70), (75, 5), (5, 5)]),
]


# Grandes instances (200 à 500 charges), une source pour ~15 charges (comme g01-g20).
# Mêmes charges par type que ci-dessus ; les sources et leurs capacités sont placées
# automatiquement (voir make_instance), seuls les seuils par type sont fixés ici.
# Comme pour g01-g20, les seuils ont été ajustés après un premier run (avec ceux des
# petites instances, aucun des 60 runs n'était faisable) pour qu'environ la moitié
# des graphes soient faisables avec la distance non pondérée.
LARGE_DEFAULTS = {
    'ville': dict(neg=0.10, conso=(0.8, 0.6, 3.0), prod=(-1.5, -0.3), p_max_depart=21,
                  l_max_depart=40, n_departs_max=(4, 6), margin=(1.125, 1.75)),
    'montagne': dict(neg=0.11, conso=(0.2, 0.7, 3.0), prod=(-3.0, -0.5), p_max_depart=6.5,
                     l_max_depart=60, n_departs_max=(3, 4), margin=(1.125, 1.75)),
    'epars': dict(neg=0.10, conso=(0.2, 0.6, 0.6), prod=(-3.0, -0.5), p_max_depart=6,
                  l_max_depart=120, n_departs_max=(4, 6), margin=(1.125, 1.75)),
    'mixte': dict(neg=0.11, conso=(0.4, 0.8, 3.0), prod=(-3.0, -0.5), p_max_depart=15,
                  l_max_depart=90, n_departs_max=(4, 6), margin=(1.125, 1.75)),
}


def large(idx, kind, n, **geo):
    """Paramètres d'une grande instance g<idx> de type `kind` à n charges."""
    prm = dict(LARGE_DEFAULTS[kind], **geo)
    n_src = n // 15
    return dict(name=f"g{idx}_{kind}_c{n}_s{n_src}", seed=idx, kind=kind, n=n,
                n_neg=int(round(prm.pop('neg') * n)), capacities=None, n_src=n_src, **prm)


INSTANCES += [
    # --- ville : 2 à 5 centres urbains -----------------------------------------
    large(21, 'ville', 200, centers=[(12, 12), (28, 20)], sigmas=[3.5, 3.0], n_transit=120),
    large(22, 'ville', 250, centers=[(10, 10), (30, 12), (20, 28)], sigmas=[3.0, 3.0, 3.0],
          n_transit=150),
    large(23, 'ville', 300, centers=[(12, 12), (32, 14), (22, 32)], sigmas=[4.0, 3.0, 3.5],
          n_transit=180),
    large(24, 'ville', 400, centers=[(10, 10), (35, 10), (10, 35), (35, 35)],
          sigmas=[3.5, 3.0, 3.0, 4.0], n_transit=240),
    large(25, 'ville', 500, centers=[(12, 12), (40, 10), (25, 30), (10, 45), (42, 42)],
          sigmas=[4.0, 3.5, 4.5, 3.0, 3.5], n_transit=300),
    # --- montagne : longue vallée et nombreuses vallées secondaires -------------
    large(26, 'montagne', 200, road_length=120, n_branches=6, branch_length=15, n_transit=340,
          width=1.8),
    large(27, 'montagne', 250, road_length=140, n_branches=7, branch_length=15, n_transit=420,
          width=1.8),
    large(28, 'montagne', 300, road_length=160, n_branches=9, branch_length=15, n_transit=510,
          width=1.8),
    large(29, 'montagne', 400, road_length=200, n_branches=11, branch_length=18, n_transit=680,
          width=1.8),
    large(30, 'montagne', 500, road_length=240, n_branches=14, branch_length=18, n_transit=850,
          width=1.8),
    # --- épars : même densité que g14 (zone de 8.5 * sqrt(n) km de côté) --------
    large(31, 'epars', 200, size=120, n_transit=140),
    large(32, 'epars', 250, size=134, n_transit=175),
    large(33, 'epars', 300, size=147, n_transit=210),
    large(34, 'epars', 400, size=170, n_transit=280),
    large(35, 'epars', 500, size=190, n_transit=350),
    # --- mixte : villes + campagne ----------------------------------------------
    large(36, 'mixte', 200, frac_city=0.5, centers=[(30, 30), (70, 65)], sigmas=[3.0, 2.5],
          size=100, n_transit=140),
    large(37, 'mixte', 250, frac_city=0.5, centers=[(30, 30), (80, 40), (45, 85)],
          sigmas=[3.0, 2.5, 2.5], size=110, n_transit=175),
    large(38, 'mixte', 300, frac_city=0.55, centers=[(35, 35), (90, 40), (50, 95)],
          sigmas=[3.5, 3.0, 2.5], size=120, n_transit=210),
    large(39, 'mixte', 400, frac_city=0.5, centers=[(35, 35), (105, 35), (35, 105), (100, 100)],
          sigmas=[3.5, 3.0, 3.0, 3.5], size=140, n_transit=280),
    large(40, 'mixte', 500, frac_city=0.55, centers=[(40, 40), (115, 40), (75, 80), (40, 120),
                                                    (120, 120)],
          sigmas=[4.0, 3.0, 3.5, 3.0, 3.0], size=155, n_transit=350),
]


if __name__ == "__main__":
    import sys
    # python data/generate_instances.py g21 g22 : seulement les instances dont le nom contient g21 ou g22
    motifs = sys.argv[1:]
    for prm in INSTANCES:
        if not motifs or any(m in prm['name'] for m in motifs):
            make_instance(**prm)
