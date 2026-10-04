"""
Vérificateur de solutions : seul juge de la faisabilité d'un ensemble d'arbres
(un par source), indépendamment de l'algorithme qui les a construits.

Une solution est un dict {source: liste d'arêtes (u, v) du graphe}. Elle est
faisable si :
  Structure
    - chaque clé est un nœud 'Source' et chaque arête existe dans G ;
    - chaque arbre est un arbre (connexe, sans cycle) qui contient sa source et
      ne passe par aucune autre source ;
    - les arbres sont disjoints en nœuds (réseau radial : un nœud n'est
      alimenté que par une source) ;
    - chaque nœud Conso / Prod est dans un arbre.
  Contraintes (départ = sous-arbre accroché à un voisin de la source)
    - L_départ <= l_max_depart  (longueur totale du sous-arbre, arête depuis la source comprise)
    - |P_départ| <= p_max_depart (somme des 'demande' du sous-arbre)
    - |P_source| <= capacité de la source ('demande'), P_source calculée sur les
      nœuds de l'arbre (pas sur la cellule de Voronoi)
    - nombre de départs <= n_departs_max

check_solution renvoie aussi l'objectif J = L_tot et PD, pour que la valeur
annoncée d'une solution soit celle du certificat.

necessary_conditions donne des preuves d'infaisabilité d'une instance
(conditions nécessaires violées), avant tout calcul de solution : si la liste
est vide, l'instance n'est pas prouvée faisable pour autant.
"""
import networkx as nx
import pandas as pd

from nx_graph import GraphData

TOL = 1e-6
TERMINAL_LABELS = ('Conso', 'Prod')
CONSTRAINTS = ('L', 'P', 'source', 'n')


def source_param(graph: GraphData, s, column: str) -> float:
    """Attribute `column` of source s in graph.nodes, inf if missing (no constraint)."""
    if column not in graph.nodes.columns or pd.isna(graph.nodes.loc[s, column]):
        return float('inf')
    return float(graph.nodes.loc[s, column])


def source_limits(graph: GraphData, s) -> dict:
    """{capacity, p_max_depart, l_max_depart, n_departs_max} of source s."""
    return {'capacity': float(graph.G.nodes[s]['conso_ete']),
            **{c: source_param(graph, s, c) for c in ('p_max_depart', 'l_max_depart', 'n_departs_max')}}


def excess(value: float, limit: float) -> float:
    """Relative excess max(0, value / limit - 1), inf if limit <= 0 < value."""
    if value <= limit + TOL:
        return 0.0
    return (value - limit) / limit if limit > 0 else float('inf')


def departs_of_tree(G: nx.Graph, source, tree_edges: list) -> list:
    """
    Départs of the tree `tree_edges` rooted at `source`: one per neighbour of
    the source in the tree. Returns a list of dicts {head, length, power, pd, n_nodes}.
    The edges must form a tree containing `source` (see check_structure).

    The power flowing on an edge (parent, n) of the tree is the sum of the
    demands of the subtree of n; pd = sum over the edges of the départ of
    |power flowing on the edge| * length of the edge (MW.km).
    """
    T = nx.Graph()
    T.add_weighted_edges_from((u, v, G[u][v]['weight']) for u, v in tree_edges)
    if source not in T:
        return []
    parent = nx.dfs_predecessors(T, source)
    flow = {}  # {n: power flowing on the edge (parent[n], n)}
    for n in nx.dfs_postorder_nodes(T, source):
        if n != source:
            flow[n] = G.nodes[n]['conso_ete'] + sum(flow[c] for c in T.neighbors(n) if c != parent[n])
    heads = list(T.neighbors(source))
    T.remove_node(source)
    departs = []
    for head in heads:
        nodes = nx.node_connected_component(T, head)
        length = G[source][head]['weight'] + T.subgraph(nodes).size(weight='weight')
        pd_ = sum(abs(flow[n]) * G[parent[n]][n]['weight'] for n in nodes)
        departs.append({'head': head, 'length': length, 'power': flow[head], 'pd': pd_,
                        'n_nodes': len(nodes)})
    return departs


def check_structure(G: nx.Graph, trees: dict) -> dict:
    """
    Structural checks of `trees` (see module docstring).
    Returns {errors: list of messages, valid: {source: tree is usable},
    uncovered: loads in no tree, owner: {node: source of its tree}}.
    """
    errors, valid, owner = [], {}, {}
    for s, edges in trees.items():
        if s not in G or G.nodes[s]['label'] != 'Source':
            errors.append(f"{s} n'est pas une source")
            continue
        ok = True
        missing = [e for e in edges if not G.has_edge(*e)]
        if missing:
            errors.append(f"Source {s} : arêtes absentes du graphe {missing}")
            ok = False
        T = nx.Graph()
        T.add_node(s)
        T.add_edges_from(edges)
        if T.number_of_edges() < len(edges):
            errors.append(f"Source {s} : arêtes en double")
        if not nx.is_connected(T):
            errors.append(f"Source {s} : arbre non connexe ({nx.number_connected_components(T)} composantes)")
            ok = False
        elif T.number_of_edges() != T.number_of_nodes() - 1:
            errors.append(f"Source {s} : l'arbre contient un cycle")
            ok = False
        others = [n for n in T if n != s and n in G and G.nodes[n]['label'] == 'Source']
        if others:
            errors.append(f"Source {s} : l'arbre passe par les sources {others}")
            ok = False
        for n in T:
            if n in owner:
                errors.append(f"Noeud {n} alimenté par les sources {owner[n]} et {s}")
                ok = False
            else:
                owner[n] = s
        valid[s] = ok
    uncovered = sorted(n for n, label in G.nodes('label') if label in TERMINAL_LABELS and n not in owner)
    if uncovered:
        errors.append(f"{len(uncovered)} charges reliées à aucune source : {uncovered}")
    return {'errors': errors, 'valid': valid, 'uncovered': uncovered, 'owner': owner}


def check_solution(graph: GraphData, trees: dict) -> dict:
    """
    Certificate of the solution `trees` = {source: list of edges (u, v)}.

    Returns
    -------
    dict with
      feasible   : True iff no structural error and no violated constraint
      errors     : structural errors (see check_structure), uncovered : loads in no tree
      sources    : {s: {capacity, power, departs, pd, n_departs, n_cell, source_ok, n_ok,
                        p_max_depart, l_max_depart, n_departs_max}},
                   each départ with {head, length, power, pd, n_nodes, L_ok, P_ok}
      violations : {'L' | 'P' | 'source' | 'n': [violated, tested]}
      excess     : sum of the relative excesses of the violated constraints
      L_tot = J, PD : objective and power x distance of the trees
    The constraints of a structurally invalid tree are not evaluated (its
    source is missing from `sources`); the solution is infeasible anyway.
    """
    G = graph.G
    struct = check_structure(G, trees)
    viol = {c: [0, 0] for c in CONSTRAINTS}
    sources, total_excess, L_tot, PD = {}, 0.0, 0.0, 0.0
    for s, edges in trees.items():
        if not struct['valid'].get(s):
            continue
        lim = source_limits(graph, s)
        departs = departs_of_tree(G, s, edges)
        for d in departs:
            e_L, e_P = excess(d['length'], lim['l_max_depart']), excess(abs(d['power']), lim['p_max_depart'])
            d['L_ok'], d['P_ok'] = e_L == 0, e_P == 0
            viol['L'][0] += not d['L_ok']
            viol['P'][0] += not d['P_ok']
            total_excess += e_L + e_P
            L_tot += d['length']
            PD += d['pd']
        viol['L'][1] += len(departs)
        viol['P'][1] += len(departs)
        tree_nodes = {n for e in edges for n in e} - {s}
        power = sum(G.nodes[n]['conso_ete'] for n in tree_nodes)
        e_src, e_n = excess(abs(power), lim['capacity']), excess(len(departs), lim['n_departs_max'])
        src = {'power': power, 'departs': departs, 'pd': sum(d['pd'] for d in departs),
               'n_departs': len(departs), 'n_cell': len(tree_nodes) + 1,
               'source_ok': e_src == 0, 'n_ok': e_n == 0, **lim}
        total_excess += e_src + e_n
        viol['source'][0] += not src['source_ok']
        viol['n'][0] += not src['n_ok']
        viol['source'][1] += 1
        viol['n'][1] += 1
        sources[s] = src
    feasible = not struct['errors'] and all(v[0] == 0 for v in viol.values())
    return {'feasible': feasible, 'errors': struct['errors'], 'uncovered': struct['uncovered'],
            'sources': sources, 'violations': viol, 'excess': total_excess,
            'L_tot': L_tot, 'J': L_tot, 'PD': PD}


def solution_key(cert: dict) -> tuple:
    """
    Sort key of certificates: feasible solutions first (by J), then the
    infeasible ones by number of errors / violations, then total excess.
    """
    n_viol = sum(v[0] for v in cert['violations'].values())
    return (len(cert['errors']), n_viol, cert['excess'], cert['J'])


def necessary_conditions(graph: GraphData) -> list:
    """
    Proofs that the instance is infeasible: list of violated necessary
    conditions (empty list: no proof, the instance may still be infeasible).
      - |sum of the demands| <= sum over the sources of min(capacity, n_departs_max * p_max_depart)
        (|P_source| <= capacity and |P_source| <= sum of its |P_départ|)
      - each load i is at distance <= l_max_depart of some source s
        (L_départ >= length of the path from s to i)
    """
    G = graph.G
    sources = [n for n, label in G.nodes('label') if label == 'Source']
    loads = [n for n, label in G.nodes('label') if label in TERMINAL_LABELS]
    lims = {s: source_limits(graph, s) for s in sources}
    reasons = []
    total = sum(G.nodes[n]['conso_ete'] for n in loads)
    supply = sum(min(l['capacity'], l['n_departs_max'] * l['p_max_depart']) for l in lims.values())
    if abs(total) > supply + TOL:
        reasons.append(f"|demande totale| = {abs(total):.2f} MW > {supply:.2f} MW"
                       " = somme des min(capacité, n_departs_max x p_max_depart)")
    reach = set()
    for s, l in lims.items():
        cutoff = l['l_max_depart'] + TOL if l['l_max_depart'] < float('inf') else None
        reach.update(nx.single_source_dijkstra_path_length(G, s, cutoff=cutoff))
    far = [n for n in loads if n not in reach]
    if far:
        reasons.append(f"charges à plus de l_max_depart de toutes les sources : {far}")
    return reasons
