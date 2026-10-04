"""
Graph Voronoi cells, where the distance can be scaled by the capacity of the sources.

Each node is attributed to the source that reaches it at the smallest cost, an
edge (v, u) costing weight(v, u) / w(source of v), with w the source weight
given by the `weighting` argument (see WEIGHTINGS):
  - 'none'     : w = 1, plain graph Voronoi cells
  - 'capacity' : w = capacity of the source
  - 'sqrt'     : w = sqrt(capacity of the source)
Edges shared by two cells are split at the point at equal distance of both sources.

`compute_voronoi_cells` links this to a GraphData (see nx_graph.py).
"""
from math import sqrt
from typing import Tuple
from heapq import heappop, heappush
from itertools import count

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString, Point

from nx_graph import GraphData
from src.outils import init_stream_logger

logger = init_stream_logger()

WEIGHTINGS = ('none', 'capacity', 'sqrt')


def source_weight(G: nx.Graph, source, weighting="capacity") -> float:
    """
    Weight w of `source` (edge costs from this source are divided by w).
    `weighting` is one of WEIGHTINGS, or a dict {source: w} of free weights
    (see voronoi_variable_distance.py).
    The capacity is the node attribute 'conso_ete' (must be > 0 unless weighting='none').
    """
    if isinstance(weighting, dict):
        return weighting[source]
    if weighting == 'none':
        return 1.0
    capacity = G.nodes[source]['conso_ete']
    if weighting == 'capacity':
        return capacity
    if weighting == 'sqrt':
        return sqrt(capacity)
    raise ValueError(f"Unknown weighting {weighting!r}, expected one of {WEIGHTINGS}")


def voronoi_cells_SB(G: nx.Graph, center_nodes, weight="weight", weighting="capacity"):
    """
    Voronoi cells centered at `center_nodes`, with edge lengths divided by the
    weight (see source_weight) of the center the path comes from.

    Returns
    -------
    dist : dict
        {node: scaled distance to its closest center}
    cells : dict
        {center: set of nodes closest to it}. Nodes reachable from no center
        are grouped under the extra key "unreachable" (if any).
    paths : dict
        {node: shortest path (list of nodes) from its closest center}
    """
    dist, paths = multi_source_dijkstra_SB(G, center_nodes, weight=weight, weighting=weighting)
    return dist, cells_from_paths(G, paths), paths


def cells_from_paths(G: nx.Graph, paths: dict) -> dict:
    """Group nodes by the center their shortest path starts from."""
    nearest = {v: p[0] for v, p in paths.items()}
    cells = nx.utils.groups(nearest)
    unreachable = set(G) - set(nearest)
    if unreachable:
        cells["unreachable"] = unreachable
    return cells


def multi_source_dijkstra_SB(G, sources, weight="weight", weighting="capacity"):
    """
    Multi-source Dijkstra where the cost of an edge is its weight divided by
    the weight (see source_weight) of the source the path starts from.
    Returns ({node: distance}, {node: path from its source}).
    """
    if not sources:
        raise ValueError("sources must not be empty")
    for s in sources:
        if s not in G:
            raise nx.NodeNotFound(f"Node {s} not found in graph")
    weight = nx.weighted._weight_function(G, weight)
    paths = {source: [source] for source in sources}
    dist = _dijkstra_multisource_SB(G, sources, weight, paths, weighting)
    return dist, paths


def _dijkstra_multisource_SB(G, sources, weight, paths, weighting="capacity"):
    """
    networkx's _dijkstra_multisource, adapted to scale edge costs by the
    weight of the closest source (lines marked "# SB"). Fills `paths`.
    """
    G_succ = G._adj
    dist = {}   # final distances
    seen = {}   # best tentative distances
    c = count()  # tie-breaker so that nodes are never compared
    fringe = []
    closest_source = {}  # SB
    for source in sources:
        seen[source] = 0
        closest_source[source] = source  # SB
        heappush(fringe, (0, next(c), source))
    while fringe:
        d, _, v = heappop(fringe)
        if v in dist:
            continue
        dist[v] = d
        closest_s = closest_source[v]  # SB
        w_s = source_weight(G, closest_s, weighting)  # SB
        for u, e in G_succ[v].items():
            vu_dist = dist[v] + weight(v, u, e) / w_s  # SB
            if u in dist:
                if vu_dist < dist[u]:
                    raise ValueError("Contradictory paths found:", "negative weights?")
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                closest_source[u] = closest_s
                heappush(fringe, (vu_dist, next(c), u))
                paths[u] = paths[v] + [u]
    return dist


def voronoi_cells_edges(G: nx.Graph, cells: dict, dist: dict,
                        weighting: str = "capacity") -> Tuple[dict, list, dict]:
    """
    Put edges in cells. An edge (u, v) with u and v in different cells is split
    at the point q equidistant (scaled, see source_weight) from both sources:
    dist[u] + t*d/w_u = dist[v] + (1-t)*d/w_v, with q at t*d from u.

    `G` is modified in place: for each split edge, a node q (label 'dummy',
    attribute pos) and edges (u, q), (v, q) (label 'dummy') are added; the
    initial edge is kept. New node ids start at len(G), so G's nodes must be
    0..len(G)-1. Work on a copy of G if the original must be preserved.

    Parameters
    ----------
    cells : dict
        {source: set of nodes}, without the "unreachable" key.
    dist : dict
        {node: scaled distance to its source}, computed with the same `weighting`.

    Returns
    -------
    edge_cells : dict
        {source: list of edges (u, v) or half edges (w, q)}, plus the key
        'split' listing the initial edges that have been split.
    nodes_to_add : list
        one (q, pos_q, u, v, length(u, q), length(v, q)) tuple per split edge.
    dist_limit_nodes : dict
        {q: {cell of u: distance to q, cell of v: distance to q}} (equal values).
    """
    if "unreachable" in cells:
        raise ValueError("Some nodes are unreachable from every source: "
                         f"{sorted(cells['unreachable'])}")
    edge_cells = {k: [] for k in cells}
    edge_cells['split'] = []
    node_cell = {n: k for k in cells for n in cells[k]}

    new_node = len(G)
    nodes_to_add = []
    dist_limit_nodes = {}
    for u, v, d in G.edges.data('weight'):
        if node_cell[u] == node_cell[v]:
            edge_cells[node_cell[u]].append((u, v))
            continue
        w_u = source_weight(G, node_cell[u], weighting)
        w_v = source_weight(G, node_cell[v], weighting)
        # relative position of q on [u, v]
        t = ((dist[v] - dist[u]) * w_v * w_u / d + w_u) / (w_v + w_u) if d != 0 else 1
        if t < -1e-9 or t > 1 + 1e-9:
            raise ValueError(f"Split point of edge ({u}, {v}) is outside the edge (t={t})")
        t = min(max(t, 0), 1)
        pos_u, pos_v = np.array(G.nodes[u]['pos']), np.array(G.nodes[v]['pos'])
        q = new_node
        new_node += 1
        nodes_to_add.append((q, pos_u + t * (pos_v - pos_u), u, v, d * t, d * (1 - t)))
        edge_cells[node_cell[u]].append((u, q))
        edge_cells[node_cell[v]].append((v, q))
        edge_cells['split'].append((u, v))
        dist_limit_nodes[q] = {node_cell[u]: dist[u] + d * t / w_u,
                               node_cell[v]: dist[v] + d * (1 - t) / w_v}
    for q, pos_q, u, v, l1, l2 in nodes_to_add:
        G.add_node(q, label='dummy', pos=pos_q)
        G.add_edge(u, q, label='dummy', weight=l1)
        G.add_edge(v, q, label='dummy', weight=l2)
    return edge_cells, nodes_to_add, dist_limit_nodes


def edges_on_shortest_paths(cells: dict, paths: dict) -> dict:
    """
    {source: list of edges (u, v), u < v, on the shortest paths from the nodes
    of its cell to the source}.
    """
    used = {k: set() for k in cells}
    for k, nodes in cells.items():
        for n in nodes:
            p = paths[n]
            used[k].update((min(a, b), max(a, b)) for a, b in zip(p, p[1:]))
    return {k: sorted(es) for k, es in used.items()}


def add_split_nodes_to_dataframes(nodes: gpd.GeoDataFrame, edges: gpd.GeoDataFrame,
                                  nodes_to_add: list) -> Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """
    Copies of `nodes` and `edges` with the 'dummy' nodes q and the half edges
    (u, q), (v, q) added, so that split edges can be plotted.
    Assumes the index of `nodes` is 0..n-1 (as the ids given by voronoi_cells_edges).
    """
    new_nodes, new_edges = [], []
    for q, pos_q, u, v, l1, l2 in nodes_to_add:
        new_nodes.append({'id_noeud': q, 'type': 'dummy', 'geometry': Point(pos_q)})
        for w, length in ((u, l1), (v, l2)):
            new_edges.append({'length_km': length, 'start_node_idx': w, 'end_node_idx': q,
                              'nb_lignes_existantes': None,
                              'geometry': LineString([nodes.geometry[w], Point(pos_q)])})
    nodes1 = pd.concat([nodes, gpd.GeoDataFrame(new_nodes, geometry='geometry')], ignore_index=True)
    edges1 = pd.concat([edges, gpd.GeoDataFrame(new_edges, geometry='geometry')], ignore_index=True)
    return nodes1, edges1


def compute_voronoi_cells(graph: GraphData, weighting="capacity") -> None:
    """
    Compute the Voronoi cells of the sources of `graph`, with the distance
    scaled according to `weighting` (one of WEIGHTINGS or {source: w}, see source_weight),
    and store the results in the GraphData:
      graph.cells, graph.edge_cells, graph.used_edges_by_cells   (see nx_graph.py)
      graph.dist, graph.paths                                    (from voronoi_cells_SB)
      graph.G_split                : copy of graph.G with the split edges
      graph.nodes_split, graph.edges_split : GeoDataFrames of G_split, for plotting
    graph.G itself is left unchanged.
    """
    sources = [n for n, label in graph.G.nodes('label') if label == 'Source']
    graph.dist, graph.cells, graph.paths = voronoi_cells_SB(graph.G, sources, weighting=weighting)
    graph.used_edges_by_cells = edges_on_shortest_paths(graph.cells, graph.paths)
    graph.G_split = graph.G.copy()
    graph.edge_cells, nodes_to_add, _ = voronoi_cells_edges(graph.G_split, graph.cells, graph.dist,
                                                            weighting)
    graph.nodes_split, graph.edges_split = add_split_nodes_to_dataframes(
        graph.nodes, graph.edges, nodes_to_add)
    logger.info("Cellules de Voronoi calculées (pondération %s) : %s",
                weighting if isinstance(weighting, str) else "variable",
                {k: len(v) for k, v in graph.cells.items()})
