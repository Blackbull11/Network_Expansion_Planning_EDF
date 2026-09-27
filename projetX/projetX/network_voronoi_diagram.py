#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 14:47:19 2026

@author: SB
"""
from typing import Tuple
import networkx as nx
import matplotlib
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import json
from pathlib import Path
import os
import geopandas as gpd

from shapely.geometry import Point, LineString
from src.outils import init_stream_logger

from heapq import heappop, heappush
from itertools import count, combinations
from collections import namedtuple

logger = init_stream_logger()

# TODO : also used in pluton/main.py, to be put in io/utils file ?
def lire_json_graphe(folderpath:str):
    """
    Read a json file containing graph data.

    Parameters
    ----------
    folderpath : str
        path of the folder containing the file graphe_maille.json
    Raises
    ------
    FileNotFoundError
        Exception if file not found.

    Returns
    -------
    nodes : geopandas.GeoDataFrame
        the nodes of the graph
    edges : geopandas.GeoDataFrame
        the edges of the graph (undirected)

    """
    # Lecture du graphe depuis un fichier JSON au format graphe_maille.json
    logger.info("Lecture du graphe maille depuis un fichier JSON.")
    graphe_filename = os.path.join(folderpath, "graphe_maille.json")
    if Path(graphe_filename).exists():
        with open(graphe_filename, "r") as f:
            graphe_data = json.load(f)
        # Conversion en DataFrame
        nodes = pd.DataFrame(graphe_data["nodes"])
        edges = pd.DataFrame(graphe_data["edges"])
        #arcs = pd.DataFrame(graphe_data["arcs"]) # arcs are directed, whereas edges are not
        # Conversion explicite en GeoDataFrame pour l'affichage
        nodes["geometry"] = nodes["geometry"].apply(lambda coords: Point(coords))
        nodes = gpd.GeoDataFrame(nodes, geometry="geometry")
        edges["geometry"] = edges["geometry"].apply(lambda coords: LineString(coords))
        edges = gpd.GeoDataFrame(edges, geometry="geometry")
        logger.info("Graphe maille chargé et converti en GeoDataFrame.")
    else:
        raise FileNotFoundError(f"Fichier <graphe_maille.json> introuvable dans le dossier {folderpath}. \n veuillez run 'generate_instances.py' ")
    return nodes, edges

def from_geodf_to_undirected_networkx(nodes: gpd.GeoDataFrame,
                                      edges: gpd.GeoDataFrame) -> nx.Graph :
    """
    Convert nodes and edges GeoDataFrames to a networkx graph.

    Parameters
    ----------
    nodes : geopandas.GeoDataFrame
        the nodes of the graph
    edges : geopandas.GeoDataFrame
        the edges of the graph (undirected)

    Returns
    -------
    G : nx.Graph
        networkx undirected Graph

    """
    G = nx.Graph()
    for idx, row in nodes.iterrows():
        point = row['geometry']
        G.add_node(idx, pos=(point.x, point.y),
                   label=row['type'],
                   conso_hiver = - row['demande_hivers'], 
                   conso_ete = row['demande']
                   )
    for line in edges.itertuples():
        start_idx = line.start_node_idx
        end_idx = line.end_node_idx
        G.add_edge(start_idx, end_idx, 
                    weight=line.length_km)
    check_reachable_nodes(G)
    return G

def check_reachable_nodes(G: nx.Graph):
    """
    Check if all nodes have connected edges, end the program otherwise.

    Parameters
    ----------
    G : nx.Graph

    Returns
    -------
    None.

    """
    unreachable = []
    for n, nbrsdict in G.adjacency():
        if len(nbrsdict)==0:
            unreachable += [n]
    print(unreachable)
    if len(unreachable) > 0:
        logger.error(f'Some nodes have no connected arcs: {unreachable}')
        exit(1)
    # print(len(G))
    # G.remove_nodes_from(unreachable)
    # print(len(G))
    # for n in unreachable:
    #     nodes.drop(n, inplace=True)
    # print(len(nodes))

def draw_graph_from_dataframes(nodes: gpd.GeoDataFrame,
                               edges: gpd.GeoDataFrame,
                               cells: dict=None,
                               edge_cells: dict=None,
                               nodes_close_to_source: dict=None,
                               used_edges_by_cells: dict=None,
                               principal_edges: list=None,
                               title: str=None
                               ) -> matplotlib.axes._axes.Axes:
    """
    Plot graph from GeoDataFrames nodes and edges, and draw nodes/edges
    in the color associated with their cell.

    Parameters
    ----------
    nodes : gpd.GeoDataFrame
        nodes of the graph (and data associated with nodes)
    edges : gpd.GeoDataFrame
        eddges of the graph (and data associated with edgess)
    cells : dict, optional
        Nodes cells in the form {cell_id: list of nodes in the cell}. The default is None.
    edge_cells : TYPE, optional
        Edges cells in the form {cell_id: list of edges in the cell}.
        Note that initial edges may be split in two if their starting and ending points
        are in different celles. In this case, a point is inserted in the initial segment,
        and two edges are created and associated with their respective cells.
        This is done in voronoi_cells_edges
        The default is None.

    Returns
    -------
    ax : matplotlib.axes._axes.Axes
        DESCRIPTION.

    """
    f, ax = plt.subplots()
    
    colors = ['b','hotpink','orange','mediumseagreen'] # /!\ only 4 colors ! TODO: improve
    
    if cells:
        inv_cells = { n : k for k in cells for n in cells[k]}
        nodes['cell'] = pd.Series(inv_cells)
        cell_color = { k : colors[k] if k != 'unreachable' else 'r' for k in cells }
        #print(cell_color)
        if not edge_cells:
            edges[edges['nb_lignes_existantes'] == 0].plot(ax=ax, color="grey", linewidth=0.2)
            #edges[edges['nb_lignes_existantes'] != 0].plot(ax=ax, color="red", linewidth=0.5)
        Plot_data = namedtuple('Plot_data', ['marker', 'markersize'])
        node_types = {'Source': Plot_data('s', 30),
                     'Transit': Plot_data('.', 5),
                     'Conso': Plot_data('o', 12),
                     'Prod': Plot_data('s', 12)}
        for k in cells:
            nodes_k = nodes[nodes['cell'] == k]
            color_k = cell_color[k]            
            for nt, nt_plot in node_types.items():
                nodes_nt = nodes_k[nodes_k['type'] == nt]
                if len(nodes_nt):
                    nodes_nt.plot(ax=ax, color=color_k, marker=nt_plot.marker, markersize=nt_plot.markersize)
            if nodes_close_to_source:
                criterion = nodes['id_noeud'].map(lambda i: i in nodes_close_to_source[k])
                nodes[criterion].plot(ax=ax, color='k', marker='x', markersize=10)

    if edge_cells:
        inv_cells = { edge : k for k in edge_cells for edge in edge_cells[k]}
        edges['node1'] = edges[['start_node_idx','end_node_idx']].min(axis=1)
        edges['node2'] = edges[['start_node_idx','end_node_idx']].max(axis=1)
        edges.set_index(['node1', 'node2'], inplace=True)
        edges['cell']=pd.Series(inv_cells) #.loc[edges['nodes']]
        if used_edges_by_cells:
            inv_used_edges = {edge : 1 for k in used_edges_by_cells for edge in used_edges_by_cells[k]}
            s = pd.Series(inv_used_edges)
            s.name='used'
            edges = pd.concat([edges, s], axis=1)
        cell_color = { k : colors[k] if k != 'split' else 'r' for k in edge_cells }
        for k in edge_cells:
            if k != 'split':
                edges_k = edges[edges['cell'] == k]
                if used_edges_by_cells:
                    edges_k[edges_k['used'] == 1].plot(ax=ax, color=cell_color[k], linewidth=1)
                    edges_k[edges_k['used'] != 1].plot(ax=ax, color=cell_color[k], linewidth=0.5, linestyle = "dotted")
                else :
                    edges_k.plot(ax=ax, color=cell_color[k], linewidth=0.2)
        nodes[nodes['type'] == 'dummy'].plot(ax=ax, color='k', marker=".", markersize=1)
    
    if principal_edges:
        edges['node1'] = edges[['start_node_idx','end_node_idx']].min(axis=1)
        edges['node2'] = edges[['start_node_idx','end_node_idx']].max(axis=1)
        edges.set_index(['node1', 'node2'], inplace=True)
        edges.loc[principal_edges].plot(ax=ax, color='black', linewidth=1, linestyle = "dotted")
    
    # Set margins for the axes so that nodes aren't clipped
    ax.plot()
    ax = plt.gca()
    ax.margins(0.20)
    plt.axis("off")
    if title:
        plt.title(title, fontsize=6)
    plt.show()
    
    return ax

# copied from networkx !
# only change : call multi_source_dijkstra_SB instead of multi_source_dijkstra_path
def voronoi_cells_SB(G: nx.Graph, center_nodes, weight="weight"):
    """Returns the Voronoi cells centered at `center_nodes` with respect
    to the shortest-path distance metric.

    If $C$ is a set of nodes in the graph and $c$ is an element of $C$,
    the *Voronoi cell* centered at a node $c$ is the set of all nodes
    $v$ that are closer to $c$ than to any other center node in $C$ with
    respect to the shortest-path distance metric. [1]_

    For directed graphs, this will compute the "outward" Voronoi cells,
    as defined in [1]_, in which distance is measured from the center
    nodes to the target node. For the "inward" Voronoi cells, use the
    :meth:`DiGraph.reverse` method to reverse the orientation of the
    edges before invoking this function on the directed graph.

    Parameters
    ----------
    G : NetworkX graph

    center_nodes : set
        A nonempty set of nodes in the graph `G` that represent the
        center of the Voronoi cells.

    weight : string or function
        The edge attribute (or an arbitrary function) representing the
        weight of an edge. This keyword argument is as described in the
        documentation for :func:`~networkx.multi_source_dijkstra_path`,
        for example.

    Returns
    -------
    dictionary
        A mapping from center node to set of all nodes in the graph
        closer to that center node than to any other center node. The
        keys of the dictionary are the element of `center_nodes`, and
        the values of the dictionary form a partition of the nodes of
        `G`.

    Examples
    --------
    To get only the partition of the graph induced by the Voronoi cells,
    take the collection of all values in the returned dictionary::

        >>> G = nx.path_graph(6)
        >>> center_nodes = {0, 3}
        >>> cells = nx.voronoi_cells(G, center_nodes)
        >>> partition = set(map(frozenset, cells.values()))
        >>> sorted(map(sorted, partition))
        [[0, 1], [2, 3, 4, 5]]

    Raises
    ------
    ValueError
        If `center_nodes` is empty.

    References
    ----------
    .. [1] Erwig, Martin. (2000),"The graph Voronoi diagram with applications."
        *Networks*, 36: 156--163.
        https://doi.org/10.1002/1097-0037(200010)36:3<156::AID-NET2>3.0.CO;2-L

    """
    # Determine the shortest paths from any one of the center nodes to
    # every node in the graph.
    #
    # This raises `ValueError` if `center_nodes` is an empty set.
    length, paths = multi_source_dijkstra_SB(G, center_nodes, weight=weight)
    # Determine the center node from which the shortest path originates.
    nearest = {v: p[0] for v, p in paths.items()}
    # Get the mapping from center node to all nodes closer to it than to
    # any other center node.
    cells = nx.utils.groups(nearest)
    # We collect all unreachable nodes under a special key, if there are any.
    unreachable = set(G) - set(nearest)
    if unreachable:
        cells["unreachable"] = unreachable
    return length, cells, paths

# copied from networkx, only change : call to _dijkstra_multisource_SB instead
# of _dijkstra_multisource, to use the capacity of sources
def multi_source_dijkstra_SB(G, sources, target=None, cutoff=None, weight="weight"):
    """Find shortest weighted paths and lengths from a given set of
    source nodes.

    Uses Dijkstra's algorithm to compute the shortest paths and lengths
    between one of the source nodes and the given `target`, or all other
    reachable nodes if not specified, for a weighted graph.

    Parameters
    ----------
    G : NetworkX graph

    sources : non-empty set of nodes
        Starting nodes for paths. If this is just a set containing a
        single node, then all paths computed by this function will start
        from that node. If there are two or more nodes in the set, the
        computed paths may begin from any one of the start nodes.

    target : node label, optional
        Ending node for path

    cutoff : integer or float, optional
        Length (sum of edge weights) at which the search is stopped.
        If cutoff is provided, only return paths with summed weight <= cutoff.

    weight : string or function
        If this is a string, then edge weights will be accessed via the
        edge attribute with this key (that is, the weight of the edge
        joining `u` to `v` will be ``G.edges[u, v][weight]``). If no
        such edge attribute exists, the weight of the edge is assumed to
        be one.

        If this is a function, the weight of an edge is the value
        returned by the function. The function must accept exactly three
        positional arguments: the two endpoints of an edge and the
        dictionary of edge attributes for that edge. The function must
        return a number or None to indicate a hidden edge.

    Returns
    -------
    distance, path : pair of dictionaries, or numeric and list
        If target is None, returns a tuple of two dictionaries keyed by node.
        The first dictionary stores distance from one of the source nodes.
        The second stores the path from one of the sources to that node.
        If target is not None, returns a tuple of (distance, path) where
        distance is the distance from source to target and path is a list
        representing the path from source to target.

    Examples
    --------
    >>> G = nx.path_graph(5)
    >>> length, path = nx.multi_source_dijkstra(G, {0, 4})
    >>> for node in [0, 1, 2, 3, 4]:
    ...     print(f"{node}: {length[node]}")
    0: 0
    1: 1
    2: 2
    3: 1
    4: 0
    >>> path[1]
    [0, 1]
    >>> path[3]
    [4, 3]

    >>> length, path = nx.multi_source_dijkstra(G, {0, 4}, 1)
    >>> length
    1
    >>> path
    [0, 1]

    Notes
    -----
    Edge weight attributes must be numerical.
    Distances are calculated as sums of weighted edges traversed.

    The weight function can be used to hide edges by returning None.
    So ``weight = lambda u, v, d: 1 if d['color']=="red" else None``
    will find the shortest red path.

    Based on the Python cookbook recipe (119466) at
    https://code.activestate.com/recipes/119466/

    This algorithm is not guaranteed to work if edge weights
    are negative or are floating point numbers
    (overflows and roundoff errors can cause problems).

    Raises
    ------
    ValueError
        If `sources` is empty.
    NodeNotFound
        If any of `sources` is not in `G`.

    See Also
    --------
    multi_source_dijkstra_path
    multi_source_dijkstra_path_length

    """
    if not sources:
        raise ValueError("sources must not be empty")
    for s in sources:
        if s not in G:
            raise nx.NodeNotFound(f"Node {s} not found in graph")
    if target in sources:
        return (0, [target])
    weight = nx.weighted._weight_function(G, weight)
    paths = {source: [source] for source in sources}  # dictionary of paths
    dist = _dijkstra_multisource_SB(
        G, sources, weight, paths=paths, cutoff=cutoff, target=target
    )
    if target is None:
        return (dist, paths)
    try:
        return (dist[target], paths[target])
    except KeyError as err:
        raise nx.NetworkXNoPath(f"No path to {target}.") from err

# copied from networkx !
# adapted to handle source capacities (lines marked with "# SB")
def _dijkstra_multisource_SB(G, sources, weight, pred=None, paths=None, cutoff=None, target=None):
    """Uses Dijkstra's algorithm to find shortest weighted paths.

    Parameters
    ----------
    G : NetworkX graph

    sources : non-empty iterable of nodes
        Starting nodes for paths. If this is just an iterable containing
        a single node, then all paths computed by this function will
        start from that node. If there are two or more nodes in this
        iterable, the computed paths may begin from any one of the start
        nodes.

    weight: function
        Function with (u, v, data) input that returns that edge's weight
        or None to indicate a hidden edge

    pred: dict of lists, optional(default=None)
        dict to store a list of predecessors keyed by that node
        If None, predecessors are not stored.

    paths: dict, optional (default=None)
        dict to store the path list from source to each node, keyed by node.
        If None, paths are not stored.

    target : node label, optional
        Ending node for path. Search is halted when target is found.

    cutoff : integer or float, optional
        Length (sum of edge weights) at which the search is stopped.
        If cutoff is provided, only return paths with summed weight <= cutoff.

    Returns
    -------
    distance : dictionary
        A mapping from node to shortest distance to that node from one
        of the source nodes.

    Raises
    ------
    NodeNotFound
        If any of `sources` is not in `G`.

    Notes
    -----
    The optional predecessor and path dictionaries can be accessed by
    the caller through the original pred and paths objects passed
    as arguments. No need to explicitly return pred or paths.

    """
    G_succ = G._adj  # For speed-up (and works for both directed and undirected graphs)

    push = heappush
    pop = heappop
    dist = {}  # dictionary of final distances
    seen = {}
    # fringe is heapq with 3-tuples (distance,c,node)
    # use the count c to avoid comparing nodes (may not be able to)
    c = count()
    fringe = []
    closest_source = {} # SB
    for source in sources:
        seen[source] = 0
        closest_source[source] = source # SB
        push(fringe, (0, next(c), source))
    while fringe:
        (d, _, v) = pop(fringe)
        if v in dist:
            continue  # already searched this node.
        dist[v] = d
        if v == target:
            break
        closest_s = closest_source[v]   # SB
        capacity = G.nodes[closest_s]['conso_ete']  # SB
        for u, e in G_succ[v].items():
            cost = weight(v, u, e)/capacity     # SB
            if cost is None:
                continue
            vu_dist = dist[v] + cost
            if cutoff is not None:
                if vu_dist > cutoff:
                    continue
            if u in dist:
                u_dist = dist[u]
                if vu_dist < u_dist:
                    raise ValueError("Contradictory paths found:", "negative weights?")
                elif pred is not None and vu_dist == u_dist:
                    pred[u].append(v)
            elif u not in seen or vu_dist < seen[u]:
                seen[u] = vu_dist
                closest_source[u] = closest_s
                push(fringe, (vu_dist, next(c), u))
                if paths is not None:
                    paths[u] = paths[v] + [u]
                if pred is not None:
                    pred[u] = [v]
            elif vu_dist == seen[u]:
                if pred is not None:
                    pred[u].append(v)

    # The optional predecessor and path dictionaries can be accessed
    # by the caller via the pred and paths objects passed as arguments.
    return dist

def voronoi_cells_edges(G: nx.Graph, cells: dict, dist: dict)-> Tuple[dict, list]:
    """
    Put edges in cells, and divide edges whose ends are in different cells, in
    the point that is at equal distance of both sources, with the shortest path
    distance divided by the capacity of the source (for now, maybe a better distance later)
    

    Parameters
    ----------
    G : nx.Graph
        DESCRIPTION.
    cells : dict
        DESCRIPTION.
    dist : dict
        DESCRIPTION.

    Returns
    -------
    Tuple[dict, list]
        a dict that contains the cells and the edges associated with each cell
        and a list of added edges for edges that have been divided

    """
    edge_cells = { k : [] for k in cells}
    edge_cells['split'] = []
    # assign nodes to Voronoi cells
    node_cell = { n : k for k in cells for n in cells[k]}
    
    new_nodes_index = len(G) # /!\ will not work if nodes have been deleted from G
    nodes_to_add = []
    dist_limit_nodes = {}
    
    # process each edge
    for u, v, d in G.edges.data('weight'):
        if node_cell[u] == node_cell[v]:
            # the entire edge is in te cell
            edge_cells[node_cell[u]] += [(u,v)]
        else:
            capacity_u = G.nodes[node_cell[u]]['conso_ete']
            capacity_v = G.nodes[node_cell[v]]['conso_ete']
            # divide the edge : find point q on segment [u,v], s.t. dist(u) + d(u,q) = dist(v) + d(v,q)
            # q is on u,v so duq + dqv = d
            # TODO : double check formula
            pos_q_in_uv = ((dist[v] - dist[u]) * capacity_v * capacity_u / d + capacity_u) / (capacity_v + capacity_u) if d != 0 else 1
            pos_u = G.nodes[u]['pos']
            pos_v = G.nodes[v]['pos']
            # print(pos_u, pos_v, pos_q_in_uv)
            if pos_q_in_uv < 0 or pos_q_in_uv > 1 :
                exit(0) # TODO mieux gerer ca
            pos_q = np.array(pos_u) + pos_q_in_uv * (np.array(pos_v) - np.array(pos_u))
            # print((u, node_cell[u], pos_u, dist[u]),
            #        (v, node_cell[v], pos_v, dist[v]), pos_q_in_uv, pos_q)
            # insert new node at q
            q = new_nodes_index
            new_nodes_index += 1
            nodes_to_add += [(q, pos_q, u, v, d * pos_q_in_uv, d * (1 - pos_q_in_uv))]
            # assign (u,q) to node_cell[u] and (q,v) to node_cell[v]
            edge_cells[node_cell[u]] += [(u,q)]
            edge_cells[node_cell[v]] += [(v,q)]
            edge_cells['split'] += [(u,v)]
            dist_limit_nodes[q] = {node_cell[u]: dist[u] + d * pos_q_in_uv / capacity_u,
                                   node_cell[v]: dist[v] + d * (1 - pos_q_in_uv) / capacity_v} # should be twice the same distance # OK
    for (q, pos_q, u, v, l1, l2) in nodes_to_add:
        G.add_node(q, label='dummy', pos=pos_q)
        G.add_edge(u, q, label='dummy', weight=l1)
        G.add_edge(v, q, label='dummy', weight=l2)
    return edge_cells, nodes_to_add, dist_limit_nodes

def plot_with_plotly(G):
    import plotly.graph_objects as go
    graph_edges = [(u,v) for u,v,d in G.edges.data('label') if d != 'dummy']
    coord = [(G.nodes[u]['pos'], G.nodes[v]['pos']) for u,v in graph_edges]
    edge_x = []
    edge_y = []
    # nodes in the middle for hover
    mnode_x = []
    mnode_y = []
    edge_label = []
    for p,q in coord:
        x_0, y_0 = p
        x_1, y_1 = q
        edge_x.extend([x_0, x_1, None])
        edge_y.extend([y_0, y_1, None])
        # nodes for hover
        mnode_x.extend([(x_0 + x_1)/2])
        mnode_y.extend([(y_0 + y_1)/2])
        mnode_x.extend([(3*x_0 + x_1)/4])
        mnode_y.extend([(3*y_0 + y_1)/4])
        mnode_x.extend([(x_0 + 3*x_1)/4])
        mnode_y.extend([(y_0 + 3*y_1)/4])
    # 3 "middle" points for each edge
    for u,v in graph_edges:
        edge_label += [f'{u}-{v}']*3
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        mode='lines')
    
    mnode_trace = go.Scatter(x = mnode_x, y = mnode_y, mode = "markers", showlegend = False,
                         #hovertemplate = "Edge %{hovertext}<extra></extra>",
                         hoverinfo='text',
                         hovertext = edge_label, marker = go.Marker(opacity = 0))
    
    graph_nodes = [n for (n,d) in G.nodes.data('label') if d != 'dummy']
    node_x = [G.nodes[n]['pos'][0] for n in graph_nodes]
    node_y = [G.nodes[n]['pos'][1] for n in graph_nodes]
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers', #"markers+text
        hoverinfo='text',
        #hovertext=node_labels,
        marker=dict(
            # showscale=True,
            # colorscale options
            #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
            #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
            #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
            # colorscale='YlGnBu',
            # reversescale=True,
            color='orange',
            size=10, #node_sizes
            # colorbar=dict(
            #     thickness=15,
            #     title=dict(
            #       text='Node Connections',
            #       side='right'
            #     ),
            #     xanchor='left',
            # ),
            line_width=2), #line=dict(width=1, color="black")
        text=list(G.nodes()),
        textposition="bottom center",
        textfont=dict(size=12),
        )
    
    fig = go.Figure(layout=go.Layout(showlegend=False,
                                     hovermode='closest',
                                     xaxis=dict(showgrid=True, zeroline=False, showticklabels=True),
                                     yaxis=dict(showgrid=True, zeroline=False, showticklabels=True)))

    fig.add_trace(edge_trace)
    fig.add_trace(node_trace)
    fig.add_trace(mnode_trace)
    fig.show()
    
def edges_on_shortest_paths(cells: dict, paths:dict) -> dict:
    """
    For each cell, compute and store edges on shortest paths, 
    from nodes on shortest paths for each node to its closest source

    Parameters
    ----------
    cells : dict
        DESCRIPTION.
    paths : dict
        DESCRIPTION.

    Returns
    -------
    dict
        DESCRIPTION.

    """
    paths_by_cells = {k: [] for k in cells}
    for k,ns in cells.items():
        for n in ns:
            present = False
            update = False
            for i, p in enumerate(paths_by_cells[k]):
                if n in p:
                    present = True
                    break
                if p[-1] in paths[n]:
                    paths_by_cells[k][i] = paths[n]
                    update = True
                    break
            if present or update:
                continue
            paths_by_cells[k] += [paths[n]]
    used_edges_by_cells = {k: [] for k in cells}
    for k, ps in paths_by_cells.items():
        for p in ps:
            for i in range(len(p)-1):
                (u,v) = (p[i], p[i+1]) if p[i] < p[i+1] else (p[i+1], p[i])
                if (u,v) not in used_edges_by_cells[k]:
                    used_edges_by_cells[k] += [(u,v)]
    return used_edges_by_cells

if __name__ == "__main__":
    folderpath = "exemple_ete_4_PS"
    # generate_instances.main(folderpath, run_new_simulation=True)
    nodes, edges = lire_json_graphe(folderpath)

    # generate_instances.affiche_graphe_maille(nodes, edges).plot()
    # figurepath = os.path.join(folderpath, "Figures/")
    # plt.savefig(figurepath + "graphe_maille.pdf", bbox_inches="tight", dpi=100, format="pdf" )
    # plt.show()

    # networkx undirected graph
    G = from_geodf_to_undirected_networkx(nodes, edges)

    # plot tests
    # pos = {n: G.nodes[n]['pos'] for n in G}
    # options = {
    #     "font_size": 8,
    #     "node_size": 150,
    #     "node_color": "white",
    #     "edgecolors": "black",
    #     "linewidths": 1,
    #     "width": 1,
    # }
    # nx.draw_networkx(G, pos, **options)
    # plot_with_plotly(G)

    center_nodes = [n for n in G if G.nodes('label')[n]=='Source']
    
    ###########################################################################
    # COMPUTE VORONOI CELLS
    ###########################################################################
    
    # networkx voronoi cells with weight = edge length
    # cells = nx.voronoi_cells(G, center_nodes, weight='weight' )#, weight=weight_function)
    # draw_graph_from_dataframes(nodes, edges, cells)
    
    # voronoi cells with weight = edge length / capacity of closest source
    dist, cells, paths = voronoi_cells_SB(G, center_nodes)
    # draw_graph_from_dataframes(nodes, edges, cells)
    
    # compute used edges (edges on the shortest paths) 
    # ideally these edges would be the same as the ones associated with lines to build
    # but here many constraints from the original problem are missing (feeders, power limit, length limit, "power*length" limit...)
    used_edges_by_cells = edges_on_shortest_paths(cells, paths)
    ###########################################################################
    # BALANCE CHARGES BETWEEN CELLS ?
    ###########################################################################
    # capa_by_source = np.array([G.nodes('conso_ete')[n] for n in center_nodes])
    # for i in range(150):
    #     print(i, [G.nodes('conso_ete')[n] for n in center_nodes]) #if G.nodes('label')[n]=='Source'])
    #     change_capa_and_iterate(G, cells, capa_by_source)
    #     length, cells, paths = voronoi_cells_SB(G, center_nodes)
    #     draw_graph_from_dataframes(nodes, edges, cells)
    
    ###########################################################################
    # COMPUTE LIMITS OF CELLS ON SHARED EDGES
    ###########################################################################
    edge_cells, nodes_to_add, dist_limit_nodes = voronoi_cells_edges(G, cells, dist)
    edges1 = edges.copy(deep=True)
    nodes1 = nodes.copy(deep=True)
    for (q, pos_q, u, v, l1, l2) in nodes_to_add:
        #print(u,v,q)
        pos_u = G.nodes[u]['pos']
        pos_v = G.nodes[v]['pos']
        edges1.loc[len(edges1)] = pd.Series({'edge_idx': len(edges), 'length_km':100,
                                           'start_node_idx': u,  'end_node_idx': q,
                                           'nb_lignes_existantes': None , 'geometry': LineString([pos_u, pos_q])})
        edges1.loc[len(edges1)] = pd.Series({'edge_idx': len(edges), 'length_km':100,
                                           'start_node_idx': v,  'end_node_idx': q,
                                           'nb_lignes_existantes': None , 'geometry': LineString([pos_v, pos_q])})
        nodes1.loc[len(nodes1)] = pd.Series({'type': 'dummy', 'geometry': Point(pos_q)})
    draw_graph_from_dataframes(nodes1, edges1, cells, edge_cells, used_edges_by_cells=used_edges_by_cells)
    
    ###########################################################################
    # test : notion de "principal"
    # assume principal are edges and nodes on shortest paths between sources for now
    # (select any shortest path if several)
    ###########################################################################
    G = from_geodf_to_undirected_networkx(nodes, edges)
    principals = []
    for (source1, source2) in list(combinations(center_nodes,2)):
        principals += [nx.shortest_path(G, source1, source2, weight='weight')]
    # set these edges weight (length) to 0
    for ppal in principals:
        for i in range(len(ppal)-1):
            u = ppal[i]
            v = ppal[i+1]
            G.edges[u,v]['weight'] = 0.2
    
    dist, cells, paths = voronoi_cells_SB(G, center_nodes)
    edge_cells, nodes_to_add, dist_limit_nodes = voronoi_cells_edges(G, cells, dist)
    used_edges_by_cells = edges_on_shortest_paths(cells, paths)
    edges1 = edges.copy(deep=True)
    nodes1 = nodes.copy(deep=True)
    for (q, pos_q, u, v, l1, l2) in nodes_to_add:
        #print(u,v,q)
        pos_u = G.nodes[u]['pos']
        pos_v = G.nodes[v]['pos']
        edges1.loc[len(edges1)] = pd.Series({'edge_idx': len(edges), 'length_km':100,
                                           'start_node_idx': u,  'end_node_idx': q,
                                           'nb_lignes_existantes': None , 'geometry': LineString([pos_u, pos_q])})
        edges1.loc[len(edges1)] = pd.Series({'edge_idx': len(edges), 'length_km':100,
                                           'start_node_idx': v,  'end_node_idx': q,
                                           'nb_lignes_existantes': None , 'geometry': LineString([pos_v, pos_q])})
        nodes1.loc[len(nodes1)] = pd.Series({'type': 'dummy', 'geometry': Point(pos_q)})
    principal_edges = [sorted((ppal[i],ppal[i+1])) for ppal in principals for i in range(len(ppal)-1)]
    draw_graph_from_dataframes(nodes1, edges1, cells, edge_cells, used_edges_by_cells=used_edges_by_cells, principal_edges=principal_edges)
    
    edges1 = edges.copy(deep=True)
    nodes1 = nodes.copy(deep=True)
    # test minimum spanning trees 
    used_edges_by_cells = {k: [] for k in cells}
    for k in center_nodes:
        G_k = G.subgraph(cells[k])
        for (p1, p2) in nx.minimum_spanning_tree(G_k).edges:
            (u,v) = (p1, p2) if p1 < p2 else (p2, p1)
            used_edges_by_cells[k] += [(u,v)]
    for (q, pos_q, u, v, l1, l2) in nodes_to_add:
        #print(u,v,q)
        pos_u = G.nodes[u]['pos']
        pos_v = G.nodes[v]['pos']
        edges1.loc[len(edges1)] = pd.Series({'edge_idx': len(edges), 'length_km':100,
                                           'start_node_idx': u,  'end_node_idx': q,
                                           'nb_lignes_existantes': None , 'geometry': LineString([pos_u, pos_q])})
        edges1.loc[len(edges1)] = pd.Series({'edge_idx': len(edges), 'length_km':100,
                                           'start_node_idx': v,  'end_node_idx': q,
                                           'nb_lignes_existantes': None , 'geometry': LineString([pos_v, pos_q])})
        nodes1.loc[len(nodes1)] = pd.Series({'type': 'dummy', 'geometry': Point(pos_q)})
    principal_edges = [sorted((ppal[i],ppal[i+1])) for ppal in principals for i in range(len(ppal)-1)]
    draw_graph_from_dataframes(nodes1, edges1, cells, edge_cells, used_edges_by_cells=used_edges_by_cells, principal_edges=principal_edges)